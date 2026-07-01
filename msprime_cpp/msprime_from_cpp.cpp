/*
 * Exemple MINIMAL de l'option B1 : appeler msprime depuis C++, sans
 * pybind11, en utilisant directement l'API Python C (plus bas niveau
 * mais sans dépendance externe supplémentaire -- Python.h est livré
 * avec Python lui-même).
 *
 * Compile avec (adapter les chemins à ton env conda) :
 *   g++ -o msprime_from_cpp msprime_from_cpp.cpp \
 *       $(python3-config --includes --ldflags) \
 *       -lpython3.11 \
 *       -O2
 *
 * Ce que ça fait :
 *   1. Initialise un interpréteur Python embarqué dans le process C++
 *   2. Importe msprime
 *   3. Simule un arbre de coalescence (2 populations, 5 ind. chacune)
 *   4. Récupère le temps de coalescence racine (TMRCA) en C++
 *   5. Finalise l'interpréteur et affiche le résultat
 *
 * C'est le mécanisme exact qu'on utiliserait dans dosimulpart() pour
 * remplacer la construction manuelle du GeneTreeC par un appel à msprime.
 *
 * AVANTAGES vs notre architecture actuelle (subprocess + fichier .snp) :
 *   - Pas d'I/O disque : les génotypes restent en mémoire (RAM)
 *   - Pas de démarrage de processus à chaque particule
 *   - Un seul interpréteur Python initialisé une fois pour tout le run
 *
 * INCONVÉNIENTS / COMPLEXITÉS :
 *   - Python.h doit être disponible au compile-time
 *   - Le GIL (Global Interpreter Lock) de Python doit être géré si on
 *     parallélise avec OpenMP -- chaque thread doit acquérir le GIL
 *     avant d'appeler Python, ce qui peut sérialiser la parallélisation
 *   - Le binaire résultant dépend d'une version précise de Python
 *     (fragile sur les machines où l'environnement conda change)
 *   - Gestion de la mémoire Python (Py_INCREF/Py_DECREF) est fastidieuse
 *     et source de leaks si mal gérée
 */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <vector>
#include <stdexcept>

/*
 * Initialise l'interpréteur Python et importe msprime.
 * À appeler UNE SEULE FOIS au démarrage du programme (équivalent de
 * Py_Initialize() dans main() du vrai DIYABC).
 */
void init_python() {
    Py_Initialize();
    if (!Py_IsInitialized()) {
        throw std::runtime_error("Impossible d'initialiser l'interpréteur Python");
    }
    std::cout << "Interpréteur Python initialisé." << std::endl;
}

/*
 * Simule un arbre de coalescence pour un seul locus et retourne le TMRCA
 * (time to most recent common ancestor) -- l'équivalent C++ d'un appel à
 * msprime.sim_ancestry().
 *
 * Dans le vrai DIYABC, cette fonction remplacerait la construction du
 * GeneTreeC dans dosimulpart() -- elle retournerait les génotypes simulés
 * (0/1 par lignée) plutôt que juste le TMRCA, mais le principe est identique.
 */
double simulate_one_locus_tmrca(
    int n_pop1,
    int n_pop2,
    double Ne,
    double t_split,
    int seed
) {
    /* Script Python minimal, généré dynamiquement depuis le C++ -- dans
     * le vrai DIYABC, ce script serait construit à partir des paramètres
     * tirés par drawparams() et des événements du scénario actif. */
    std::string python_code = R"(
import msprime

demography = msprime.Demography()
demography.add_population(name="pop1", initial_size=)" + std::to_string((int)Ne) + R"()
demography.add_population(name="pop2", initial_size=)" + std::to_string((int)Ne) + R"()
demography.add_population(name="ANC",  initial_size=)" + std::to_string((int)Ne) + R"()
demography.add_population_split(
    time=)" + std::to_string((int)t_split) + R"(,
    derived=["pop1", "pop2"],
    ancestral="ANC"
)

ts = msprime.sim_ancestry(
    samples={"pop1": )" + std::to_string(n_pop1) + R"(, "pop2": )" + std::to_string(n_pop2) + R"(},
    demography=demography,
    sequence_length=1,
    random_seed=)" + std::to_string(seed) + R"(,
    ploidy=2,
)

tree = ts.first()
result_tmrca = tree.time(tree.root)
)";

    /* Exécute le script dans l'interpréteur embarqué */
    PyObject* main_module = PyImport_AddModule("__main__");
    PyObject* global_dict = PyModule_GetDict(main_module);

    PyObject* result = PyRun_String(
        python_code.c_str(),
        Py_file_input,
        global_dict,
        global_dict
    );

    if (!result) {
        PyErr_Print();
        throw std::runtime_error("Erreur lors de l'exécution du script msprime");
    }
    Py_DECREF(result);

    /* Récupère la variable result_tmrca depuis le namespace Python */
    PyObject* tmrca_obj = PyDict_GetItemString(global_dict, "result_tmrca");
    if (!tmrca_obj) {
        throw std::runtime_error("Variable result_tmrca non trouvée dans le namespace Python");
    }

    double tmrca = PyFloat_AsDouble(tmrca_obj);
    return tmrca;
}

void finalize_python() {
    Py_Finalize();
    std::cout << "Interpréteur Python finalisé." << std::endl;
}

int main() {
    try {
        init_python();

        /* Paramètres de la simulation -- dans le vrai DIYABC, ces valeurs
         * viendraient de drawparams() et du scénario tiré par drawscenario() */
        int    n_pop1   = 5;       /* individus dans pop1 */
        int    n_pop2   = 5;       /* individus dans pop2 */
        double Ne       = 10000.0; /* taille effective des populations */
        double t_split  = 2000.0;  /* temps de séparation (générations) */

        std::cout << "\nSimulation de 5 loci indépendants :" << std::endl;
        std::cout << "  Ne=" << Ne << "  t_split=" << t_split << std::endl;
        std::cout << "  n_pop1=" << n_pop1 << "  n_pop2=" << n_pop2 << std::endl;

        for (int locus = 0; locus < 5; locus++) {
            double tmrca = simulate_one_locus_tmrca(
                n_pop1, n_pop2, Ne, t_split,
                locus + 1  /* seed = locus + 1, msprime rejette seed=0 */
            );
            std::cout << "  Locus " << locus
                      << " : TMRCA = " << tmrca << " générations" << std::endl;
        }

        finalize_python();

    } catch (const std::exception& e) {
        std::cerr << "ERREUR : " << e.what() << std::endl;
        return 1;
    }
    return 0;
}