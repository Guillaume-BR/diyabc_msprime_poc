"""
Vérifie que scenario_parser produit, sur le vrai header.txt du dataset
human, exactement les événements qu'on a décortiqués à la main avec le
mentor (voir notes/exploration.md) pour le scénario 1.
"""

import os
import shutil
import struct
from pathlib import Path

import msprime
import pytest

from bridge.ancestry_simulation import (
    build_samples_argument,
    simulate_independent_loci,
    simulate_snp_genotypes,
)
from bridge.demography_builder import (
    build_demography,
    evaluate_expression,
    extract_referenced_names,
    get_parameter_names_used_by_scenario,
)
from bridge.loci_parser import parse_loci_description
from bridge.observed_data import count_samples_per_population, population_index_to_name
from bridge.parameter_sampling import (
    _draw_one_value,
    draw_parameter_values,
    draw_scenario,
)
from bridge.pipeline import (
    build_random_demography_for_scenario_index,
    compute_summary_statistics,
    read_header_text,
    run_poc_for_directory,
)
from bridge.prior_parser import is_constant_prior, parse_priors
from bridge.reftable_loop import (
    run_reftable_simulation,
    simulate_reference_directory,
    write_reftable_bin,
    write_reftable_txt,
)
from bridge.scenario_parser import parse_header_scenarios
from bridge.scenario_types import (
    MergeEvent,
    OrderConstraint,
    Prior,
    SampleEvent,
    Scenario,
    SplitEvent,
    VarNeEvent,
)
from bridge.snp_writer import write_snp_file
from bridge.stats_group_parser import parse_requested_statistic_names

REFERENCE_DIR = Path(__file__).parent.parent / "reference" / "human"
GENERAL_BINARY_PATH = os.environ.get("DIYABC_GENERAL_PATH")


@pytest.fixture
def header_text() -> str:
    return (REFERENCE_DIR / "header.txt").read_text()


def test_unimplemented_scenarios_are_skipped_with_warning(header_text):
    """Les scénarios 2,3,5,6 utilisent 'split' (pas encore implémenté) :
    ils doivent être ignorés avec un avertissement, pas faire planter le
    parsing des autres scénarios."""
    # with pytest.warns(UserWarning, match="split"):
    scenarios = parse_header_scenarios(header_text)

    # Scénarios gérables aujourd'hui : 1 (merge/varNe) et 4 (idem,
    # numérotation t11..t44). Les 4 autres utilisent split -> ignorés.
    found_indices = {s.index for s in scenarios}
    assert found_indices == {1, 2, 3, 4, 5, 6}


def test_scenario1_metadata(header_text):
    scenarios = parse_header_scenarios(header_text)
    scenario1 = scenarios[0]
    assert scenario1.index == 1
    assert scenario1.weight == pytest.approx(0.16667)
    assert scenario1.initial_pop_size_exprs == ["N1", "N2", "N3", "N4"]


def test_scenario1_events(header_text):
    scenarios = parse_header_scenarios(header_text)
    scenario1 = scenarios[0]

    expected = [
        SampleEvent(time_expr="0", pop=1),
        SampleEvent(time_expr="0", pop=2),
        SampleEvent(time_expr="0", pop=3),
        SampleEvent(time_expr="0", pop=4),
        MergeEvent(time_expr="t1", ancestral_pop=2, derived_pop=1),
        VarNeEvent(time_expr="t2-d3", pop=3, new_size_expr="Nbn3"),
        VarNeEvent(time_expr="t2-d4", pop=4, new_size_expr="Nbn4"),
        MergeEvent(time_expr="t2", ancestral_pop=3, derived_pop=4),
        VarNeEvent(time_expr="t2", pop=3, new_size_expr="N34"),
        VarNeEvent(time_expr="t3-d34", pop=3, new_size_expr="Nbn34"),
        MergeEvent(time_expr="t3", ancestral_pop=2, derived_pop=3),
        VarNeEvent(time_expr="t4", pop=2, new_size_expr="Na"),
    ]

    assert scenario1.events == expected


def test_scenario4_events(header_text):
    scenarios = parse_header_scenarios(header_text)
    scenario4 = next(s for s in scenarios if s.index == 4)

    expected = [
        SampleEvent(time_expr="0", pop=1),
        SampleEvent(time_expr="0", pop=2),
        SampleEvent(time_expr="0", pop=3),
        SampleEvent(time_expr="0", pop=4),
        MergeEvent(time_expr="t11", ancestral_pop=2, derived_pop=1),
        VarNeEvent(time_expr="t22-d3", pop=3, new_size_expr="Nbn3"),
        MergeEvent(time_expr="t22", ancestral_pop=2, derived_pop=3),
        VarNeEvent(time_expr="t33-d4", pop=4, new_size_expr="Nbn4"),
        MergeEvent(time_expr="t33", ancestral_pop=2, derived_pop=4),
        VarNeEvent(time_expr="t44", pop=2, new_size_expr="Na"),
    ]

    assert scenario4.events == expected


def test_scenario2_events(header_text):
    """Vérifie qu'un événement 'split' (admixture) est correctement
    interprété : 't1 split 1 4 2 ra' -> pop 1 disparaît, chaque lignée
    part vers pop 4 avec probabilité 'ra', sinon vers pop 2 -- sémantique
    vérifiée contre history.cpp/particuleC.cpp (voir docstring de
    SplitEvent)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    assert scenario2.events[4] == SplitEvent(
        time_expr="t1",
        derived_pop=1,
        ancestral_pop1=4,
        ancestral_pop2=2,
        admixture_rate="ra",
    )


def test_build_demography_scenario2_admixture(header_text):
    """Vérifie que build_demography traduit un SplitEvent en événement
    msprime Admixture avec les bonnes populations et proportions
    (rate, 1-rate)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
        "ra": 0.3,
    }

    demography = build_demography(scenario2, values)

    admixtures = [
        e for e in demography.events if isinstance(e, msprime.demography.Admixture)
    ]
    assert len(admixtures) == 1
    admixture = admixtures[0]
    assert admixture.time == 10
    assert admixture.derived == "pop1"
    assert admixture.ancestral == ["pop4", "pop2"]
    assert admixture.proportions == pytest.approx([0.3, 0.7])


def test_get_parameter_names_used_by_scenario2(header_text):
    """Vérifie que le taux d'admixture 'ra' est bien inclus dans les
    paramètres référencés par un scénario qui contient un SplitEvent --
    sinon il serait exclu à tort des colonnes du reftable.bin pour ce
    scénario (même bug que celui déjà corrigé pour t11..t44)."""
    scenarios = parse_header_scenarios(header_text)
    scenario2 = next(s for s in scenarios if s.index == 2)

    used_names = get_parameter_names_used_by_scenario(scenario2)

    assert "ra" in used_names
    expected = {
        "N1",
        "N2",
        "N3",
        "N4",
        "t1",
        "ra",
        "t2",
        "d3",
        "Nbn3",
        "d4",
        "Nbn4",
        "N34",
        "t3",
        "d34",
        "Nbn34",
        "t4",
        "Na",
    }
    assert used_names == expected


def test_priors_and_constraints(header_text):
    priors, constraints = parse_priors(header_text)

    assert len(priors) == 21
    assert len(constraints) == 4

    priors_by_name = {p.name: p for p in priors}
    assert priors_by_name["N1"].category == "N"
    assert priors_by_name["N1"].law == "UN"
    assert priors_by_name["N1"].bounds == (1000.0, 100000.0, 0.0, 0.0)

    assert priors_by_name["t1"].category == "T"
    assert priors_by_name["t1"].bounds == (1.0, 30.0, 0.0, 0.0)

    assert OrderConstraint(param1="t4", operator=">", param2="t3") in constraints
    assert OrderConstraint(param1="t3", operator=">", param2="t2") in constraints
    assert OrderConstraint(param1="t44", operator=">", param2="t33") in constraints
    assert OrderConstraint(param1="t44", operator=">", param2="t22") in constraints


def test_draw_parameter_values(header_text):
    """Vérifie que draw_parameter_values tire bien une valeur pour chaque
    prior, et que le tirage retourné respecte toutes les contraintes
    d'ordre (t4>t3, t3>t2, t44>t33, t44>t22)."""
    priors, constraints = parse_priors(header_text)

    seed = 42
    values = draw_parameter_values(priors, constraints, seed)

    # Toutes les valeurs ont bien été tirées
    assert set(values.keys()) == {p.name for p in priors}

    # Toutes les contraintes sont respectées par ce tirage
    for constraint in constraints:
        assert constraint.is_satisfied(values)


def test_draw_parameter_values_reproducible(header_text):
    """Même graine -> même tirage (déterminisme attendu pour la
    reproductibilité scientifique)."""
    priors, constraints = parse_priors(header_text)

    values1 = draw_parameter_values(priors, constraints, seed=123)
    values2 = draw_parameter_values(priors, constraints, seed=123)

    assert values1 == values2


def test_draw_parameter_values_rounds_N_and_T_not_A(header_text):
    """DIYABC arrondit à l'entier le plus proche les priors de catégorie
    N (taille) et T (temps) juste après le tirage, mais laisse le taux
    d'admixture (catégorie A, ex: 'ra') continu -- vérifié contre
    particuleC.cpp ("if (category<2) value = floor(0.5+value)", avec
    category 0=N, 1=T, 2=A dans header.cpp)."""
    priors, constraints = parse_priors(header_text)
    priors_by_name = {p.name: p for p in priors}

    non_integer_ra_seen = False
    for seed in range(1, 50):
        values = draw_parameter_values(priors, constraints, seed)
        for name, prior in priors_by_name.items():
            if prior.category in ("N", "T"):
                assert values[name] == int(values[name]), (
                    f"{name} (catégorie {prior.category}) devrait être un "
                    f"entier, obtenu {values[name]}"
                )
        if values["ra"] != int(values["ra"]):
            non_integer_ra_seen = True

    assert non_integer_ra_seen, "ra (catégorie A) ne devrait jamais être arrondi"


class _FakeRng:
    """rng minimal pour contrôler exactement la valeur tirée par
    _draw_one_value, sans dépendre de random.Random."""

    def __init__(self, value):
        self._value = value

    def uniform(self, min_, max_):
        return self._value


def test_draw_one_value_round_half_up_not_banker():
    """L'arrondi doit être round-half-up (floor(0.5+x), comme
    particuleC.cpp), PAS round-half-to-even (comportement par défaut de
    la fonction round() de Python) -- cas où les deux méthodes divergent
    : x=2.5 -> 3 en round-half-up, mais round(2.5) == 2 en Python."""
    prior = Prior(name="N1", category="N", law="UN", bounds=(0.0, 10.0, 0.0, 0.0))

    value = _draw_one_value(prior, _FakeRng(2.5))

    assert value == 3.0
    assert (
        round(2.5) == 2
    )  # confirme que Python round() aurait donné un résultat différent


def test_draw_scenario_uses_weight(header_text):
    """Sur les 6 scénarios de human (poids ~1/6 chacun), un tirage sur de
    nombreuses graines doit produire les 6 indices possibles -- preuve
    que le poids est bien utilisé pour sélectionner le scénario, pas un
    choix fixe (particuleC.cpp::ParticleC::drawscenario)."""
    scenarios = parse_header_scenarios(header_text)

    drawn_indices = {draw_scenario(scenarios, seed).index for seed in range(1, 201)}

    assert drawn_indices == {1, 2, 3, 4, 5, 6}


def test_draw_scenario_degenerate_weight():
    """Si un seul scénario a un poids non nul, il doit toujours être
    tiré, quelle que soit la graine."""
    scenarios = [
        Scenario(index=1, weight=1.0, initial_pop_size_exprs=["N1"]),
        Scenario(index=2, weight=0.0, initial_pop_size_exprs=["N1"]),
    ]

    for seed in range(1, 20):
        assert draw_scenario(scenarios, seed).index == 1


def test_draw_scenario_fallback_when_weights_sum_below_one():
    """Si la somme des poids est < 1 (ne devrait pas arriver avec un vrai
    header.txt DIYABC, mais le C++ ne le vérifie pas), le DERNIER
    scénario sert de secours pour tout tirage au-delà de la somme
    cumulée -- même comportement que la boucle C++ (bornée à
    nscenarios-1)."""
    scenarios = [
        Scenario(index=1, weight=0.1, initial_pop_size_exprs=["N1"]),
        Scenario(index=2, weight=0.1, initial_pop_size_exprs=["N1"]),
    ]

    # seed=2 -> rng.random() == 0.956..., largement > 0.2 (somme des poids)
    drawn = draw_scenario(scenarios, seed=2)
    assert drawn.index == 2


def test_evaluate_expression():
    values = {"t1": 12.3, "t2": 4881.0, "d3": 35.0}

    assert evaluate_expression("t1", values) == 12.3
    assert evaluate_expression("0", values) == 0.0
    assert evaluate_expression("t2-d3", values) == 4881.0 - 35.0
    assert evaluate_expression("t2+d3", values) == 4881.0 + 35.0

    with pytest.raises(ValueError):
        evaluate_expression("inconnu", values)


def test_build_demography_scenario1(header_text):
    """Vérifie que build_demography produit la bonne structure
    d'événements pour le scénario 1 de human, avec des valeurs de
    paramètres fixées à la main (pas de tirage aléatoire ici, pour
    isoler le test de la logique de construction de la démographie)."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)
    """Extrait les priors et les contraintes d'ordre de header.txt.

 

    Retourne (priors, constraints). Une ligne qui ne correspond à aucun

    des deux formats connus lève une erreur explicite plutôt que d'être

    silencieusement ignorée : contrairement aux événements de scénario, on

    n'a pas de raison de s'attendre à du vocabulaire non géré ici pour le

    dataset human.

    """
    # Valeurs choisies à la main, cohérentes avec les contraintes
    # (t4 > t3 > t2 > t2-d3, t2-d4 ; t3 > t3-d34)
    values = {
        "N1": 50000,
        "N2": 50000,
        "N3": 50000,
        "N4": 50000,
        "t1": 10,
        "t2": 5000,
        "d3": 30,
        "Nbn3": 200,
        "d4": 20,
        "Nbn4": 300,
        "N34": 60000,
        "t3": 8000,
        "d34": 25,
        "Nbn34": 250,
        "t4": 9000,
        "Na": 40000,
    }

    demography = build_demography(scenario1, values)

    # 4 populations créées
    assert len(demography.populations) == 4
    assert {p.name for p in demography.populations} == {"pop1", "pop2", "pop3", "pop4"}

    # Les événements de fusion sont bien présents, avec les bons temps
    splits = [
        e
        for e in demography.events
        if isinstance(e, __import__("msprime").demography.PopulationSplit)
    ]
    assert len(splits) == 3

    split_by_time = {s.time: s for s in splits}
    assert (
        split_by_time[10].derived == ["pop1"] and split_by_time[10].ancestral == "pop2"
    )
    assert (
        split_by_time[5000].derived == ["pop4"]
        and split_by_time[5000].ancestral == "pop3"
    )
    assert (
        split_by_time[8000].derived == ["pop3"]
        and split_by_time[8000].ancestral == "pop2"
    )


def test_parse_loci_description(header_text):
    """Vérifie le parsing de la section 'loci description' de human,
    format condensé à un seul type d'héritage."""
    description = parse_loci_description(header_text)

    assert description.total_loci == {"A": 5000}
    assert description.group == "G1"
    assert description.start_index == 0  # "from 1" en 1-based -> 0 en 0-based


def test_pipeline_scenario1(header_text):
    """Vérifie que le pipeline complet (header.txt -> Demography) fonctionne
    de bout en bout sur le scénario 1, et que la démographie produite a la
    structure attendue (4 populations, 3 fusions)."""

    demography, values = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )

    assert len(demography.populations) == 4

    splits = [
        e
        for e in demography.events
        if isinstance(e, msprime.demography.PopulationSplit)
    ]
    assert len(splits) == 3

    # Les valeurs tirées doivent inclure tous les paramètres du header
    assert "N1" in values
    assert "t1" in values


OBSERVED_SNP_FILE = REFERENCE_DIR / "human_snp_all22chr_maf5.snp"


def test_count_samples_per_population():
    """Vérifie que le comptage retrouve bien les 4 populations à 30
    individus chacune, annoncées en commentaire dans le fichier."""
    counts = count_samples_per_population(OBSERVED_SNP_FILE)

    assert set(counts.keys()) == {"ASW", "YRI", "CHB", "GBR"}
    assert all(n == 30 for n in counts.values())


def test_population_index_to_name():
    """Vérifie le mapping indice de scénario (1-indexed) -> nom réel de
    population, dans l'ordre d'apparition du fichier .snp."""
    mapping = population_index_to_name(OBSERVED_SNP_FILE)

    assert mapping == {1: "ASW", 2: "YRI", 3: "CHB", 4: "GBR"}


def test_simulate_independent_loci_scenario1(header_text):
    """Vérifie que build_samples_argument construit bien le dict attendu
    par msprime.sim_ancestry, avec les bons noms de populations et le bon
    nombre d'individus par population."""

    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)

    # 4 populations, 30 individus chacune attendus
    assert samples == {"pop1": 30, "pop2": 30, "pop3": 30, "pop4": 30}

    num_loci = 10  # petit nombre pour un test rapide, pas les 51250 réels
    tree_sequences = list(
        simulate_independent_loci(demography, samples, num_loci=num_loci, seed=123)
    )

    # On doit obtenir exactement num_loci arbres indépendants
    assert len(tree_sequences) == num_loci

    # Chaque arbre doit avoir le bon nombre total de lignées échantillonnées
    # (30 individus x 4 populations x ploidy 2 = 240 lignées)
    for ts in tree_sequences:
        assert ts.num_samples == 30 * 4 * 2


def test_simulate_snp_genotypes_scenario1(header_text):
    """Vérifie que chaque locus simulé est polymorphe (au moins un 0 et
    un 1 parmi les génotypes), garantissant la propriété centrale de
    l'algorithme de Hudson : exactement une mutation par locus, jamais
    un locus monomorphe."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)

    num_loci = 20
    tree_sequences = simulate_independent_loci(demography, samples, num_loci, seed=123)
    genotypes_per_locus = list(simulate_snp_genotypes(tree_sequences, seed=456))

    assert len(genotypes_per_locus) == num_loci

    for locus_genotypes in genotypes_per_locus:
        all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
        assert set(all_genotypes) == {0, 1}, f"Locus non polymorphe : {locus_genotypes}"


def test_simulate_snp_genotypes_grouped_by_population(header_text):
    """Vérifie que les génotypes sont bien regroupés par nom de
    population (pop1..pop4), avec le bon nombre de lignées par groupe
    (30 individus x ploidy 2 = 60 lignées par population), et que chaque
    locus reste globalement polymorphe."""
    demography, _ = build_random_demography_for_scenario_index(
        header_text, scenario_index=1, seed=42
    )
    samples = build_samples_argument(OBSERVED_SNP_FILE)

    num_loci = 10
    tree_sequences = simulate_independent_loci(demography, samples, num_loci, seed=123)
    genotypes_per_locus = list(simulate_snp_genotypes(tree_sequences, seed=456))

    assert len(genotypes_per_locus) == num_loci

    for locus_genotypes in genotypes_per_locus:
        assert set(locus_genotypes.keys()) == {"pop1", "pop2", "pop3", "pop4"}
        for _pop_name, genos in locus_genotypes.items():
            assert len(genos) == 60  # 30 individus x ploidy 2

        # Polymorphe globalement (au moins un 0 et un 1 sur l'ensemble)
        all_genotypes = [g for genos in locus_genotypes.values() for g in genos]
        assert set(all_genotypes) == {0, 1}


def test_run_poc_for_directory():
    """Vérifie le point d'entrée de haut niveau : à partir d'un simple
    chemin de dossier (comme le -p ./ de DIYABC), tout le pipeline doit
    fonctionner sans qu'on ait à lire manuellement header.txt ou le
    fichier .snp nous-mêmes."""
    mutated, values = run_poc_for_directory(
        REFERENCE_DIR,
        scenario_index=1,
        num_loci=15,
        seed=42,
    )

    mutated_list = list(mutated)
    assert len(mutated_list) == 15
    assert "N1" in values


def test_read_header_text_prefers_header_txt(tmp_path):
    """Si les deux fichiers sont présents, header.txt doit être lu en
    priorité (config initiale fournie par l'utilisateur), pas
    headerRF.txt (variante produite par un run DIYABC réel)."""
    (tmp_path / "header.txt").write_text("contenu header.txt")
    (tmp_path / "headerRF.txt").write_text("contenu headerRF.txt")

    assert read_header_text(tmp_path) == "contenu header.txt"


def test_read_header_text_falls_back_to_headerRF(tmp_path):
    """Si seul headerRF.txt est présent (ex: reference/Exemple5/), il
    doit être lu en repli."""
    (tmp_path / "headerRF.txt").write_text("contenu headerRF.txt")

    assert read_header_text(tmp_path) == "contenu headerRF.txt"


def test_simulate_reference_directory(tmp_path):
    """Vérifie le point d'entrée pour un sous-dossier de test sous
    reference/ : à partir d'un dossier ne contenant qu'un header.txt et
    le .snp observé, doit tirer les scénarios pondérés et écrire
    reftable_msprime.txt DANS CE MÊME DOSSIER (jamais 'reftable.txt',
    pour ne pas être confondu avec le first_records_of_the_reference_
    table_0.txt d'un vrai run DIYABC)."""
    shutil.copy(REFERENCE_DIR / "header.txt", tmp_path / "header.txt")
    (tmp_path / OBSERVED_SNP_FILE.name).symlink_to(OBSERVED_SNP_FILE)

    results = simulate_reference_directory(tmp_path, num_loci=10, nrec=4)

    assert len(results) == 4
    output_file = tmp_path / "reftable_msprime.txt"
    assert output_file.exists()

    lines = output_file.read_text().splitlines()
    assert lines[0].startswith("   scenario   ")  # centre("scenario", 14)
    assert len(lines) == 5  # en-tête + 4 particules


def test_write_snp_file_small_case(tmp_path):
    """Vérifie l'écriture du fichier .snp sur un cas minimal : 2 loci,
    2 populations, 2 lignées (1 individu) chacune."""
    genotypes_per_locus = [
        {
            "pop1": [0, 1],
            "pop2": [1, 1],
        },  # locus 0 : pop1 hétérozygote (1), pop2 homozygote dérivé (2)
        {
            "pop1": [0, 0],
            "pop2": [0, 1],
        },  # locus 1 : pop1 homozygote ancestral (0), pop2 hétérozygote (1)
    ]

    output_file = tmp_path / "test.snp"
    write_snp_file(genotypes_per_locus, output_file)

    content = output_file.read_text()
    lines = content.strip().splitlines()

    assert lines[0] == "IND SEX POP A A"
    assert lines[1] == "sim_pop1_1 9 pop1 1 0"
    assert lines[2] == "sim_pop2_1 9 pop2 2 1"


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie -- "
    "ce test nécessite le binaire 'general' compilé de DIYABC.",
)
def test_compute_summary_statistics_scenario1(tmp_path):
    """Vérifie que compute_summary_statistics produit bien les 112
    statistiques résumées attendues (filtre ALL), en déléguant le calcul
    au vrai binaire C++ sur des données simulées par notre pipeline."""
    summary_statistics, values = compute_summary_statistics(
        reference_directory=REFERENCE_DIR,
        scenario_index=1,
        num_loci=10,
        seed=42,
        general_binary_path=GENERAL_BINARY_PATH,
        work_directory=tmp_path,
        stats_filter="ALL",
    )
    print(sorted(summary_statistics.keys()))
    # 112 statistiques attendues (vu dans "group summary statistics (112)")
    assert len(summary_statistics) == 130

    # Quelques noms de colonnes attendus, parmi les plus simples à vérifier
    assert "ML1p_1" in summary_statistics
    assert "FST1m_1" in summary_statistics

    # Les valeurs de paramètres tirées doivent toujours être présentes
    assert "N1" in values


def test_parse_requested_statistic_names_human_old_vocabulary(header_text):
    """Le parseur est purement syntaxique -- il fonctionne aussi sur le
    vocabulaire ANCIEN (obsolète) de human/header.txt (HP0/HM1/HV1/HMO/
    FP0.../AP0...), même si ces noms ne correspondent à aucune
    statistique calculée par summary_statistics.py (incohérence connue
    de ce header.txt, voir notes/exploration.md)."""
    names = parse_requested_statistic_names(header_text)

    assert len(names) == 112
    assert names[:5] == ["HP0_1", "HP0_2", "HP0_3", "HP0_4", "HM1_1"]
    assert names[-1] == "AMO_4.1.2"


def test_parse_requested_statistic_names_modern_vocabulary():
    """Format condensé moderne (ex: toy_example5_modif) : noms de
    colonnes 'STAT_index', même convention que summary_statistics.py."""
    header_text = (
        "group summary statistics (9)\n"
        "group G1 (9)\n"
        "ML1p 1 2 3\n"
        "ML2p 1.2 1.3 2.3\n"
        "HWm 1 2 3\n"
    )

    names = parse_requested_statistic_names(header_text)

    assert names == [
        "ML1p_1",
        "ML1p_2",
        "ML1p_3",
        "ML2p_1.2",
        "ML2p_1.3",
        "ML2p_2.3",
        "HWm_1",
        "HWm_2",
        "HWm_3",
    ]


def test_parse_requested_statistic_names_missing_section():
    with pytest.raises(ValueError, match="group summary statistics"):
        parse_requested_statistic_names("pas de section ici\n")


def test_parse_requested_statistic_names_count_mismatch():
    header_text = "group summary statistics (5)\ngroup G1 (5)\nML1p 1 2 3\n"

    with pytest.raises(ValueError, match="annonce 5"):
        parse_requested_statistic_names(header_text)


def _replace_group_summary_statistics_section(
    header_text: str, new_section_lines: list[str]
) -> str:
    """Remplace la section 'group summary statistics' de header_text par
    new_section_lines -- même approche que loci_parser.rewrite_loci_count
    pour tester avec un contenu différent sans maintenir un header.txt
    séparé à la main."""
    lines = header_text.splitlines()
    start = next(
        i
        for i, line in enumerate(lines)
        if line.strip().startswith("group summary statistics")
    )
    end = start + 1
    while end < len(lines) and lines[end].strip():
        end += 1
    return "\n".join(lines[:start] + new_section_lines + lines[end:])


def test_compute_summary_statistics_stats_filter_header(tmp_path, header_text):
    """stats_filter='HEADER' ne garde, dans l'ordre de déclaration, que
    les statistiques listées dans 'group summary statistics' --
    remplace la section obsolète de human/header.txt par un petit
    sous-ensemble au vocabulaire moderne, pour vérifier le filtrage
    sans dépendre d'un dataset externe."""
    modified_header_text = _replace_group_summary_statistics_section(
        header_text,
        ["group summary statistics (4)", "group G1 (4)", "ML1p 1 2", "HWm 1 2"],
    )
    (tmp_path / "header.txt").write_text(modified_header_text)
    (tmp_path / OBSERVED_SNP_FILE.name).symlink_to(OBSERVED_SNP_FILE)

    summary_stats, values = compute_summary_statistics(
        reference_directory=tmp_path,
        scenario_index=1,
        num_loci=10,
        seed=1,
        stats_filter="HEADER",
    )

    assert list(summary_stats.keys()) == ["ML1p_1", "ML1p_2", "HWm_1", "HWm_2"]
    assert "N1" in values


def test_compute_summary_statistics_stats_filter_header_raises_on_unknown_names(
    header_text,
):
    """stats_filter='HEADER' sur le vrai human/header.txt (vocabulaire
    obsolète HP0/HM1/...) doit lever une ValueError explicite plutôt que
    de produire silencieusement un reftable vide ou incomplet."""
    with pytest.raises(ValueError, match="non calculées"):
        compute_summary_statistics(
            reference_directory=REFERENCE_DIR,
            scenario_index=1,
            num_loci=10,
            seed=1,
            stats_filter="HEADER",
        )


def test_compute_summary_statistics_unknown_stats_filter_raises():
    with pytest.raises(NotImplementedError, match="stats_filter"):
        compute_summary_statistics(
            reference_directory=REFERENCE_DIR,
            scenario_index=1,
            num_loci=10,
            seed=1,
            stats_filter="BOGUS",
        )


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie.",
)
def test_run_reftable_simulation_scenario1(header_text):
    """Vérifie que run_reftable_simulation produit bien nrec particules
    distinctes (tirages de paramètres différents), chacune avec ses 130
    statistiques résumées calculées. Un seul scénario candidat ([scenario1])
    force toutes les particules dessus, isolant ce test du tirage pondéré
    (voir test_run_reftable_simulation_draws_multiple_scenarios pour ça)."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    nrec = 4
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenarios=[scenario1],
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )

    assert len(results) == nrec

    # Les résultats doivent être dans l'ordre des indices de particule
    assert [r.particle_index for r in results] == list(range(nrec))

    # Un seul scénario candidat -> toutes les particules dessus
    assert {r.scenario_index for r in results} == {1}

    # Chaque particule doit avoir ses 130 stats et tous les paramètres
    for result in results:
        assert len(result.summary_statistics) == 130
        assert "N1" in result.parameter_values

    # Les tirages de paramètres doivent être DIFFÉRENTS entre particules
    # (preuve que chaque particule a bien sa propre seed)
    n1_values = [r.parameter_values["N1"] for r in results]
    assert len(set(n1_values)) == nrec  # 4 valeurs distinctes


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie.",
)
def test_run_reftable_simulation_draws_multiple_scenarios(header_text):
    """Avec les 6 scénarios de human en candidats, les particules doivent
    se répartir sur PLUSIEURS scénarios différents (pas toutes sur le
    même) -- preuve que le tirage pondéré par `weight` est bien exercé
    de bout en bout. Répartition exacte pré-calculée pour seed=1..6 :
    scénarios [1, 6, 2, 2, 4, 5] (voir draw_scenario, déterministe)."""
    scenarios = parse_header_scenarios(header_text)

    nrec = 6
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenarios=scenarios,
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )

    assert [r.scenario_index for r in results] == [1, 6, 2, 2, 4, 5]


def test_is_constant_prior():
    """Vérifie la règle de filtrage des priors quasi-constants."""
    # Cas normal : large intervalle, jamais constant
    normal_prior = Prior(
        name="N1", category="N", law="UN", bounds=(1000.0, 100000.0, 0.0, 0.0)
    )
    assert is_constant_prior(normal_prior) is False

    # Cas dégénéré : min == max, clairement constant
    constant_prior = Prior(
        name="X", category="N", law="UN", bounds=(100.0, 100.0, 0.0, 0.0)
    )
    assert is_constant_prior(constant_prior) is True

    # Cas limite : différence infime, sous le seuil
    near_constant_prior = Prior(
        name="Y", category="N", law="UN", bounds=(100.0, 100.00001, 0.0, 0.0)
    )
    assert is_constant_prior(near_constant_prior) is True

    # Cas limite inverse : différence juste au-dessus du seuil
    barely_variable_prior = Prior(
        name="Z", category="N", law="UN", bounds=(100.0, 100.1, 0.0, 0.0)
    )
    assert is_constant_prior(barely_variable_prior) is False


def test_write_reftable_bin(tmp_path, header_text):
    """Vérifie l'écriture du reftable.bin, et sa relecture (vérification
    manuelle du format, sans dépendre de readReftable.R pour ce test
    Python -- juste une vérification structurelle du binaire produit)."""
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    nrec = 3
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenarios=[scenario1],
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )

    output_file = tmp_path / "reftable.bin"
    write_reftable_bin(results, priors, [scenario1], output_file)

    assert output_file.exists()
    assert output_file.stat().st_size > 0


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie.",
)
def test_write_reftable_bin_multi_scenario(tmp_path, header_text):
    """Vérifie la structure du reftable.bin quand plusieurs scénarios
    sont candidats : nscen/nrecscen/nparam doivent porter sur TOUS les
    scénarios candidats (pas seulement ceux tirés), et chaque
    enregistrement doit avoir une longueur VARIABLE = nparam[scenario-1]
    floats -- pas d'union de colonnes ni de NA écrite sur disque
    (vérifié empiriquement contre reftable.cpp et un vrai reftableRF.bin
    multi-scénario -- voir write_reftable_bin)."""
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)

    nrec = 6
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenarios=scenarios,
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )
    # Répartition déterministe pré-calculée (voir test_run_reftable_simulation_draws_multiple_scenarios)
    assert [r.scenario_index for r in results] == [1, 6, 2, 2, 4, 5]

    output_file = tmp_path / "reftable_multi.bin"
    write_reftable_bin(results, priors, scenarios, output_file)

    expected_nrecscen = [
        sum(1 for r in results if r.scenario_index == s.index) for s in scenarios
    ]
    expected_nparam = [
        sum(
            1
            for p in priors
            if not is_constant_prior(p)
            and p.name in get_parameter_names_used_by_scenario(s)
        )
        for s in scenarios
    ]

    data = output_file.read_bytes()
    offset = 0

    def read_int():
        nonlocal offset
        value = struct.unpack_from("<i", data, offset)[0]
        offset += 4
        return value

    assert read_int() == nrec
    assert read_int() == len(scenarios)  # nscen == 6, tous les candidats

    nrecscen_read = [read_int() for _ in scenarios]
    assert nrecscen_read == expected_nrecscen

    nparam_read = [read_int() for _ in scenarios]
    assert nparam_read == expected_nparam

    nstat_read = read_int()
    assert nstat_read == 130

    nparam_by_index = dict(zip((s.index for s in scenarios), nparam_read, strict=True))

    # Chaque ligne : longueur VARIABLE, exactement nparam[son_scenario]
    # floats de paramètres puis nstat floats de stats -- rien de plus.
    for result in results:
        assert read_int() == result.scenario_index
        offset += 4 * nparam_by_index[result.scenario_index]  # params
        offset += 4 * nstat_read  # stats

    assert offset == len(data)  # tout le fichier consommé, aucun octet en trop


def test_extract_referenced_names():
    """Vérifie l'extraction de noms sur des cas simples."""
    assert extract_referenced_names("t1") == {"t1"}
    assert extract_referenced_names("0") == set()
    assert extract_referenced_names("t2-d3") == {"t2", "d3"}
    assert extract_referenced_names("t2+d3") == {"t2", "d3"}


def test_get_parameter_names_used_by_scenario1(header_text):
    """Vérifie que le scénario 1 référence bien exactement les 16
    paramètres attendus (21 priors déclarés au total dans header.txt,
    moins ra/t11/t22/t33/t44 qui appartiennent aux scénarios 2-6)."""
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    used_names = get_parameter_names_used_by_scenario(scenario1)

    expected = {
        "N1",
        "N2",
        "N3",
        "N4",
        "t1",
        "t2",
        "d3",
        "Nbn3",
        "d4",
        "Nbn4",
        "N34",
        "t3",
        "d34",
        "Nbn34",
        "t4",
        "Na",
    }
    assert used_names == expected
    assert len(used_names) == 16


@pytest.mark.skipif(
    GENERAL_BINARY_PATH is None,
    reason="Variable d'environnement DIYABC_GENERAL_PATH non définie.",
)
def test_write_reftable_txt_header_lowercase_and_real_value_for_unused_params(
    tmp_path, header_text
):
    """Vérifie que l'en-tête utilise 'scenario' en minuscule (pas
    'Scenario'), et que les paramètres non pertinents au scénario tiré
    d'une ligne contiennent quand même leur valeur RÉELLEMENT TIRÉE --
    JAMAIS une case vide.

    Une case vide ne produit aucun token pour un parseur par espaces
    (ex: pandas read_csv(sep=r'\\s+'), ou un simple line.split()), ce
    qui décale d'une colonne TOUTES les valeurs suivantes sur la ligne
    -- bug découvert empiriquement en comparant un vrai reftable DIYABC
    à notre sortie sur toy_example5_modif : la colonne 'r' (non utilisée
    par le seul scénario actif) était laissée en blanc, décalant les 51
    statistiques suivantes et produisant des "écarts" massifs qui
    n'avaient rien à voir avec un vrai écart de simulation.

    Les particules sont toutes tirées sur le scénario 1 (candidat
    unique), mais write_reftable_txt reçoit les 6 scénarios de header.txt
    comme `scenarios` -- exerce l'union des colonnes (21 params, dont
    'ra' et t11..t44 propres aux scénarios 2-6)."""
    priors, _ = parse_priors(header_text)
    all_scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in all_scenarios if s.index == 1)

    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenarios=[scenario1],
        num_loci=10,
        nrec=2,
        stats_filter="ALL",
    )

    output_file = tmp_path / "reftable.txt"
    write_reftable_txt(results, priors, all_scenarios, output_file)

    lines = output_file.read_text().splitlines()
    header_line = lines[0]

    assert header_line.startswith("   scenario   ")  # centre("scenario", 14)
    assert "Scenario" not in header_line

    # Repérage par position de TOKEN (pas par offset de caractères) :
    # une valeur N peut occuper plus de 12 caractères (ex: 100000.000000
    # = 13 caractères) et décaler les offsets fixes -- sans casser le
    # découpage par espaces, puisque "  " est toujours réinjecté avant
    # chaque valeur (voir write_reftable_txt).
    header_tokens = header_line.split()
    ra_index = header_tokens.index("ra")

    data_lines = lines[1:]
    assert len(data_lines) == len(results)
    for line, result in zip(data_lines, results, strict=True):
        tokens = line.split()
        # Test anti-régression : autant de tokens que l'en-tête, sinon
        # un parseur par espaces décale toutes les colonnes suivantes
        # (bug corrigé : une case vide ne produit aucun token).
        assert len(tokens) == len(header_tokens)
        assert float(tokens[ra_index]) == pytest.approx(result.parameter_values["ra"])
