"""
Parseur du format statobsRF.txt / statobs.txt produit par
HeaderC::calstatobs (binaire C++ `general`) -- format texte simple :
une ligne d'en-têtes (noms de colonnes alignés), une ligne de valeurs
(notation scientifique), toutes deux séparées par des espaces.

Référence : observé empiriquement en lançant `general -R "FST1;ML1" -r 1`
sur un dossier de test (voir notes/exploration.md pour la validation de
cette approche par délégation au C++).
"""


def parse_statobs(statobs_text: str) -> dict[str, float]:
    """Parse le contenu d'un fichier statobsRF.txt/statobs.txt en dict
    {nom_colonne: valeur}.

    Lève ValueError si le fichier ne contient pas exactement 2 lignes
    non vides, ou si le nombre de noms et de valeurs ne correspond pas.
    """
    lines = [line for line in statobs_text.splitlines() if line.strip()]
    if len(lines) != 2:
        raise ValueError(
            f"Format inattendu : {len(lines)} ligne(s) non vide(s) trouvée(s), "
            f"2 attendues (en-têtes puis valeurs)."
        )

    names = lines[0].split()
    values = [float(v) for v in lines[1].split()]

    if len(names) != len(values):
        raise ValueError(
            f"Nombre de noms ({len(names)}) différent du nombre de "
            f"valeurs ({len(values)})."
        )

    return dict(zip(names, values))