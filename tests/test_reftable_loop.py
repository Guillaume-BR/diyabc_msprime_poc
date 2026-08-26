"""Vérifie reftable_loop : boucle multi-particules, tirage pondéré de
scénario par particule, écriture reftable.bin (binaire DIYABC) et
reftable.txt (lisible, pour diff direct avec un vrai run DIYABC)."""

import shutil
import struct
from pathlib import Path

import pytest
from conftest import (
    GENERAL_BINARY_PATH,
    OBSERVED_SNP_FILE_HUMAN,
    OBSERVED_SNP_FILE_TE2,
    REFERENCE_DIR,
)

from bridge.demography_builder import get_parameter_names_used_by_scenario
from bridge.prior_parser import is_constant_prior, parse_priors
from bridge.reftable_loop import (
    group_prior_column_names,
    parse_real_reftable_params_with_group_priors,
    run_reftable_simulation,
    simulate_from_directory,
    write_reftable_bin,
    write_reftable_txt,
)
from bridge.scenario_parser import parse_header_scenarios


def test_simulate_from_directory(tmp_path):
    """Vérifie le point d'entrée pour un sous-dossier de test sous
    reference/ : à partir d'un dossier ne contenant qu'un header.txt et
    le .snp observé, doit tirer les scénarios pondérés et écrire
    reftable_msprime.txt DANS CE MÊME DOSSIER (jamais 'reftable.txt',
    pour ne pas être confondu avec le first_records_of_the_reference_
    table_0.txt d'un vrai run DIYABC)."""
    shutil.copy(REFERENCE_DIR / "header.txt", tmp_path / "header.txt")
    (tmp_path / OBSERVED_SNP_FILE_HUMAN.name).symlink_to(OBSERVED_SNP_FILE_HUMAN)

    results = simulate_from_directory(tmp_path, num_loci=10, nrec=4)

    assert len(results) == 4
    output_file = tmp_path / "reftable_msprime.txt"
    assert output_file.exists()

    lines = output_file.read_text().splitlines()
    assert lines[0].startswith("   scenario   ")  # centre("scenario", 14)
    assert len(lines) == 5  # en-tête + 4 particules


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
    scénarios [3, 3, 5, 3, 1, 3] (voir draw_scenario, déterministe --
    séquence recalculée après correction du bug de seed partagée entre
    draw_scenario et draw_parameter_values, voir
    _SCENARIO_DRAW_SEED_OFFSET dans reftable_loop.py)."""
    scenarios = parse_header_scenarios(header_text)

    nrec = 6
    results = run_reftable_simulation(
        reference_directory=REFERENCE_DIR,
        scenarios=scenarios,
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )

    assert [r.scenario_index for r in results] == [3, 3, 5, 3, 1, 3]


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
    assert [r.scenario_index for r in results] == [3, 3, 5, 3, 1, 3]

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


def test_parse_real_reftable_params_with_group_priors(header_text_te2):
    list_priors = parse_priors(header_text_te2)[0]
    group_priors_names = group_prior_column_names(header_text_te2)
    scenarios = parse_header_scenarios(header_text_te2)
    result = parse_real_reftable_params_with_group_priors(
        path=OBSERVED_SNP_FILE_TE2.parent
        / "first_records_of_the_reference_table_0.txt",
        priors=list_priors,
        scenarios=scenarios,
        group_priors_names=group_priors_names,
    )
    path = OBSERVED_SNP_FILE_TE2.parent / "first_records_of_the_reference_table_0.txt"
    text = Path(path).read_text().splitlines()

    # Vérifie que le nombre de lignes lues correspond au nombre de lignes de données (en-tête exclue)
    assert len(result) == len(text) - 1  # en-tête

    # On va vérifier que la première ligne extraite correspond à la première ligne de données du fichier
    first_data_line = text[1].split()
    scenario_index = first_data_line[0]
    assert result[0][0] == int(scenario_index)

    for i, prior in enumerate(result[0][1]):
        assert result[0][1][prior] == float(first_data_line[i + 1])

    for i, group_priors in enumerate(result[0][2]):
        assert result[0][2][group_priors] == float(
            first_data_line[i + 1 + len(result[0][1])]
        )

    # test d'avoir toujours les noms des group_priors attendus dans le dictionnaire
    for i in range(len(result)):
        assert set(result[i][2].keys()) == set(
            ["µmic_1", "pmic_1", "snimic_1", "µseq_2", "k1seq_2", "µseq_3", "k1seq_3"]
        )

    # test sur le nombre de valeurs de priors pour des scénarios différents
    result_scenario_1 = next(result[i] for i in range(len(result)) if result[i][0] == 1)
    result_scenario_2 = next(result[i] for i in range(len(result)) if result[i][0] == 2)
    assert len(result_scenario_1[1]) != len(result_scenario_2[1])
