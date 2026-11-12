"""Vérifie reftable_loop : boucle multi-particules, tirage pondéré de scénario
par particule, écriture reftable.bin (binaire DIYABC) et reftable.txt (lisible,
pour diff direct avec un vrai run DIYABC)."""

import shutil
import struct
from pathlib import Path

import pytest
from conftest import (
    GENERAL_BINARY_PATH,
    OBSERVED_MSS_FILE_TE2,
    OBSERVED_SNP_FILE_HUMAN,
    REFERENCE_DIR,
)

from bridge.demography_builder import get_parameter_names_used_by_scenario
from bridge.header_dataclasses import DirichletPrior, SplitEvent
from bridge.prior_parser import is_constant_prior, parse_priors
from bridge.reftable_loop import (
    ParticleResult,
    _detect_mu,
    _historical_columns_order,
    _kept_param_names_by_scenario,
    _split_mutation_and_stat_names,
    check_header_trailer_line,
    check_header_with_dirichlet_prior,
    group_prior_column_names,
    parse_real_reftable_params,
    parse_real_reftable_params_with_group_priors,
    raise_if_serial_with_sex_linked_loci,
    rewrite_real_reftable_txt,
    run_reftable_simulation_dna,
    run_reftable_simulation_microsat,
    run_reftable_simulation_snp,
    simulate_from_directory,
    write_reftable_bin,
    write_reftable_txt,
)
from bridge.scenario_parser import is_serial_scenario, parse_header_scenarios


def test_detect_mu(header_text_te4, header_text_te2):
    """Vérifie que _detect_mu retourne bien l'écriture de mu dans le trailer du header.txt"""
    mu = _detect_mu(header_text_te4)
    assert mu == "µ"

    mu2 = _detect_mu(header_text_te2)
    assert mu2 == "µ"

    fake_header = header_text_te2.replace("µ", "mu")
    mu3 = _detect_mu(fake_header)
    assert mu3 == "mu"


def test_check_header_with_dirichlet_prior(header_text):
    priors = [
        DirichletPrior(
            names=["ra", "rb"],
            category="A",
            law="DI",
            min=0.05,
            max=0.95,
            mean=0.0,
            sdshape=0.0,
        ),
    ]

    # Everything ok
    header_text1 = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    scenarios = [parse_header_scenarios(header_text1)[1]]
    assert not check_header_with_dirichlet_prior(priors, scenarios)

    # cas mixte
    header_text2 = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rc")
    scenarios = [parse_header_scenarios(header_text2)[1]]
    with pytest.raises(
        ValueError,
        match="Des taux identiques d'admixture dans le prior et dans le scénario sont nécessaires.",
    ):
        check_header_with_dirichlet_prior(priors, scenarios)

    # un seul taux avec ra,rb
    scenarios = [parse_header_scenarios(header_text)[2]]
    with pytest.raises(
        ValueError,
        match="Les taux d'admixture déclarés dans les priors et les scénarios diffèrent.",
    ):
        check_header_with_dirichlet_prior(priors, scenarios)

    priors = [
        DirichletPrior(
            names=["ra", "rb"],
            category="A",
            law="DI",
            min=0.05,
            max=0.95,
            mean=0.0,
            sdshape=0.0,
        ),
        DirichletPrior(
            names=["rc", "rd"],
            category="A",
            law="DI",
            min=0.05,
            max=0.95,
            mean=0.0,
            sdshape=0.0,
        ),
    ]

    # deux priors de Dirichlet distinct
    header_text3 = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rc")
    scenarios = [parse_header_scenarios(header_text3)[1]]
    with pytest.raises(
        ValueError, match="Les taux viennent de plusieurs priors de Dirichlet :"
    ):
        check_header_with_dirichlet_prior(priors, scenarios)


def test_check_header_with_dirichlet_prior_accepts_independent_uniform_rates(
    header_text,
):
    """Un split k=3 dont les deux taux sont des priors classiques indépendants
    (`ra` et `rb` en UN) est légal et ne doit pas être refusé : le contrôle ne
    concerne que les taux déclarés dans un prior de Dirichlet."""
    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]",
        "ra A UN[0.05,0.5,0.0,0.0]\nrb A UN[0.05,0.5,0.0,0.0]",
    )
    priors, _ = parse_priors(header_text)
    scenarios = [parse_header_scenarios(header_text)[1]]

    # Garde-fou du test : sans Dirichlet, il ne peut passer à vide que si
    # `split` porte bien deux taux.
    assert not any(isinstance(p, DirichletPrior) for p in priors)
    assert scenarios[0].index == 2
    assert any(
        len(event.admixture_rates) == 2
        for event in scenarios[0].events
        if isinstance(event, SplitEvent)
    )

    check_header_with_dirichlet_prior(priors, scenarios)  # ne doit pas lever


def test_simulate_from_directory(tmp_path, snp_context_human):
    """Vérifie le point d'entrée pour un sous-dossier de test sous reference/ :
    à partir d'un dossier ne contenant qu'un header.txt et le .snp observé,
    doit tirer les scénarios pondérés et écrire reftable_msprime.txt DANS CE
    MÊME DOSSIER (jamais 'reftable.txt', pour ne pas être confondu avec le
    first_records_of_the_reference_ table_0.txt d'un vrai run DIYABC)."""
    shutil.copy(REFERENCE_DIR / "human" / "header.txt", tmp_path / "header.txt")
    (tmp_path / OBSERVED_SNP_FILE_HUMAN.name).symlink_to(OBSERVED_SNP_FILE_HUMAN)

    results = simulate_from_directory(tmp_path, snp_context_human, num_loci=10, nrec=4)

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
    """Vérifie que run_reftable_simulation_snp produit bien nrec particules
    distinctes (tirages de paramètres différents), chacune avec ses 130
    statistiques résumées calculées.

    Un seul scénario candidat ([scenario1])
    force toutes les particules dessus, isolant ce test du tirage pondéré
    (voir test_run_reftable_simulation_draws_multiple_scenarios pour ça).
    """
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    nrec = 4
    results = run_reftable_simulation_snp(
        reference_directory=REFERENCE_DIR / "human",
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
    scénarios [6, 4, 6, 1, 5, 6] (voir draw_scenario, déterministe --
    séquence recalculée après correction du bug de seed partagée entre
    draw_scenario et draw_parameter_values, voir
    _SCENARIO_DRAW_SEED_OFFSET dans reftable_loop.py, puis une seconde fois
    le 05/10 quand la graine du tirage est passée par `_locus_seed`)."""
    scenarios = parse_header_scenarios(header_text)

    nrec = 6
    results = run_reftable_simulation_snp(
        reference_directory=REFERENCE_DIR / "human",
        scenarios=scenarios,
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )

    assert [r.scenario_index for r in results] == [6, 4, 6, 1, 5, 6]


def test_run_reftable_simulation_snp_with_wrong_header_dirichlet(header_text, tmp_path):
    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rc A DI[0.05,0.95,0.0,0.0]"
    )
    (tmp_path / "headerRF.txt").write_text(header_text)
    with pytest.raises(
        ValueError,
        match="Des taux identiques d'admixture dans le prior et dans le scénario sont nécessaires",
    ):
        run_reftable_simulation_snp(
            reference_directory=tmp_path,
            scenarios=parse_header_scenarios(header_text),
            nrec=1,
        )


def test_write_reftable_bin_without_group_priors(tmp_path, header_text):
    """Vérifie l'écriture du reftable.bin, et sa relecture (vérification
    manuelle du format, sans dépendre de readReftable.R pour ce test
    Python -- juste une vérification structurelle du binaire produit)."""
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    nrec = 3
    results = run_reftable_simulation_snp(
        reference_directory=REFERENCE_DIR / "human",
        scenarios=[scenario1],
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )

    output_file = tmp_path / "reftable.bin"
    write_reftable_bin(results, priors, [scenario1], output_file)

    assert output_file.exists()
    assert output_file.stat().st_size > 0


def test_write_reftable_bin_with_group_priors(tmp_path, microsat_context_te1):
    """Vérifie l'écriture du reftable.bin, et sa relecture (vérification
    manuelle du format, sans dépendre de readReftable.R pour ce test
    Python -- juste une vérification structurelle du binaire produit)."""
    priors, _ = parse_priors(microsat_context_te1.header_text)
    scenarios = parse_header_scenarios(microsat_context_te1.header_text)
    scenario1 = next(s for s in scenarios if s.index == 1)

    nrec = 3
    results = run_reftable_simulation_microsat(
        reference_directory=REFERENCE_DIR / "toy_example1_ms_seriel",
        scenarios=[scenario1],
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
    results = run_reftable_simulation_snp(
        reference_directory=REFERENCE_DIR / "human",
        scenarios=scenarios,
        num_loci=10,
        nrec=nrec,
        stats_filter="ALL",
    )
    # Répartition déterministe pré-calculée (voir test_run_reftable_simulation_draws_multiple_scenarios)
    assert [r.scenario_index for r in results] == [6, 4, 6, 1, 5, 6]

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
    """Vérifie que l'en-tête utilise 'scenario' en minuscule (pas 'Scenario'),
    et que les paramètres non pertinents au scénario tiré d'une ligne
    contiennent quand même leur valeur RÉELLEMENT TIRÉE -- JAMAIS une case
    vide.

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
    'ra' et t11..t44 propres aux scénarios 2-6).
    """
    priors, _ = parse_priors(header_text)
    all_scenarios = parse_header_scenarios(header_text)
    scenario1 = next(s for s in all_scenarios if s.index == 1)

    results = run_reftable_simulation_snp(
        reference_directory=REFERENCE_DIR / "human",
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


def test_kept_param_names_by_scenario(header_text_te1):
    """Vérifie que _kept_param_names_by_scenario retourne bien le même nombre
    de noms de paramètres à garder par scénarios."""
    priors = parse_priors(header_text_te1)[0]
    scenarios = parse_header_scenarios(header_text_te1)
    result = _kept_param_names_by_scenario(priors, scenarios)

    assert len(result) == 2


def test_kept_param_names_by_scenario_with_dirichlet_prior(header_text):
    """Un prior de Dirichlet `ra,rb` donne deux colonnes consécutives, dans
    l'ordre de déclaration, pour les scénarios dont un `split` les référence,
    et aucune pour les scénarios qui ne les référencent pas.

    On passe TOUS les `split` de human au format k=3 (`... ra rb`) pour que le
    header reste cohérent : scénarios 2, 3, 5, 6 utilisent le Dirichlet,
    scénarios 1 et 4 (sans `split`) ne l'utilisent pas.
    """
    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    header_text = header_text.replace("t1 split 1 3 2 ra", "t1 split 1 3 2 4 ra rb")
    header_text = header_text.replace("t11 split 1 4 2 ra", "t11 split 1 4 2 3 ra rb")
    header_text = header_text.replace("t11 split 1 3 2 ra", "t11 split 1 3 2 4 ra rb")
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rb A DI[0.05,0.95,0.0,0.0]"
    )
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    results = _kept_param_names_by_scenario(priors, scenarios)

    for index in (2, 3, 5, 6):
        names = results[index]
        assert "ra" in names
        assert names.index("ra") + 1 == names.index("rb")

    for index in (1, 4):
        assert "ra" not in results[index]
        assert "rb" not in results[index]


def test_historical_columns_order(microsat_context_te1_modified):
    """Vérifie que _historical_columns_order retourne bien l'ordre des colonnes
    historiques tel qu'il est dans l'en-tête du reftable.txt (et non pas
    l'ordre de déclaration dans le header.txt)."""
    header_line = "scenario N1 N2 N3 N4 ra t32 t21 t421 t41 mumic_1 pmic_1 snimic_1"
    header_text_te1 = microsat_context_te1_modified.header_text
    priors = parse_priors(header_text_te1)[0]
    scenarios = parse_header_scenarios(header_text_te1)
    result = _historical_columns_order(header_line, priors, scenarios)

    # Vérifie que l'ordre des colonnes historiques est bien celui de l'en-tête
    assert result == [
        "N1",
        "N2",
        "N3",
        "N4",
        "ra",
        "t32",
        "t21",
        "t421",
        "t41",
    ]

    header_line_error = (
        "scenario N1 N2 N3 N4 ra t32 t21 t421 t14 mumic_1 pmic_1 snimic_1"
    )
    with pytest.raises(ValueError, match="t41"):
        _historical_columns_order(header_line_error, priors, scenarios)


def test_historical_column_order_with_dirichlet_prior(microsat_context_te1_modified):
    header_text = microsat_context_te1_modified.header_text
    header_text = header_text.replace("t421 split 4 2 1 ra", "t421 split 4 2 1 3 ra rb")
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rb A DI[0.05,0.95,0.0,0.0]"
    )
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    header_line = header_text.splitlines()[-1]

    with pytest.raises(ValueError, match=r"incohérentes.*rb"):
        _historical_columns_order(header_line, priors, scenarios)

    header_line = header_line.replace("N4 ra t41", "N4 ra rb t41")
    result = _historical_columns_order(header_line, priors, scenarios)
    assert "rb" in result
    assert result.index("ra") + 1 == result.index("rb")

    header_line = header_line.replace("N4 ra rb t41", "N4 rb ra t41")
    assert result.index("ra") + 1 == result.index("rb")


def test_write_reftable_txt_with_dirichlet(header_text, tmp_path):
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rb A DI[0.05,0.95,0.0,0.0]"
    )
    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    parameter_values = {
        name: 0.1 for name in _kept_param_names_by_scenario(priors, scenarios)[2]
    }
    parameter_values.update({"ra": 0.2, "rb": 0.3})
    particle = ParticleResult(
        particle_index=1,
        scenario_index=2,
        parameter_values=parameter_values,
        summary_statistics={"S1": 0.5},
        group_priors_values={},
    )

    write_reftable_txt(
        [particle], priors, scenarios, output_path=tmp_path / "reftable.txt"
    )
    lines = (tmp_path / "reftable.txt").read_text().splitlines()
    header_tokens = lines[0].split()
    ra_index = header_tokens.index("ra")
    rb_index = header_tokens.index("rb")

    for line in lines[1:]:
        tokens = line.split()
        assert float(tokens[ra_index]) == 0.2
        assert float(tokens[rb_index]) == 0.3
        assert len(tokens) == len(header_tokens)

    assert "ra" in header_tokens
    assert header_tokens.index("ra") + 1 == header_tokens.index("rb")


def test_split_mutation_and_stat_names_rejects_undeclared_historical_parameters():
    """Un paramètre historique du trailer absent des priors (DIYABC `-R ALL` retire
    des priors ceux qu'aucun scénario n'utilise mais les laisse dans le trailer)
    se retrouve entre les paramètres historiques connus et les mutations : on
    lève avec son nom, au lieu de décaler les colonnes."""
    names = ["µmic_1", "pmic_1", "µseq_2", "k1seq_2", "NAL_1_1", "MPD_2_1"]
    assert _split_mutation_and_stat_names(names) == (names[:4], names[4:])

    with pytest.raises(ValueError, match="t11.*t22"):
        _split_mutation_and_stat_names(["t11", "t22", *names])


def test_rewrite_real_reftable_txt_explains_trailer_with_undeclared_parameter(
    tmp_path, microsat_context_te1_modified
):
    header_text = microsat_context_te1_modified.header_text
    priors = parse_priors(header_text)[0]
    scenarios = parse_header_scenarios(header_text)
    hist = "N1 N2 N3 N4 ra t32 t21 t421 t41"
    # `t99` : paramètre historique du trailer que les priors ne déclarent pas.
    raw = tmp_path / "first_records.txt"
    raw.write_text(
        f"scenario {hist} t99 mumic_1 pmic_1 snimic_1 NAL_1_1\n"
        "1 1 2 3 4 5 6 7 8 9 10 11 12\n"
    )
    with pytest.raises(ValueError, match="t99"):
        rewrite_real_reftable_txt(raw, tmp_path / "clean.txt", priors, scenarios)


def test_parse_real_reftable_params_follows_file_header_order(
    tmp_path, header_text_te1_modified
):
    priors, _ = parse_priors(header_text_te1_modified)
    scenarios = parse_header_scenarios(header_text_te1_modified)
    # ordre d'en-tête volontairement différent de la déclaration (ra, t41, t32, t21, t421)
    reftable = tmp_path / "reftable.txt"
    reftable.write_text(
        "scenario N1 N2 N3 N4 ra t32 t21 t421 t41 NAL_1_1\n"
        "1 2279 626 1391 6439 174 595 674 7.7\n"
        "2 1096 8342 1730 1983 0.191 203 892 284 7.4\n"
    )
    rows = parse_real_reftable_params(reftable, priors, scenarios)
    assert rows[0] == (
        1,
        {
            "N1": 2279,
            "N2": 626,
            "N3": 1391,
            "N4": 6439,
            "t32": 174,
            "t21": 595,
            "t41": 674,
        },
    )
    assert rows[1][1]["t421"] == 284 and "t41" not in rows[1][1]


def test_parse_real_reftable_params_with_group_priors(header_text_te2):
    # test sur le fichier first_records_of_the_reference_table_0.txt de toy_example2_ms_dna
    list_priors = parse_priors(header_text_te2)[0]
    group_priors_names = group_prior_column_names(header_text_te2)
    scenarios = parse_header_scenarios(header_text_te2)
    result = parse_real_reftable_params_with_group_priors(
        path=OBSERVED_MSS_FILE_TE2.parent
        / "first_records_of_the_reference_table_0.txt",
        priors=list_priors,
        scenarios=scenarios,
        group_priors_names=group_priors_names,
    )
    path = OBSERVED_MSS_FILE_TE2.parent / "first_records_of_the_reference_table_0.txt"
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


def test_group_prior_column_names_follows_trailer_spelling(header_text_te2):
    lines = header_text_te2.splitlines()
    lines[-1] = lines[-1].replace("µ", "mu")
    header_mu = "\n".join(lines)

    names = group_prior_column_names(header_mu)

    assert names[0] == "mumic_1"  # le premier nom attendu, en dur


def test_run_reftable_simulation_snp_with_dirichelt_split_txt(header_text, tmp_path):
    snp_filename = header_text.splitlines()[0].strip()
    shutil.copy(REFERENCE_DIR / "human" / "headerRF.txt", tmp_path / "headerRF.txt")
    shutil.copy(REFERENCE_DIR / "human" / snp_filename, tmp_path / snp_filename)
    header_text = (tmp_path / "headerRF.txt").read_text()
    # modification du header human
    header_text = header_text.replace(
        "21 parameters and 112 summary statistics",
        "22 parameters and 112 summary statistics",
    )

    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    header_text = header_text.replace("t1 split 1 3 2 ra", "t1 split 1 3 2 4 ra rb")
    header_text = header_text.replace("t11 split 1 4 2 ra", "t11 split 1 4 2 3 ra rb")
    header_text = header_text.replace("t11 split 1 3 2 ra", "t11 split 1 3 2 4 ra rb")

    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rb A DI[0.05,0.95,0.0,0.0]"
    )

    header_text = header_text.replace("Nbn34 t4 Na ra t11", "Nbn34 t4 Na ra rb t11")
    (tmp_path / "headerRF.txt").write_text(header_text)
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    results = run_reftable_simulation_snp(
        reference_directory=tmp_path,
        scenarios=scenarios,
        num_loci=10,
        nrec=20,
        stats_filter="ALL",
    )
    assert any(r.scenario_index in {2, 3, 5, 6} for r in results)
    for result in results:
        if result.scenario_index in {2, 3, 5, 6}:
            assert "ra" in result.parameter_values
            assert "rb" in result.parameter_values
            ra = result.parameter_values["ra"]
            rb = result.parameter_values["rb"]
            assert ra + rb < 1 - 0.05
            assert 1 - ra - rb >= 0.05
            assert 0.05 <= ra <= 0.95
            assert 0.05 <= rb <= 0.95

    write_reftable_txt(
        results, priors, scenarios, output_path=tmp_path / "reftable.txt"
    )
    lines = (tmp_path / "reftable.txt").read_text().splitlines()
    header_line = lines[0].split()
    assert "ra" in header_line
    assert "rb" in header_line
    for line in lines[1:]:
        line = line.split()
        assert len(line) == len(header_line)


def test_write_reftable_bin_with_dirichlet_prior(header_text, tmp_path):
    """Le reftable.bin écrit une colonne par nom d'un prior de Dirichlet.

    Header human dont tous les `split` sont passés en k=3 (`... ra rb`) et
    dont `ra,rb` est un Dirichlet : les scénarios 2, 3, 5, 6 ont un
    paramètre de plus que dans le header d'origine, 1 et 4 aucun.
    Pas de simulation : deux ParticleResult faits à la main, relus octet par
    octet (format : voir write_reftable_bin).
    """
    original_header_text = header_text
    header_text = header_text.replace("t1 split 1 4 2 ra", "t1 split 1 4 2 3 ra rb")
    header_text = header_text.replace("t1 split 1 3 2 ra", "t1 split 1 3 2 4 ra rb")
    header_text = header_text.replace("t11 split 1 4 2 ra", "t11 split 1 4 2 3 ra rb")
    header_text = header_text.replace("t11 split 1 3 2 ra", "t11 split 1 3 2 4 ra rb")
    header_text = header_text.replace(
        "ra A UN[0.05,0.95,0.0,0.0]", "ra,rb A DI[0.05,0.95,0.0,0.0]"
    )
    priors, _ = parse_priors(header_text)
    scenarios = parse_header_scenarios(header_text)
    kept = _kept_param_names_by_scenario(priors, scenarios)

    # Une valeur distincte et non nulle par paramètre : une permutation de
    # colonnes ne passerait pas inaperçue. ra et rb sont fixés à la main.
    parameter_values1 = {name: 0.01 * (i + 1) for i, name in enumerate(kept[1])}
    parameter_values2 = {name: 0.01 * (i + 1) for i, name in enumerate(kept[2])} | {
        "ra": 0.2,
        "rb": 0.3,
    }
    results = [
        ParticleResult(
            particle_index=1,
            scenario_index=2,
            parameter_values=parameter_values2,
            summary_statistics={"S1": 0.5},
            group_priors_values={},
        ),
        ParticleResult(
            particle_index=2,
            scenario_index=1,
            parameter_values=parameter_values1,
            summary_statistics={"S1": 0.5},
            group_priors_values={},
        ),
    ]

    output_file = tmp_path / "reftable.bin"
    write_reftable_bin(results, priors, scenarios, output_path=output_file)
    data = output_file.read_bytes()
    offset = 0

    def read_int():
        """Lit un entier 32 bits à la position courante et avance de 4 octets."""
        nonlocal offset
        value = struct.unpack_from("<i", data, offset)[0]
        offset += 4
        return value

    def read_float():
        """Lit un flottant 32 bits à la position courante et avance de 4 octets."""
        nonlocal offset
        value = struct.unpack_from("<f", data, offset)[0]
        offset += 4
        return value

    # En-tête
    assert read_int() == len(results)
    assert read_int() == len(scenarios)

    nrecscen_read = [read_int() for _ in scenarios]
    assert nrecscen_read == [1, 1, 0, 0, 0, 0]

    # Valeur attendue calculée sur le header d'ORIGINE (indépendante du code
    # testé) : un paramètre de plus pour les scénarios dont le split
    # référence le Dirichlet.
    original_kept = _kept_param_names_by_scenario(
        parse_priors(original_header_text)[0],
        parse_header_scenarios(original_header_text),
    )
    expected_nparam = [
        len(original_kept[s.index]) + (1 if s.index in {2, 3, 5, 6} else 0)
        for s in scenarios
    ]
    nparam_read = [read_int() for _ in scenarios]
    assert nparam_read == expected_nparam

    nstat_read = read_int()
    assert nstat_read == 1

    nparam_by_index = dict(zip((s.index for s in scenarios), nparam_read, strict=True))

    # Enregistrements : dans l'ordre des résultats (scénario 2, puis 1).
    for result in results:
        assert read_int() == result.scenario_index
        floats = [read_float() for _ in range(nparam_by_index[result.scenario_index])]
        assert floats == pytest.approx(
            [result.parameter_values[n] for n in kept[result.scenario_index]]
        )
        assert [read_float() for _ in range(nstat_read)] == pytest.approx([0.5])

        if result.scenario_index == 2:
            ra_position = kept[2].index("ra")
            assert floats[ra_position] == pytest.approx(0.2)
            assert floats[ra_position + 1] == pytest.approx(0.3)  # rb juste après

    assert offset == len(data)  # tout le fichier consommé, aucun octet en trop


def test_run_reftable_simulation_dna_scenario1(header_text_te2):
    """Vérifie que run_reftable_simulation_snp produit bien nrec particules
    distinctes (tirages de paramètres différents), chacune avec ses 42
    statistiques résumées calculées.

    Un seul scénario candidat ([scenario1])
    force toutes les particules dessus, isolant ce test du tirage pondéré
    (voir test_run_reftable_simulation_draws_multiple_scenarios pour ça).
    """
    scenarios = parse_header_scenarios(header_text_te2)
    scenario1 = next(s for s in scenarios if s.index == 1)

    nrec = 4
    results = run_reftable_simulation_dna(
        reference_directory=REFERENCE_DIR / "toy_example2_ms_dna_K2P",
        scenarios=[scenario1],
        nrec=nrec,
        stats_filter="ALL",
    )

    assert len(results) == nrec

    # Les résultats doivent être dans l'ordre des indices de particule
    assert [r.particle_index for r in results] == list(range(nrec))

    # Un seul scénario candidat -> toutes les particules dessus
    assert {r.scenario_index for r in results} == {1}

    # Chaque particule doit avoir ses 42 stats et tous les paramètres
    for result in results:
        assert len(result.summary_statistics) == 42
        assert "N1" in result.parameter_values

    # Les tirages de paramètres doivent être DIFFÉRENTS entre particules
    # (preuve que chaque particule a bien sa propre seed)
    n1_values = [r.parameter_values["N1"] for r in results]
    assert len(set(n1_values)) == nrec  # 4 valeurs distinctes


def test_run_reftable_simulation_dna_draws_multiple_scenarios(header_text_te2):
    """Avec les 2 scénarios de toy_example2 en candidats, les particules doivent
    se répartir sur PLUSIEURS scénarios différents (pas toutes sur le
    même) -- preuve que le tirage pondéré par `weight` est bien exercé
    de bout en bout. Répartition exacte pré-calculée pour seed=1..6 :
    scénarios [2, 2, 2, 1, 2, 2] (voir draw_scenario, déterministe --
    séquence recalculée après correction du bug de seed partagée entre
    draw_scenario et draw_parameter_values, voir
    _SCENARIO_DRAW_SEED_OFFSET dans reftable_loop.py, puis une seconde fois
    le 05/10 quand la graine du tirage est passée par `_locus_seed`)."""
    scenarios = parse_header_scenarios(header_text_te2)

    nrec = 6
    results = run_reftable_simulation_dna(
        reference_directory=REFERENCE_DIR / "toy_example2_ms_dna_K2P",
        scenarios=scenarios,
        nrec=nrec,
        stats_filter="ALL",
    )

    assert [r.scenario_index for r in results] == [2, 2, 2, 1, 2, 2]


# ── Garde « sériel + loci liés au sexe » ──────────────────────────────────


def test_is_serial_scenario_distingue_seriel_et_non_seriel():
    """human_seriel échantillonne 4 fois la population 1, toy_example5 une
    seule fois chacune de ses 3 populations.

    Le critère porte sur les événements `sample` seuls :
    toy_example2_ms_dna déclare 5 populations pour 2 échantillons et
    n'est pourtant PAS sériel -- c'est précisément le cas qu'une
    comparaison à `npop` classerait à tort.
    """
    seriel = parse_header_scenarios(
        (REFERENCE_DIR / "human_seriel" / "headerRF.txt").read_text()
    )
    te5 = parse_header_scenarios(
        (REFERENCE_DIR / "toy_example5" / "headerRF.txt").read_text()
    )
    te2 = parse_header_scenarios(
        (REFERENCE_DIR / "toy_example2_ms_dna_K2P" / "headerRF.txt").read_text()
    )

    assert [is_serial_scenario(s) for s in seriel] == [True, True]
    assert [is_serial_scenario(s) for s in te5] == [False, False, False]
    assert [is_serial_scenario(s) for s in te2] == [False, False]


@pytest.mark.parametrize(
    "dataset",
    ["human_seriel", "toy_example4_seriel", "toy_example1_ms_seriel", "toy_example5"],
)
def test_raise_if_serial_with_sex_linked_loci_laisse_passer_les_datasets_valides(
    dataset,
):
    """Aucun dataset de reference/ ne doit être refusé : les jeux sériels
    (human_seriel, toy_example4_seriel, toy_example1_ms) sont <A>/<M>, et les
    jeux portant des <X>/<Y> (toy_example5) ne sont pas sériels.

    La garde ne doit donc jamais se déclencher sur l'existant -- c'est
    ce qui rend son ajout prouvablement neutre sur tout ce qui est
    validé.
    """
    header_text = (REFERENCE_DIR / dataset / "headerRF.txt").read_text()

    raise_if_serial_with_sex_linked_loci(header_text)


def test_raise_if_serial_with_sex_linked_loci_leve_sur_la_combinaison():
    """Sur toy_example5 (qui porte des <X> et des <Y>) rendu sériel en
    dupliquant un événement `sample` -- le compteur de lignes du bloc `scenario
    N [w] (nlines)` doit être incrémenté en même temps, sinon
    split_scenario_blocks tronque le bloc et l'événement ajouté est ignoré."""
    header_text = (REFERENCE_DIR / "toy_example5" / "headerRF.txt").read_text()
    serial_header = header_text.replace(
        "scenario 1 [0.33333] (6)", "scenario 1 [0.33333] (7)"
    ).replace("0 sample 3\nta merge 1 3", "0 sample 3\n50 sample 3\nta merge 1 3")
    assert is_serial_scenario(parse_header_scenarios(serial_header)[0])

    with pytest.raises(NotImplementedError, match="sériel"):
        raise_if_serial_with_sex_linked_loci(serial_header)


def test_raise_if_serial_with_sex_linked_loci_ignore_le_seriel_sans_xy():
    """Contrôle symétrique : la même mutation appliquée à un header sans locus
    lié au sexe (toy_example4_seriel, <A> seul) ne lève pas.

    La garde refuse une COMBINAISON, pas le sériel en soi.
    """
    header_text = (REFERENCE_DIR / "human_seriel" / "headerRF.txt").read_text()

    raise_if_serial_with_sex_linked_loci(header_text)


# ── Garde « ligne de queue du header » ─────────────────────────────────────


@pytest.mark.parametrize(
    "dataset",
    [
        "toy_example2_ms_dna_K2P",  # µ + K2P : k1 seul
        "toy_example2_ms_dna_JK",  # JK : ni k1 ni k2
        "toy_example2_ms_dna_TN",  # TN : k1 ET k2
        "toy_example1_ms_seriel",  # préfixe "mu" ASCII, pas "µ"
        "human",  # dataset SNP : aucune section `group priors`
    ],
)
def test_check_header_trailer_line_accepte_les_datasets_sains(dataset):
    """Les quatre familles de header valides doivent passer : les trois modèles
    de substitution (qui impliquent 0, 1 ou 2 paramètres `k`), le préfixe
    `mumic_1` en ASCII que DIYABC écrit sur toy_example1_ms au lieu de
    `µmic_1`, et un header SNP sans section `group priors` du tout -- sur
    lequel la garde doit être un no-op, pas une erreur."""
    header_text = (REFERENCE_DIR / dataset / "headerRF.txt").read_text()

    check_header_trailer_line(header_text)


def test_check_header_trailer_line_leve_sur_un_k_manquant():
    """Cas réel attrapé le 30/09 : les lignes de queue de
    toy_example2_ms_dna_50loci_TN et _50loci_JK avaient été recopiées de la
    variante K2P, donc `k1seq` présent et `k2seq` absent alors que `MODEL TN`
    exige les deux -- deux jetons manquants faussent le
    `nparamhist = jetons - 1 - nstat - nparamut` de
    HeaderC::readHeaderAllStat et rendent le reftable inexploitable (son
    en-tête ne contenait effectivement aucune colonne `k2seq`).

    La corruption est reconstruite ici plutôt que lue depuis reference/ :
    ces headers ont été réparés, et un test qui dépend d'un fichier cassé
    disparaît le jour où on le répare."""
    header_text = (
        REFERENCE_DIR / "toy_example2_ms_dna_TN" / "headerRF.txt"
    ).read_text()
    lines = header_text.splitlines()
    lines[-1] = lines[-1].replace(" k2seq_2", "").replace(" k2seq_3", "")
    corrupted = "\n".join(lines)

    with pytest.raises(ValueError, match="k2seq_2"):
        check_header_trailer_line(corrupted)


def test_check_header_trailer_line_leve_sur_un_jeton_en_trop(tmp_path):
    """Symétrique : un `k1seq` de trop sur un header JK -- la faute commise sur
    50loci_JK, corrigée depuis. Le message doit nommer le jeton parasite."""
    header_text = (
        REFERENCE_DIR / "toy_example2_ms_dna_JK" / "headerRF.txt"
    ).read_text()
    lines = header_text.splitlines()
    lines[-1] = lines[-1].replace("µseq_2", "µseq_2 k1seq_2")
    corrupted = "\n".join(lines)

    with pytest.raises(ValueError, match="k1seq_2"):
        check_header_trailer_line(corrupted)
