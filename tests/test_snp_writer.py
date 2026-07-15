"""Vérifie snp_writer : écriture du fichier .snp au format DIYABC (chemin
subprocess déprécié, gardé pour la validation croisée avec le binaire
réel, voir CLAUDE.md)."""

from bridge.snp_writer import write_snp_file


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
