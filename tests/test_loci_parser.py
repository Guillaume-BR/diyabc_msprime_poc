"""Vérifie loci_parser : parsing de la section 'loci description' de
header.txt."""

from bridge.loci_parser import parse_loci_description


def test_parse_loci_description(header_text):
    """Vérifie le parsing de la section 'loci description' de human,
    format condensé à un seul type d'héritage."""
    description = parse_loci_description(header_text)

    assert description.total_loci == {"A": 5000}
    assert description.group == "G1"
    assert description.start_index == 0  # "from 1" en 1-based -> 0 en 0-based


def test_parse_loci_description_detailed(header_text_te2):
    """Vérifie le parsing de la section 'loci description' de toy_example2,
    format détaillé avec une seul type d'héritage avec des microsat et des
    séquences dna."""
    list_loci = parse_loci_description(header_text_te2)

    assert len(list_loci) == 20
    assert list_loci[0].name == "Locus_M_A_1_"
    assert list_loci[0].heritage == "A"
    assert list_loci[0].ms_or_seq == "M"
    assert list_loci[0].group == "G1"
    assert list_loci[0].motif_size == 2
    assert list_loci[0].motif_range == 40
    assert list_loci[19].dnalength == 100
