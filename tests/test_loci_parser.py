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
