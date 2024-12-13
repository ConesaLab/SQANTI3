import pytest,sys, os
from collections import namedtuple

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.classification_utils import (
    calc_overlap, calc_splicesite_agreement, calc_exon_overlap
)
def test_calc_overlap():
    assert calc_overlap(1, 5, 3, 7) == 3
    assert calc_overlap(3, 7, 1, 5) == 3
    assert calc_overlap(1, 3, 4, 6) == 0
    assert calc_overlap(1, 5, 1, 5) == 5
    assert calc_overlap('NA', 5, 3, 7) == 0
    assert calc_overlap(1, 5, 'NA', 7) == 0

Exon = namedtuple('Exon', ['start', 'end'])

def test_calc_splicesite_agreement():
    query_exons = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    ref_exons = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    assert calc_splicesite_agreement(query_exons, ref_exons) == 6

    ref_exons = [Exon(100, 200), Exon(300, 450), Exon(550, 600)]
    assert calc_splicesite_agreement(query_exons, ref_exons) == 4

    ref_exons = [Exon(150, 250), Exon(350, 450), Exon(550, 650)]
    assert calc_splicesite_agreement(query_exons, ref_exons) == 0

def test_calc_exon_overlap():
    query_exons = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    ref_exons = [Exon(100, 200), Exon(300, 400), Exon(500, 600)]
    assert calc_exon_overlap(query_exons, ref_exons) == 300

    ref_exons = [Exon(150, 250), Exon(350, 450), Exon(550, 650)]
    assert calc_exon_overlap(query_exons, ref_exons) == 150
    calc_exon_overlap(query_exons, ref_exons)
    ref_exons = [Exon(700, 800)]
    assert calc_exon_overlap(query_exons, ref_exons) == 0

    query_exons = []
    ref_exons = [Exon(100, 200)]
    assert calc_exon_overlap(query_exons, ref_exons) == 0
