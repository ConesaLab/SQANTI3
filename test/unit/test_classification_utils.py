import pytest,sys, os
from collections import namedtuple

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.classification_utils import (
    calc_overlap, calc_splicesite_agreement, calc_exon_overlap, categorize_incomplete_matches, full_splice_match_subtype,
    get_diff_tss_tts, get_gene_diff_tss_tts
)
from src.qc_classes import myQueryTranscripts
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

## get_diff_tss_tts ##


TranscriptRecord = namedtuple('TranscriptRecord', ['txStart', 'txEnd', 'strand'])

def test_get_diff_tss_tts_positive_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='+')
    ref = TranscriptRecord(txStart=900, txEnd=2100, strand='+')

    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
    assert diff_tss == -100
    assert diff_tts == -100

def test_get_diff_tss_tts_negative_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='-')
    ref = TranscriptRecord(txStart=900, txEnd=2100, strand='-')

    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
    assert diff_tss == -100
    assert diff_tts == -100

def test_get_diff_tss_tts_zero_diff_positive_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='+')
    ref = TranscriptRecord(txStart=1000, txEnd=2000, strand='+')

    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
    assert diff_tss == 0
    assert diff_tts == 0

def test_get_diff_tss_tts_zero_diff_negative_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='-')
    ref = TranscriptRecord(txStart=1000, txEnd=2000, strand='-')

    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
    assert diff_tss == 0
    assert diff_tts == 0

def test_get_diff_tss_tts_positive_diff_positive_strand():
    trec = TranscriptRecord(txStart=1100, txEnd=2000, strand='+')
    ref = TranscriptRecord(txStart=1000, txEnd=1900, strand='+')

    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
    assert diff_tss == -100
    assert diff_tts == 100

def test_get_diff_tss_tts_positive_diff_negative_strand():
    trec = TranscriptRecord(txStart=1100, txEnd=2000, strand='-')
    ref = TranscriptRecord(txStart=1000, txEnd=1900, strand='-')

    diff_tss, diff_tts = get_diff_tss_tts(trec, ref)
    assert diff_tss == 100
    assert diff_tts == -100

## get_gene_diff_tss_tts ##

TRec = namedtuple('TRec', ['txStart', 'txEnd', 'strand'])
@pytest.fixture
def isoform_hit():
    return myQueryTranscripts(id="gene1", tts_diff="NA", tss_diff="NA",
                       num_exons=0,
                       length=0,
                       str_class="",
                       chrom="chr1",
                       strand="+",
                       subtype="no_subcategory",
                       percAdownTTS=0,
                       seqAdownTTS=0,
                       genes=[])

# Test cases
def test_positive_strand(isoform_hit):
    isoform_hit.add_gene('gene1')
    trec = TRec(txStart=1000, txEnd=2000, strand='+')
    start_ends_by_gene = {
        'gene1': {
            'begin': [900, 1100],
            'end': [1900, 2100]
        }
    }

    get_gene_diff_tss_tts(isoform_hit, trec, start_ends_by_gene)

    assert isoform_hit.tss_gene_diff == -100
    assert isoform_hit.tts_gene_diff == 100

def test_negative_strand(isoform_hit):
    isoform_hit.add_gene('gene1')
    trec = TRec(txStart=1000, txEnd=2000, strand='-')
    start_ends_by_gene = {
        'gene1': {
            'begin': [900, 1100],
            'end': [1900, 2100]
        }
    }

    get_gene_diff_tss_tts(isoform_hit, trec, start_ends_by_gene)

    assert isoform_hit.tss_gene_diff == 100
    assert isoform_hit.tts_gene_diff == -100

def test_multiple_genes(isoform_hit):
    isoform_hit.add_gene('gene1')
    isoform_hit.add_gene('gene2')
    trec = TRec(txStart=1000, txEnd=2000, strand='+')
    start_ends_by_gene = {
        'gene1': {
            'begin': [900, 1100],
            'end': [1900, 2100]
        },
        'gene2': {
            'begin': [950, 1050],
            'end': [1950, 2050]
        }
    }

    get_gene_diff_tss_tts(isoform_hit, trec, start_ends_by_gene)

    assert isoform_hit.tss_gene_diff == -50
    assert isoform_hit.tts_gene_diff == 50

def test_no_valid_difference(isoform_hit):
    isoform_hit.add_gene('gene1')
    trec = TRec(txStart=1000, txEnd=2000, strand='+')
    start_ends_by_gene = {
        'gene1': {
            'begin': [],
            'end': []
        }
    }

    get_gene_diff_tss_tts(isoform_hit, trec, start_ends_by_gene)

    assert isoform_hit.tss_gene_diff == 'NA'
    assert isoform_hit.tts_gene_diff == 'NA'

def test_exact_match(isoform_hit):
    isoform_hit.add_gene('gene1')
    trec = TRec(txStart=1000, txEnd=2000, strand='+')
    start_ends_by_gene = {
        'gene1': {
            'begin': [1000],
            'end': [2000]
        }
    }

    get_gene_diff_tss_tts(isoform_hit, trec, start_ends_by_gene)

    assert isoform_hit.tss_gene_diff == 0
    assert isoform_hit.tts_gene_diff == 0


### categorize_incomplete_matches ###
# Mock classes
class Exon:
    def __init__(self, start, end):
        self.start = start
        self.end = end

class Transcript:
    def __init__(self, exons, junctions, strand):
        self.exons = exons
        self.junctions = junctions
        self.strand = strand

# Fixtures
@pytest.fixture
def ref_transcript():
    return Transcript(
        exons=[Exon(100, 200), Exon(300, 400), Exon(500, 600)],
        junctions=[(200, 300), (400, 500)],
        strand='+'
    )

# Parametrized test cases
@pytest.mark.parametrize("trec_data, expected_category", [
    (
        {"exons": [Exon(100, 400)], "junctions": [], "strand": "+"},
        "intron_retention"
    ),
    (
        {"exons": [Exon(100, 200), Exon(300, 400), Exon(500, 600)], "junctions": [(200, 300), (400, 500)], "strand": "+"},
        "complete"
    ),
    (
        {"exons": [Exon(300, 400), Exon(500, 600)], "junctions": [(),(400, 500)], "strand": "+"},
        "3prime_fragment"
    ),
    (
        {"exons": [Exon(100, 200), Exon(300, 400)], "junctions": [(200, 300),()], "strand": "+"},
        "5prime_fragment"
    ),
    (
        {"exons": [Exon(300, 400)], "junctions": [(),()], "strand": "+"},
        "internal_fragment"
    ),
])
def test_categorize_incomplete_matches(ref_transcript, trec_data, expected_category):
    trec = Transcript(**trec_data)
    assert categorize_incomplete_matches(trec, ref_transcript) == expected_category

# Test for strand-specific behavior
@pytest.mark.parametrize("strand, expected_category", [
    ("+", "3prime_fragment"),
    ("-", "5prime_fragment"),
])
def test_strand_specific_categorization(ref_transcript, strand, expected_category):
    ref_transcript.strand = strand
    trec = Transcript(
        exons=[Exon(300, 400), Exon(500, 600)],
        junctions=[(400, 500)],
        strand=strand
    )
    assert categorize_incomplete_matches(trec, ref_transcript) == expected_category

# Edge case: empty transcript
def test_empty_transcript(ref_transcript):
    trec = Transcript(exons=[], junctions=[], strand="+")
    with pytest.raises(Exception):  # Adjust the exception type as needed
        categorize_incomplete_matches(trec, ref_transcript)

# Test with different reference transcript
@pytest.fixture
def complex_ref_transcript():
    return Transcript(
        exons=[Exon(100, 200), Exon(300, 400), Exon(500, 600), Exon(700, 800)],
        junctions=[(200, 300), (400, 500), (600, 700)],
        strand='+'
    )

def test_with_complex_reference(complex_ref_transcript):
    trec = Transcript(
        exons=[Exon(300, 400), Exon(500, 600)],
        junctions=[(400, 500)],
        strand='+'
    )
    assert categorize_incomplete_matches(trec, complex_ref_transcript) == "internal_fragment"

## FSM_subtype ##

@pytest.mark.parametrize("diff_tss, diff_tts, expected_subtype", [
(0, 0, 'reference_match'),
(50, 50, 'reference_match'),
(-50, -50, 'reference_match'),
(25, -25, 'reference_match'),

(0, 51, 'alternative_3end'),
(50, 100, 'alternative_3end'),
(-25, -75, 'alternative_3end'),

(51, 0, 'alternative_5end'),
(100, 50, 'alternative_5end'),
(-75, -25, 'alternative_5end'),

(51, 51, 'alternative_3end5end'),
(100, 100, 'alternative_3end5end'),
(-75, -75, 'alternative_3end5end'),
])
    
def test_full_splice_match_subtype(diff_tss, diff_tts, expected_subtype):
    assert full_splice_match_subtype(diff_tss, diff_tts) == expected_subtype

def test_edge_cases():
    assert full_splice_match_subtype(50, 50) == 'reference_match'
    assert full_splice_match_subtype(50, 51) == 'alternative_3end'
    assert full_splice_match_subtype(51, 50) == 'alternative_5end'
    assert full_splice_match_subtype(51, 51) == 'alternative_3end5end'

def test_large_values():
    assert full_splice_match_subtype(1000, 1000) == 'alternative_3end5end'
    assert full_splice_match_subtype(-1000, -1000) == 'alternative_3end5end'

def test_mixed_signs():
    assert full_splice_match_subtype(-25, 25) == 'reference_match'
    assert full_splice_match_subtype(-75, 75) == 'alternative_3end5end'
    assert full_splice_match_subtype(-75, 25) == 'alternative_5end'
    assert full_splice_match_subtype(-25, 75) == 'alternative_3end'

def test_zero_values():
    assert full_splice_match_subtype(0, 0) == 'reference_match'
    assert full_splice_match_subtype(0, 51) == 'alternative_3end'
    assert full_splice_match_subtype(51, 0) == 'alternative_5end'