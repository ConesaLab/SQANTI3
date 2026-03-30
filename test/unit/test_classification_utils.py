import pytest,sys, os
from collections import namedtuple

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.classification_utils import (
    calc_overlap, calc_splicesite_agreement, calc_exon_overlap, categorize_incomplete_matches, categorize_full_matches,
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

class TranscriptRecord:
    def __init__(self, txStart, txEnd, strand):
        self.txStart = txStart
        self.txEnd = txEnd
        self.strand = strand
        self.exons = [Exon(txStart, txStart+100), Exon(txStart+200, txEnd)] if txEnd-txStart > 300 else [Exon(txStart, txEnd)]
        self.exonCount = len(self.exons)

def test_get_diff_tss_tts_positive_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='+')
    ref = TranscriptRecord(txStart=900, txEnd=2100, strand='+')

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    assert diff_tss == -100
    assert diff_tts == -100
    assert diff_tss_genomic == -100
    assert diff_tts_genomic == -100

def test_get_diff_tss_tts_negative_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='-')
    ref = TranscriptRecord(txStart=900, txEnd=2100, strand='-')

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    assert diff_tss == -100
    assert diff_tts == -100

def test_get_diff_tss_tts_zero_diff_positive_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='+')
    ref = TranscriptRecord(txStart=1000, txEnd=2000, strand='+')

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    assert diff_tss == 0
    assert diff_tts == 0

def test_get_diff_tss_tts_zero_diff_negative_strand():
    trec = TranscriptRecord(txStart=1000, txEnd=2000, strand='-')
    ref = TranscriptRecord(txStart=1000, txEnd=2000, strand='-')

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    assert diff_tss == 0
    assert diff_tts == 0

def test_get_diff_tss_tts_positive_diff_positive_strand():
    trec = TranscriptRecord(txStart=1100, txEnd=2000, strand='+')
    ref = TranscriptRecord(txStart=1000, txEnd=1900, strand='+')

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    assert diff_tss == -100
    assert diff_tts == 100

def test_get_diff_tss_tts_positive_diff_negative_strand():
    trec = TranscriptRecord(txStart=1100, txEnd=2000, strand='-')
    ref = TranscriptRecord(txStart=1000, txEnd=1900, strand='-')

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    assert diff_tss == 100
    assert diff_tts == -100


# New test: spliced and genomic distances differ
def test_get_diff_tss_tts_spliced_vs_genomic():
    # trec: two exons, ref: two exons, but with a gap in the middle
    class CustomTranscript:
        def __init__(self, txStart, txEnd, strand, exons):
            self.txStart = txStart
            self.txEnd = txEnd
            self.strand = strand
            self.exons = exons
            self.exonCount = len(exons)

    # Genomic: trec starts at 1000, ref at 900, so genomic diff_tss = -100
    # But exons are arranged so that the spliced distance is different
    trec_exons = [Exon(1000, 1100), Exon(1200, 1300)]
    ref_exons = [Exon(900, 950), Exon(1000,1100), Exon(1200, 1300)]
    trec = CustomTranscript(1000, 1300, '+', trec_exons)
    ref = CustomTranscript(900, 1300, '+', ref_exons)

    diff_tss, diff_tts, diff_tss_genomic, diff_tts_genomic = get_diff_tss_tts(trec, ref)
    # Genomic TSS: 900-1000 = -100
    # Spliced TSS: measure along ref exons from 900 to 1000 (since diff_tss_genomic < 0)
    # Only the first ref exon covers 900-950, so spliced = -(1000-900-(1000-950)) = -50
    # But our _exonic_len_between sums overlap: min(950,1000)-max(900,900)=50 (900-950), min(1300,1000)-max(1200,1000)=0
    # So spliced = -50
    assert diff_tss_genomic == -100
    assert diff_tss == -50
    # TTS: trec ends at 1300, ref at 1300, so both are 0
    assert diff_tts_genomic == 0
    assert diff_tts == 0

## get_gene_diff_tss_tts ##

TRec = namedtuple('TRec', ['txStart', 'txEnd', 'strand'])
@pytest.fixture
def isoform_hit():
    return myQueryTranscripts(isoform="gene1", 
                       exons=0,
                       length=0,
                       chrom="chr1",
                       strand="+")

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

    assert isoform_hit.diff_to_gene_TSS == -100
    assert isoform_hit.diff_to_gene_TTS == 100

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

    assert isoform_hit.diff_to_gene_TSS == 100
    assert isoform_hit.diff_to_gene_TTS == -100

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

    assert isoform_hit.diff_to_gene_TSS == -50
    assert isoform_hit.diff_to_gene_TTS == 50

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

    assert isoform_hit.diff_to_gene_TSS is None
    assert isoform_hit.diff_to_gene_TTS is None

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

    assert isoform_hit.diff_to_gene_TSS == 0
    assert isoform_hit.diff_to_gene_TTS == 0


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
    
def test_categorize_full_matches(diff_tss, diff_tts, expected_subtype):
    assert categorize_full_matches(diff_tss, diff_tts) == expected_subtype

def test_edge_cases():
    assert categorize_full_matches(50, 50) == 'reference_match'
    assert categorize_full_matches(50, 51) == 'alternative_3end'
    assert categorize_full_matches(51, 50) == 'alternative_5end'
    assert categorize_full_matches(51, 51) == 'alternative_3end5end'

def test_large_values():
    assert categorize_full_matches(1000, 1000) == 'alternative_3end5end'
    assert categorize_full_matches(-1000, -1000) == 'alternative_3end5end'

def test_mixed_signs():
    assert categorize_full_matches(-25, 25) == 'reference_match'
    assert categorize_full_matches(-75, 75) == 'alternative_3end5end'
    assert categorize_full_matches(-75, 25) == 'alternative_5end'
    assert categorize_full_matches(-25, 75) == 'alternative_3end'

def test_zero_values():
    assert categorize_full_matches(0, 0) == 'reference_match'
    assert categorize_full_matches(0, 51) == 'alternative_3end'
    assert categorize_full_matches(51, 0) == 'alternative_5end'