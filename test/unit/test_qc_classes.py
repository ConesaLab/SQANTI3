import pytest, sys, os
from unittest.mock import mock_open, patch,MagicMock
from bx.intervals import Interval # type: ignore
from collections import defaultdict

# If the path to where the main sqanti3 directory is not in the system path, our modules wont be loaded
main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)
from src.qc_classes import genePredReader, genePredRecord, myQueryTranscripts, myQueryProteins, CAGEPeak, PolyAPeak


@pytest.fixture
def mock_file():
    mock = mock_open(read_data="mocked line\n")
    with patch("builtins.open", mock):
        yield mock


### genePredReader tests ###
def test_genePredReader_initialization(mock_file):
    filename = "test_file.txt"
    reader = genePredReader(filename)
    mock_file.assert_called_once_with(filename)
    assert reader.f is not None

def test_genePredReader_iter(mock_file):
    filename = "test_file.txt"
    reader = genePredReader(filename)
    assert iter(reader) is reader

def test_genePredReader_next_valid_line(mock_file):
    filename = "test_file.txt"
    with patch("src.qc_classes.genePredRecord.from_line") as mock_from_line:
        mock_from_line.return_value = MagicMock(name="genePredRecord")
        mock_file.return_value.read.return_value = "mocked line\n"
        
        reader = genePredReader(filename)
        next_item = next(reader)
        mock_from_line.assert_called_once_with("mocked line")
        assert isinstance(next_item, MagicMock)

def test_genePredReader_next_end_of_file(mock_file):
    filename = "test_file.txt"
    mock_file.return_value.readline.return_value = ""
    reader = genePredReader(filename)
    with pytest.raises(StopIteration):
        next(reader)

def test_genePredReader_multiple_lines(mock_file):
    filename = "test_file.txt"
    mock_file.return_value.readline.side_effect = ["mocked line 1\n", "mocked line 2\n", ""]
    with patch("src.qc_classes.genePredRecord.from_line") as mock_from_line:
        mock_from_line.return_value = MagicMock(name="genePredRecord")
        
        reader = genePredReader(filename)
        item1 = next(reader)
        mock_from_line.assert_called_with("mocked line 1")
        
        item2 = next(reader)
        mock_from_line.assert_called_with("mocked line 2")
        
        with pytest.raises(StopIteration):
            next(reader)

def test_genePredReader_file_not_found():
    filename = "nonexistent_file.txt"
    with pytest.raises(FileNotFoundError):
        reader = genePredReader(filename)

 ### genePredRecord tests ###
def test_genePredRecord_initialization():
    record = genePredRecord(
        id="gene1",
        chrom="chr1",
        strand="+",
        txStart=100,
        txEnd=200,
        cdsStart=120,
        cdsEnd=180,
        exonCount=2,
        exonStarts=[100, 150],
        exonEnds=[120, 200]
    )

    assert record.id == "gene1"
    assert record.chrom == "chr1"
    assert record.strand == "+"
    assert record.txStart == 100
    assert record.txEnd == 200
    assert record.cdsStart == 120
    assert record.cdsEnd == 180
    assert record.exonCount == 2
    assert record.exonStarts == [100, 150]
    assert record.exonEnds == [120, 200]
    assert record.length == 70  # 120-100 + 200-150
    assert len(record.exons) == 2
    assert isinstance(record.exons[0], Interval)
    assert record.junctions == [(120, 150)]

def test_genePredRecord_segments():
    record = genePredRecord(
        id="gene1",
        chrom="chr1",
        strand="+",
        txStart=100,
        txEnd=200,
        cdsStart=120,
        cdsEnd=180,
        exonCount=2,
        exonStarts=[100, 150],
        exonEnds=[120, 200]
    )

    segments = record.segments
    assert len(segments) == 2
    assert isinstance(segments[0], Interval)
    assert segments[0].start == 100
    assert segments[0].end == 120
    assert segments[1].start == 150
    assert segments[1].end == 200

def test_genePredRecord_from_line():
    line = "gene1\tchr1\t+\t100\t200\t120\t180\t2\t100,150,\t120,200,\t.\tGene1"
    
    record = genePredRecord.from_line(line)

    assert record.id == "gene1"
    assert record.chrom == "chr1"
    assert record.strand == "+"
    assert record.txStart == 100
    assert record.txEnd == 200
    assert record.cdsStart == 120
    assert record.cdsEnd == 180
    assert record.exonCount == 2
    assert record.exonStarts == [100, 150]
    assert record.exonEnds == [120, 200]
    assert record.gene == "Gene1"

from Bio import SeqRecord

def test_genePredRecord_get_splice_site():
    # Create a mock genome dictionary with a sequence
    genome_dict = {
        "chr1": SeqRecord.SeqRecord(seq="AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC")
    }
    record = genePredRecord(
        id="gene1",
        chrom="chr1",
        strand="+",
        txStart=10,
        txEnd=20,
        cdsStart=12,
        cdsEnd=18,
        exonCount=2,
        exonStarts=[10, 15],
        exonEnds=[12, 20]
    )

    # Test the splice site for the first junction (index 0)
    splice_site = record.get_splice_site(genome_dict, 0)
    assert splice_site == "AGGC"

## ERROR CASES

# Init methods
def test_genePredRecord_negative_exonCount():
    with pytest.raises(ValueError, match="Exon count must be a positive integer."):
        genePredRecord("gene1", "chr1", "+", 100, 200, 50, 150, -1, [100, 150], [200, 250])

def test_genePredRecord_exon_start_after_end():
    with pytest.raises(ValueError):
        genePredRecord("gene1", "chr1", "+", 100, 200, 50, 150, 2, [200, 150], [250, 300])

def test_genePredRecord_exon_count_mismatch():
    with pytest.raises(ValueError):
        genePredRecord("gene1", "chr1", "+", 100, 200, 50, 150, 2, [100], [200, 250])

def test_genePredRecord_transcript_exons_mismatch():
    with pytest.raises(ValueError):  # Chromosome 'chr2' does not exist in genome_dict
        genePredRecord("gene1", "chr2", "+", 100, 200, 50, 150, 2, [100, 150], [200, 250])

# get_splice_site
def test_genePredRecord_get_splice_site_invalid_index():
    gene_pred = genePredRecord("gene1", "chr1", "+", 100, 250, 50, 150, 2, [100, 150], [200, 250])
    genome_dict = {"chr1": MagicMock(seq="ATGC" * 100)}  # Mock a genome sequence
    with pytest.raises(AssertionError):  # Junction index out of range
        gene_pred.get_splice_site(genome_dict, 2)

def test_genePredRecord_get_splice_site_chrom_not_found():
    gene_pred = genePredRecord("gene1", "chr2", "+", 100, 250, 50, 150, 2, [100, 150], [200, 250])
    genome_dict = {"chr1": MagicMock(seq="ATGC" * 100)}  # Mock a genome sequence for chr1
    with pytest.raises(KeyError):  # Chromosome 'chr2' does not exist in genome_dict
        gene_pred.get_splice_site(genome_dict, 0)

### myQueryTranscript tests ###
## Class creation tests

def test_myQueryTranscripts_init():
    obj = myQueryTranscripts(
        isoform="transcript1",
        diff_to_TSS=10,
        diff_to_TTS=20,
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1", "gene2"],
        transcripts=["transcript1"],
        chrom="chr1",
        strand="+"
    )

    assert obj.isoform == "transcript1"
    assert obj.diff_to_TSS == 10
    assert obj.diff_to_TTS == 20
    assert obj.exons == 3
    assert obj.length == 1000
    assert obj.structural_category == "full-length"
    assert obj.genes == ["gene1", "gene2"]
    assert obj.transcripts == ["transcript1"]
    assert obj.chrom == "chr1"
    assert obj.strand == "+"


def test_get_total_diff():
    obj = myQueryTranscripts(isoform="transcript1", diff_to_TSS=10, diff_to_TTS=-20, exons=3, length=1000, structural_category="full-length")
    assert obj.get_total_diff() == 30


def test_update():
    obj = myQueryTranscripts(isoform="transcript1", diff_to_TSS=10, diff_to_TTS=20, exons=3, length=1000, structural_category="full-length")
    obj.update({
        "transcripts":["ref_trans1"],
        "genes":["gene1"],
        "diff_to_TSS":5,
        "diff_to_TTS":15,
        "ref_length":900,
        "ref_exons":2
    }
    )

    assert obj.transcripts == ["ref_trans1"]
    assert obj.genes == ["gene1"]
    assert obj.diff_to_TSS == 5
    assert obj.diff_to_TTS == 15
    assert obj.ref_length == 900
    assert obj.ref_exons == 2

def test_geneName_single():
    obj = myQueryTranscripts(
        isoform="transcript1",
        diff_to_TSS=10,
        diff_to_TTS=20,
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1"]
    )
    assert obj.geneName() == "gene1"

def test_geneName_multi():
    obj = myQueryTranscripts(
        isoform="transcript1",
        diff_to_TSS=10,
        diff_to_TTS=20,
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1", "gene2", "gene1"]
    )
    assert obj.geneName() == "gene1_gene2"


def test_ratioExp():
    obj = myQueryTranscripts(
        isoform="transcript1",
        diff_to_TSS=10,
        diff_to_TTS=20,
        exons=3,
        length=1000,
        structural_category="full-length",
        iso_exp=10,
        gene_exp=50
    )
    assert obj.ratioExp() == 0.2

    obj.gene_exp = 0
    assert obj.ratioExp() == None


def test_CDSlen():
    obj = myQueryTranscripts(
        isoform="transcript1",
        diff_to_TSS=10,
        diff_to_TTS=20,
        CDS_genomic_start=100,
        CDS_genomic_end=300,
        exons=3,
        length=1000,
        structural_category="full-length",
        CDS_start=100,
        CDS_end=300,
        coding="coding"
    )
    assert obj.get_orf_size() == 201

    obj.coding = "non_coding"
    assert obj.get_orf_size() == None


def test_as_dict():
    obj = myQueryTranscripts(
        isoform="transcript1",
        diff_to_TSS=10,
        diff_to_TTS=20,
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1", "gene2"],
        transcripts=["transcript1"],
        chrom="chr1",
        strand="+",
        iso_exp=10,
        gene_exp=50
    )
    d = obj.as_dict()

    assert d["isoform"] == "transcript1"
    assert d["chrom"] == "chr1"
    assert d["strand"] == "+"
    assert d["length"] == 1000
    assert d["exons"] == 3
    assert d["associated_gene"] == "gene1_gene2"
    assert d["ratio_exp"] == 0.2


def test_as_dict_single_sample_FL():
    """Test as_dict with single-sample FL count"""
    obj = myQueryTranscripts(
        isoform="transcript1",
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1"],
        transcripts=["transcript1"],
        FL=100
    )
    d = obj.as_dict()
    
    # Single sample: FL should be the value, no FL.{sample} columns
    assert d["FL"] == 100
    assert not any(key.startswith("FL.") for key in d.keys())


def test_as_dict_multi_sample_FL():
    """Test as_dict with multi-sample FL counts"""
    obj = myQueryTranscripts(
        isoform="transcript1",
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1"],
        transcripts=["transcript1"],
        FL_dict={"BioSample_1": 3, "BioSample_2": 5, "BioSample_3": 2}
    )
    d = obj.as_dict()
    
    # Multi-sample: should have FL.{sample} columns and FL as sum
    assert d["FL.BioSample_1"] == 3
    assert d["FL.BioSample_2"] == 5
    assert d["FL.BioSample_3"] == 2
    assert d["FL"] == 10  # sum of all samples
    assert "FL_dict" not in d  # FL_dict should be removed


def test_as_dict_multi_sample_FL_empty():
    """Test as_dict with empty FL_dict"""
    obj = myQueryTranscripts(
        isoform="transcript1",
        exons=3,
        length=1000,
        structural_category="full-length",
        genes=["gene1"],
        transcripts=["transcript1"]
    )
    d = obj.as_dict()
    
    # Empty FL_dict: FL should be "NA" (from None conversion)
    assert d["FL"] == "NA"
    assert not any(key.startswith("FL.") for key in d.keys())


## ERROR cases
def test_invalid_data_types():
    with pytest.raises(TypeError):
        myQueryTranscripts(
            isoform=123
        )

def test_empty_required_fields():
    with pytest.raises(ValueError):
        myQueryTranscripts(
            isoform="",  # Empty id
            diff_to_TSS=10,
            diff_to_TTS=20,
            exons=3,
            length=1000,
            structural_category="full-length"
        )

## TODO: Check if this depends on the strand
def test_invalid_CDS():
    with pytest.raises(ValueError):
        obj = myQueryTranscripts(
            isoform="transcript1",
            diff_to_TSS=10,
            diff_to_TTS=20,
            exons=3,
            length=1000,
            strand="+",
            structural_category="full-length",
            CDS_start=300,
            CDS_end=100,
            coding="coding"
        )

### myQueryProteins ###


def test_myQueryProteins_validate_input():
    with pytest.raises(ValueError):
        myQueryProteins(cds_start=-1, cds_end=100, protein_length=99)  # Negative CDS start

    with pytest.raises(ValueError):
        myQueryProteins(cds_start=100, cds_end=50, protein_length=49)  # CDS end < CDS start

    with pytest.raises(ValueError):
        myQueryProteins(cds_start=1, cds_end=100, protein_length=-50)  # Negative ORF length

    obj = myQueryProteins(cds_start=1, cds_end=100, protein_length=100)
    assert obj.cds_start == 1
    assert obj.cds_end == 100
    assert obj.cds_length == 100

### CAGEPeak ###


# Mock data for testing
@pytest.fixture
def mock_bed_file():
    return f"{sys.path[0]}/test/test_data/bedfiles/mock_peaks.bed"

# Tests for CAGEPeak
def test_cagepeak_init(mock_bed_file):
    cage = CAGEPeak(mock_bed_file)
    assert isinstance(cage.cage_peaks, defaultdict)
    assert (("chr1", "+") in cage.cage_peaks)
    #assert len(cage.cage_peaks[("chr1", "+")]) == 2

def test_cagepeak_find(mock_bed_file):
    cage = CAGEPeak(mock_bed_file)
    assert cage.find("chr1", "+", 17) == ("TRUE", 0)  # Query within peak in TSS
    assert cage.find("chr1", "+", 20) == ("TRUE", -3)  # Within peak but downstream TSS
    assert cage.find("chr1", "+", 14) == ("TRUE", +3)  # Within peak but upstream TSS (degradation)
    assert cage.find("chr1", "+", 8) == ("FALSE","NA")  # Outside peak uptream of TSS
    assert cage.find("chr1", "+", 25) == ("FALSE", -8)  # Outside peak downstream of TSS
    assert cage.find("chr1", "-", 17) == ("FALSE", -15)  # Inside peak but on opposite strand (close to another TSS)
    assert cage.find("chr1","+",105) == ("FALSE",-59) # This is to check that it is chromosme specific
    assert cage.find("chr2","+",201) == ("FALSE",-1) # Cage peak with size 1
    assert cage.find("chr2","+",200) == ("TRUE",0) # Cage peak with size 1
    assert cage.find("chr2","-",200) == ("FALSE",-1) # Cage peak with size 1
    assert cage.find("chr2","-",200001) == ("FALSE","NA") # Cage peak with size 1


def test_cagepeak_invalid_file():
    with pytest.raises(FileNotFoundError):
        CAGEPeak("non_existent_file.bed")

# Tests for PolyAPeak
def test_polyapeak_init(mock_bed_file):
    polya = PolyAPeak(mock_bed_file)
    assert isinstance(polya.polya_peaks, defaultdict)
    assert (("chr1", "+") in polya.polya_peaks)

def test_polyapeak_find(mock_bed_file):
    polya = PolyAPeak(mock_bed_file)
    assert polya.find("chr1", "+", 13) == ("TRUE", 0)  # Query within peak in 5'
    assert polya.find("chr1", "+", 16) == ("TRUE", -3)  # Within peak but donwstream 5'
    assert polya.find("chr1", "+", 8) == ("FALSE",5)  # Outside peak upstream of TSS
    assert polya.find("chr1", "+", 25) == ("FALSE", -12)  # Outside peak downstream of TSS
    assert polya.find("chr1", "-", 17) == ("FALSE", -25)  # Inside peak but on opposite strand (close to another TSS)
    assert polya.find("chr1","+",105) == ("FALSE",-64) # This is to check that it is chromosme specific

def test_polyapeak_invalid_file():
    with pytest.raises(FileNotFoundError):
        PolyAPeak("non_existent_file.bed")

# # Additional edge case tests
def test_cagepeak_empty_bed_file(tmp_path):
    empty_file = tmp_path / "empty.bed"
    empty_file.write_text("")
    cage = CAGEPeak(str(empty_file))
    assert len(cage.cage_peaks) == 0

def test_polyapeak_empty_bed_file(tmp_path):
    empty_file = tmp_path / "empty.bed"
    empty_file.write_text("")
    polya = PolyAPeak(str(empty_file))
    assert len(polya.polya_peaks) == 0


### prevalence tests ###
# Number of samples in which each transcript is detected.
# A sample expresses a transcript when its count is above 0 (--min_expression 0,
# the default) or reaches --min_expression; see the contract tests below.


def _transcript_with_counts(fl_dict):
    """Build a minimal transcript carrying a multi-sample FL count dict."""
    obj = myQueryTranscripts(
        isoform="transcript1",
        length=1000,
        exons=3,
        structural_category="full-splice_match",
        genes=["gene1"],
        transcripts=["transcript1"],
    )
    obj.FL_dict = fl_dict
    return obj


def test_prevalence_counts_samples_with_expression():
    """prevalence counts how many samples reach the detection threshold."""
    d = _transcript_with_counts({"s1": 10, "s2": 0, "s3": 5, "s4": 0, "s5": 1}).as_dict()
    assert d["prevalence"] == 3
    assert d["FL"] == 16


def test_prevalence_all_samples():
    """A transcript detected everywhere gets prevalence equal to sample count."""
    d = _transcript_with_counts({"s1": 3, "s2": 7, "s3": 2}).as_dict()
    assert d["prevalence"] == 3


def test_prevalence_no_expression():
    """No expression in any sample gives prevalence 0."""
    d = _transcript_with_counts({"s1": 0, "s2": 0, "s3": 0}).as_dict()
    assert d["prevalence"] == 0
    assert d["FL"] == 0


def test_prevalence_fractional_counts_are_expression():
    """Contract: by default any count above 0 is expression, fractional or not.

    EM-based quantifiers (bambu, salmon, kallisto, RSEM) distribute ambiguous
    reads among isoforms, so an expressed isoform can have less than one read
    in a sample. The default keeps those counts; --min_expression raises the bar.
    """
    d = _transcript_with_counts({"s1": 0.4, "s2": 0.9, "s3": 0, "s4": 1.5}).as_dict()
    assert d["prevalence"] == 3


def test_prevalence_min_expression_is_inclusive():
    """With a min_expression other than 0, a count equal to it is expression."""
    obj = _transcript_with_counts({"s1": 1.0, "s2": 1.5, "s3": 0.7})
    obj.min_expression = 1
    assert obj.as_dict()["prevalence"] == 2


def test_is_expressed():
    """0 means any count above 0; any other value is a minimum count, inclusive."""
    from src.utils import is_expressed
    assert [is_expressed(c, 0) for c in (0, 0.01, 1)] == [False, True, True]
    assert [is_expressed(c, 0.5) for c in (0.4, 0.5, 2)] == [False, True, True]


def test_prevalence_single_sample():
    """With a single sample prevalence can only be 0 or 1."""
    d = _transcript_with_counts({"s1": 5}).as_dict()
    assert d["prevalence"] == 1


def test_prevalence_absent_without_fl_counts():
    """Without --fl_count there is no FL_dict, so prevalence reports NA.

    This is the default path for most users, who do not supply abundances.
    The column must not break the classification table.
    """
    obj = myQueryTranscripts(
        isoform="transcript1",
        length=1000,
        exons=3,
        structural_category="full-splice_match",
        genes=["gene1"],
        transcripts=["transcript1"],
    )
    assert obj.as_dict()["prevalence"] == "NA"


def test_prevalence_empty_fl_counts():
    """An empty FL_dict is treated the same as no FL_dict at all."""
    assert _transcript_with_counts({}).as_dict()["prevalence"] == "NA"


def test_prevalence_uses_shared_expression_threshold():
    """prevalence and prevalence_<group> count expression with the same threshold."""
    obj = _transcript_with_counts({"s1": 0.4, "s2": 0.9, "s3": 1.0})
    obj.counts_design = {"A": ["s1", "s2"], "B": ["s3"]}
    obj.min_expression = 0.9
    d = obj.as_dict()
    assert (d["prevalence"], d["prevalence_A"], d["prevalence_B"]) == (2, 1, 1)


def test_default_min_expression_is_zero():
    """Contract: the default keeps every count above 0 (see test_prevalence_fractional_counts_are_expression)."""
    from src.config import MIN_EXPRESSION
    assert MIN_EXPRESSION == 0
    assert _transcript_with_counts({"s1": 1}).min_expression == 0


def _transcript_with_design(fl_dict, design):
    obj = _transcript_with_counts(fl_dict)
    obj.counts_design = design
    return obj


def test_group_prevalence_per_group():
    """Each group counts detections only over its own samples."""
    d = _transcript_with_design({"K1": 3, "K2": 0, "K3": 1, "B1": 0, "B2": 0.7},
                                {"K": ["K1", "K2", "K3"], "B": ["B1", "B2"]}).as_dict()
    assert d["prevalence_K"] == 2
    assert d["prevalence_B"] == 1
    assert d["prevalence"] == 3


def test_group_prevalence_ignores_samples_out_of_design():
    """Samples not assigned to a group count for prevalence but not for any group."""
    d = _transcript_with_design({"K1": 3, "B1": 2, "MIX": 9},
                                {"K": ["K1"], "B": ["B1"]}).as_dict()
    assert (d["prevalence"], d["prevalence_K"], d["prevalence_B"]) == (3, 1, 1)


def test_group_prevalence_na_without_counts():
    """An isoform missing from --fl_count has NA in every prevalence column, global and per group."""
    d = _transcript_with_design({}, {"K": ["K1"], "B": ["B1"]}).as_dict()
    assert (d["prevalence"], d["prevalence_K"], d["prevalence_B"]) == ("NA", "NA", "NA")


def test_group_prevalence_absent_without_design():
    d = _transcript_with_counts({"K1": 3, "B1": 2}).as_dict()
    assert not [k for k in d if k.startswith("prevalence_")]
    assert "counts_design" not in d
    assert "min_expression" not in d
