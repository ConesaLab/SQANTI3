import pytest, sys, os
from unittest.mock import mock_open, patch,MagicMock
from bx.intervals import Interval
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
    print(genome_dict)
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
        id="transcript1",
        tss_diff=10,
        tts_diff=20,
        num_exons=3,
        length=1000,
        str_class="full-length",
        genes=["gene1", "gene2"],
        transcripts=["transcript1"],
        chrom="chr1",
        strand="+"
    )

    assert obj.id == "transcript1"
    assert obj.tss_diff == 10
    assert obj.tts_diff == 20
    assert obj.num_exons == 3
    assert obj.length == 1000
    assert obj.str_class == "full-length"
    assert obj.genes == ["gene1", "gene2"]
    assert obj.transcripts == ["transcript1"]
    assert obj.chrom == "chr1"
    assert obj.strand == "+"


def test_get_total_diff():
    obj = myQueryTranscripts(id="transcript1", tss_diff=10, tts_diff=-20, num_exons=3, length=1000, str_class="full-length")
    assert obj.get_total_diff() == 30


def test_modify():
    obj = myQueryTranscripts(id="transcript1", tss_diff=10, tts_diff=20, num_exons=3, length=1000, str_class="full-length")
    obj.modify(
        ref_transcript="ref_trans1",
        ref_gene="gene1",
        tss_diff=5,
        tts_diff=15,
        refLen=900,
        refExons=2
    )

    assert obj.transcripts == ["ref_trans1"]
    assert obj.genes == ["gene1"]
    assert obj.tss_diff == 5
    assert obj.tts_diff == 15
    assert obj.refLen == 900
    assert obj.refExons == 2

def test_geneName_single():
    obj = myQueryTranscripts(
        id="transcript1",
        tss_diff=10,
        tts_diff=20,
        num_exons=3,
        length=1000,
        str_class="full-length",
        genes=["gene1"]
    )
    assert obj.geneName() == "gene1"

def test_geneName_multi():
    obj = myQueryTranscripts(
        id="transcript1",
        tss_diff=10,
        tts_diff=20,
        num_exons=3,
        length=1000,
        str_class="full-length",
        genes=["gene1", "gene2", "gene1"]
    )
    assert obj.geneName() == "gene1_gene2"


def test_ratioExp():
    obj = myQueryTranscripts(
        id="transcript1",
        tss_diff=10,
        tts_diff=20,
        num_exons=3,
        length=1000,
        str_class="full-length",
        isoExp=10,
        geneExp=50
    )
    assert obj.ratioExp() == 0.2

    obj.geneExp = 0
    assert obj.ratioExp() == "NA"


def test_CDSlen():
    obj = myQueryTranscripts(
        id="transcript1",
        tss_diff=10,
        tts_diff=20,
        num_exons=3,
        length=1000,
        str_class="full-length",
        CDS_start=100,
        CDS_end=300,
        coding="coding"
    )
    assert obj.CDSlen() == "201"

    obj.coding = "non_coding"
    assert obj.CDSlen() == "NA"


def test_as_dict():
    obj = myQueryTranscripts(
        id="transcript1",
        tss_diff=10,
        tts_diff=20,
        num_exons=3,
        length=1000,
        str_class="full-length",
        genes=["gene1", "gene2"],
        transcripts=["transcript1"],
        chrom="chr1",
        strand="+",
        isoExp=10,
        geneExp=50
    )
    d = obj.as_dict()

    assert d["isoform"] == "transcript1"
    assert d["chrom"] == "chr1"
    assert d["strand"] == "+"
    assert d["length"] == 1000
    assert d["exons"] == 3
    assert d["associated_gene"] == "gene1_gene2"
    assert d["ratio_exp"] == 0.2


## ERROR cases
def test_invalid_data_types():
    with pytest.raises(ValueError):
        myQueryTranscripts(
            id=123,  # Invalid type for id
            tss_diff="10",  # Invalid type for tss_diff
            tts_diff=20,
            num_exons=3,
            length=1000,
            str_class="full-length"
        )

def test_empty_required_fields():
    with pytest.raises(ValueError):
        myQueryTranscripts(
            id="",  # Empty id
            tss_diff=10,
            tts_diff=20,
            num_exons=3,
            length=1000,
            str_class="full-length"
        )

## TODO: Check if this depends on the strand
def test_invalid_CDS():
    with pytest.raises(ValueError):
        obj = myQueryTranscripts(
            id="transcript1",
            tss_diff=10,
            tts_diff=20,
            num_exons=3,
            length=1000,
            str_class="full-length",
            CDS_start=300,
            CDS_end=100,
            coding="coding"
        )

### myQueryProteins ###


def test_myQueryProteins_validate_input():
    with pytest.raises(ValueError):
        myQueryProteins(cds_start=-1, cds_end=100, orf_length=99)  # Negative CDS start

    with pytest.raises(ValueError):
        myQueryProteins(cds_start=100, cds_end=50, orf_length=49)  # CDS end < CDS start

    with pytest.raises(ValueError):
        myQueryProteins(cds_start=1, cds_end=100, orf_length=-50)  # Negative ORF length

    obj = myQueryProteins(cds_start=1, cds_end=100, orf_length=100)
    assert obj.cds_start == 1
    assert obj.cds_end == 100
    assert obj.orf_length == 100

### CAGEPeak ###


# Mock data for testing
@pytest.fixture
def mock_bed_file():
    return f"{sys.path[0]}/test/test_data/mock_peaks.bed"

# Tests for CAGEPeak
def test_cagepeak_init(mock_bed_file):
    cage = CAGEPeak(mock_bed_file)
    assert isinstance(cage.cage_peaks, defaultdict)
    assert (("chr1", "+") in cage.cage_peaks)
    #assert len(cage.cage_peaks[("chr1", "+")]) == 2

def test_cagepeak_find(mock_bed_file):
    cage = CAGEPeak(mock_bed_file)
    assert cage.find("chr1", "+", 17) == ("TRUE", 0)  # Query within peak in TSS
    assert cage.find("chr1", "+", 20) == ("TRUE", -3)  # Within peak but donwstream TSS
    assert cage.find("chr1", "+", 14) == ("TRUE", +3)  # Within peak but downstream TSS (degradation)
    assert cage.find("chr1", "+", 8) == ("FALSE","NA")  # Outside peak uptream of TSS
    assert cage.find("chr1", "+", 25) == ("FALSE", -8)  # Outside peak downstream of TSS
    assert cage.find("chr1", "-", 17) == ("FALSE", -39)  # Inside peak but on opposite strand (close to another TSS)
    assert cage.find("chr1","+",105) == ("FALSE",-59) # This is to check that it is chromosme specific

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