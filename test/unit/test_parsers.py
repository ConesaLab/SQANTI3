import logging
import re
import sys,os,pytest
from unittest.mock import mock_open, patch
from Bio import SeqIO

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)
from src.parsers import isoforms_parser, parse_GMST, parse_corrORF, reference_parser
from bx.intervals.intersection import IntervalTree

### reference_parser ###

@pytest.fixture
def reference_parser_input():
    data = {
        "gtf": os.path.join(main_path, "test/test_data/test_reference.gtf"),
        "outdir": os.path.join(main_path, "test/test_data"),
        "prefix": "test",
        "genome_dict": {"chr22": "Sequence"},
        "gene_name": False,
        "isoAnnot": False
    }
    return data
# Create a fixture for the tester logger
@pytest.fixture
def tester_logger():
    # Define a logger for testing
    logger = logging.getLogger("tester_logger")

    logger.setLevel(logging.INFO)
    # Return the logger
    return logger

# Test if the genePred file is created
def test_reference_parser_genePred(reference_parser_input):
    genePred_annot = os.path.join(main_path,reference_parser_input["outdir"],
                                  f"refAnnotation_{reference_parser_input['prefix']}.genePred")
    if os.path.exists(genePred_annot):
        print(f"Removing {genePred_annot}")
        os.remove(genePred_annot)
    reference_parser(*list(reference_parser_input.values()))
    assert os.path.exists(genePred_annot)

# Test if the genePred file is detected
def test_reference_parser_genePred_found(reference_parser_input, caplog,tester_logger):
    genePred_annot = os.path.join(main_path, reference_parser_input["outdir"],
                                   f"refAnnotation_{reference_parser_input['prefix']}.genePred")
    reference_parser(*list(reference_parser_input.values()), logger=tester_logger)
    # Capture the printed output    # Check if the specific print statement is in stdout
    assert f"{genePred_annot} already exists. Using it." in caplog.text
    assert os.path.exists(genePred_annot)

def test_reference_parser_notInGenome(reference_parser_input,caplog, tester_logger):
    reference_parser_input["genome_dict"] = {"chr21": "Sequence"}
    reference_parser(*list(reference_parser_input.values()), logger=tester_logger)
    # Check for the warning message in stderr
    expected_warning = "Reference annotation contains chromosomes not in genome:"
    caplog.set_level(logging.WARNING)
    assert expected_warning in caplog.text

def test_reference_parser_correctOutput_length(reference_parser_input):
    refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene = reference_parser(*list(reference_parser_input.values()))
    # Check if the output is correct
    assert len(refs_1exon_by_chr.keys()) == 1
    assert len(refs_1exon_by_chr["chr22"].find(0,50000000)) == 486
    assert len(refs_exons_by_chr["chr22"].find(0,50000000)) == 1824
    assert len(junctions_by_chr["chr22"]) == 4
    assert len(junctions_by_gene.keys()) == 913
    assert len(start_ends_by_gene) == 1388

def test_reference_parser_correct1exon(reference_parser_input):
    refs_1exon_by_chr, _, _, _, _ = reference_parser(*list(reference_parser_input.values()))
    # Check if the output is correct
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].id == "ENST00000615943.1"
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].chrom == "chr22"
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].strand == "-"
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].txStart == 10736170
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].txEnd == 10736283
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].cdsStart == 10736283
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].cdsEnd == 10736283
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].exonCount == 1
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].exonStarts == [10736170]
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].exonEnds == [10736283]
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].gene == "ENSG00000277248.1"
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].length == 113
    assert refs_1exon_by_chr["chr22"].find(0,50000000)[0].junctions == []

def test_reference_parser_correctExons(reference_parser_input):
    _, refs_exons_by_chr, _, _, _ = reference_parser(*list(reference_parser_input.values()))
    # Check if the output is correct
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].id == "ENST00000635667.1"
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].chrom == "chr22"
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].strand == "-"
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].txStart == 10939387
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].txEnd == 10961338
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].cdsStart == 10961338
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].cdsEnd == 10961338
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].exonCount == 9
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].exonStarts == [10939387, 10940596, 10941690, 10944966, 10947303, 10949211, 10950048, 10959066, 10961282]
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].exonEnds == [10939423, 10940707, 10941780, 10945053, 10947418, 10949269, 10950174, 10959136, 10961338] 
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].gene == "ENSG00000283047.1"
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].length == 749
    assert refs_exons_by_chr["chr22"].find(0,50000000)[0].junctions == [(10939423, 10940596), (10940707, 10941690), (10941780, 10944966), (10945053, 10947303), (10947418, 10949211), (10949269, 10950048), (10950174, 10959066), (10959136, 10961282)]

def test_reference_parser_correctJunctionsChr(reference_parser_input):
    _, _, junctions_by_chr, _, _ = reference_parser(*list(reference_parser_input.values()))
    assert len(junctions_by_chr["chr22"]) == 4
    assert len(junctions_by_chr["chr22"]["donors"]) == 5550
    assert junctions_by_chr["chr22"]["donors"][0] == 10939423
    assert len(junctions_by_chr["chr22"]["acceptors"]) == 5576
    assert junctions_by_chr["chr22"]["acceptors"][0] == 10940596
    assert junctions_by_chr["chr22"]["da_pairs"]['+'][0] == (11066515, 11067984)
    assert junctions_by_chr["chr22"]["da_pairs"]['-'][0] == (10939423, 10940596)

def test_reference_parser_correctJunctionsGene(reference_parser_input):
    _, _, _, junctions_by_gene, _ = reference_parser(*list(reference_parser_input.values()))
    assert len(junctions_by_gene.keys()) == 913
    assert len(junctions_by_gene["ENSG00000206195.11"]) == 11
    assert junctions_by_gene["ENSG00000206195.11"].pop() == (15791152, 15815475)

def test_reference_parser_correctStartEnds(reference_parser_input):
    _, _, _, _, start_ends_by_gene = reference_parser(*list(reference_parser_input.values()))
    assert len(start_ends_by_gene) == 1388
    assert start_ends_by_gene["ENSG00000206195.11"]["begin"] == {15784958, 15784962, 15784976, 15784991, 15787699}
    assert start_ends_by_gene["ENSG00000206195.11"]["end"] == {15829984, 15827434, 15790573, 15791387, 15827708}


### isoforms_parser ###

@pytest.fixture
def input_file():
    return os.path.join(main_path, "test/test_data/test_isoforms.genePred")

def count_lines(filename):
    with open(filename) as f:
        return sum(1 for line in f)

def test_isoforms_parser(input_file):
    isoforms_by_chr = isoforms_parser(input_file)
    assert len(isoforms_by_chr.keys()) == 1
    assert len(isoforms_by_chr["chr22"]) == count_lines(input_file)

def test_isoforms_parser_sorted(input_file):
    isoforms_by_chr = isoforms_parser(input_file)
    # Check if the isoforms are sorted by txStart
    txStarts = [isoform.txStart for isoform in isoforms_by_chr["chr22"]]
    assert txStarts == sorted(txStarts)


### parse_corrORF ###

@pytest.fixture
def corrORF_file():
    return os.path.join(main_path, "test/test_data/corrected_ORF_test.fasta")

@pytest.fixture
def gmst_rex():
    return re.compile(r'(\S+\t\S+\|GeneMark.hmm)\|(\d+)_aa\|(\S)\|(\d+)\|(\d+)')

def test_parse_corrORF(corrORF_file,gmst_rex):
    corrORF = parse_corrORF(corrORF_file,gmst_rex)
    assert len(corrORF) == 3
    assert corrORF["PB.124736.1"].orf_length == 210
    assert corrORF["PB.124736.1"].cds_start == 2305
    assert corrORF["PB.124736.1"].cds_end == 2937
    assert corrORF["PB.124736.1"].orf_seq == "MSRAGSHPKPAIPGRGGEKLALLLAMLGGRCCPAGGPESHTPRKPWPPIRRERPGLRSPKPREEAGPRADSKAASPALFMGPGSRGRRRPPSWARRRARGGRGSCPEPHVRAGPRFFIALAAQAARGLRVSSQGPGLWRGPAARKEGVARVGVAAGLRLRPRGRGRKTTPALSGHPARFPLQPGDSRERSRNHRALELTWQPRGAKAGGA"
    assert corrORF["PB.124736.1"].proteinID == "PB.124736.1"
    
def test_parse_corrORF_empty_file(gmst_rex):
    with patch('builtins.open', mock_open(read_data="")):
        result = parse_corrORF('dummy_file.fasta', gmst_rex)
    
    assert len(result) == 0

def test_parse_corrORF_invalid_format(gmst_rex):
    invalid_fasta = """>Invalid_format
MDEGTYIHALNNGLFTLGAPHKEVDEGPSPPEQFTAVKLSDSRITLKSGYGKYLGINSDELVVGHSDAIGPREQWEPVFKNGKMAFSASNSRFIRCSAKSKTAGEEEMIKIRSCAERETKEKDDIPEEDKGNIKQCEI"""
    
    with patch('builtins.open', mock_open(read_data=invalid_fasta)):
        with pytest.raises(SystemExit):
            parse_corrORF('dummy_file.fasta', gmst_rex)

def test_parse_corrORF_missing_fields(gmst_rex):
    missing_fields_fasta = """>PB.83093.1\tgene_1|GeneMark.hmm|138_aa|+|575
MDEGTYIHALNNGLFTLGAPHKEVDEGPSPPEQFTAVKLSDSRITLKSGYGKYLGINSDELVVGHSDAIGPREQWEPVFKNGKMAFSASNSRFIRCSAKSKTAGEEEMIKIRSCAERETKEKDDIPEEDKGNIKQCEI"""
    
    with patch('builtins.open', mock_open(read_data=missing_fields_fasta)):
        with pytest.raises(SystemExit):
            parse_corrORF('dummy_file.fasta', gmst_rex)

def test_parse_corrORF_non_integer_fields(gmst_rex):
    non_integer_fasta = """>PB.83093.1\tgene_1|GeneMark.hmm|138_aa|+|575|not_an_integer
MDEGTYIHALNNGLFTLGAPHKEVDEGPSPPEQFTAVKLSDSRITLKSGYGKYLGINSDELVVGHSDAIGPREQWEPVFKNGKMAFSASNSRFIRCSAKSKTAGEEEMIKIRSCAERETKEKDDIPEEDKGNIKQCEI"""
    
    with patch('builtins.open', mock_open(read_data=non_integer_fasta)):
        with pytest.raises(SystemExit):
            parse_corrORF('dummy_file.fasta', gmst_rex)

### parse_GMST ###

@pytest.fixture
def gmst_file():
    return os.path.join(main_path, "test/test_data/GMST_test")

@pytest.fixture
def corrORF_gmst_file():
    return os.path.join(main_path, "test/test_data/corrORF_gmst_test.fasta")

def test_parse_GMST_goodORF(corrORF_gmst_file, corrORF_file, gmst_rex, gmst_file):
    # Ensure the output file doesn't exist before running the function
    assert os.path.exists(corrORF_gmst_file) == False, f"File {corrORF_gmst_file} already exists. Abort!"
    
    orfDict = parse_GMST(corrORF_gmst_file, gmst_rex, gmst_file)
    
    # Check if the correct number of ORFs were parsed
    assert len(orfDict) == 3, "Expected 3 ORFs, but got a different number"

    # Check properties of a good ORF (PB.124736.1)
    assert orfDict["PB.124736.1"].orf_length == 210, "ORF length mismatch for PB.124736.1"
    assert orfDict["PB.124736.1"].cds_start == 2305, "CDS start position mismatch for PB.124736.1"
    assert orfDict["PB.124736.1"].cds_end == 2937, "CDS end position mismatch for PB.124736.1"
    assert orfDict["PB.124736.1"].orf_seq == "MSRAGSHPKPAIPGRGGEKLALLLAMLGGRCCPAGGPESHTPRKPWPPIRRERPGLRSPKPREEAGPRADSKAASPALFMGPGSRGRRRPPSWARRRARGGRGSCPEPHVRAGPRFFIALAAQAARGLRVSSQGPGLWRGPAARKEGVARVGVAAGLRLRPRGRGRKTTPALSGHPARFPLQPGDSRERSRNHRALELTWQPRGAKAGGA", "ORF sequence mismatch for PB.124736.1"
    assert orfDict["PB.124736.1"].proteinID == "PB.124736.1", "Protein ID mismatch for PB.124736.1"

    # Check properties of a fixed ORF (PB.96068.1)
    assert orfDict["PB.96068.1"].orf_length == 111, "ORF length mismatch for PB.96068.1"
    assert orfDict["PB.96068.1"].cds_start == 40, "CDS start position mismatch for PB.96068.1"
    assert orfDict["PB.96068.1"].cds_end == 375, "CDS end position mismatch for PB.96068.1"
    assert orfDict["PB.96068.1"].orf_seq == "MSLRVGARAKRNPWASGDPGGPDQCPLVVGADAWAHCGRAGPEVQVPAVDPGGGWENRRGVPAVKRILEAQEQLCFQCPLGVSKSNKKRINLWVPQKSPIFKSSVYESTDS", "ORF sequence mismatch for PB.96068.1"
    assert orfDict["PB.96068.1"].proteinID == "PB.96068.1", "Protein ID mismatch for PB.96068.1"

    # Compare output file with expected file
    expected_records = list(SeqIO.parse(corrORF_file, 'fasta'))
    actual_records = list(SeqIO.parse(corrORF_gmst_file, 'fasta'))
    
    # Check if the number of records in both files match
    assert len(expected_records) == len(actual_records), "Number of records in output file doesn't match expected file"
    
    # Compare each record in the output file with the expected file
    for expected, actual in zip(expected_records, actual_records):
        assert str(expected.seq) == str(actual.seq), f"Sequence mismatch for record {expected.id}"
        assert expected.id == actual.id, f"ID mismatch for record {expected.id}"
        assert expected.description == actual.description, f"Description mismatch for record {expected.id}"

    # Clean up: remove the output file after the test
    os.remove(corrORF_gmst_file)

