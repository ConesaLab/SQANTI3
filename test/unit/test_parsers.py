import logging
import re
import sys,os,pytest
from unittest.mock import mock_open, patch
from Bio import SeqIO

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)
from src.parsers import isoforms_parser, parse_corrORF, parse_TD2, reference_parser, FLcount_parser
from bx.intervals.intersection import IntervalTree

class TestReferenceParser:
    @pytest.fixture
    def reference_parser_input(self):
        data = {
            "gtf": os.path.join(main_path, "test/test_data/reference/test_reference.gtf"),
            "outdir": os.path.join(main_path, "test/test_data"),
            "prefix": "test",
            "genome_dict": {"chr22": "Sequence"},
            "gene_name": False,
            "isoAnnot": False
        }
        return data
    
    @pytest.fixture
    def tester_logger(self):
        # Define a logger for testing
        logger = logging.getLogger("tester_logger")
        logger.setLevel(logging.INFO)
        # Return the logger
        return logger

    def test_reference_parser_genePred(self, reference_parser_input):
        # Test if the genePred file is created
        genePred_annot = os.path.join(main_path,reference_parser_input["outdir"],
                                      f"refAnnotation_{reference_parser_input['prefix']}.genePred")
        if os.path.exists(genePred_annot):
            print(f"Removing {genePred_annot}")
            os.remove(genePred_annot)
        reference_parser(*list(reference_parser_input.values()))
        assert os.path.exists(genePred_annot)

    def test_reference_parser_genePred_found(self, reference_parser_input, caplog, tester_logger):
        # Test if the genePred file is detected
        genePred_annot = os.path.join(main_path, reference_parser_input["outdir"],
                                       f"refAnnotation_{reference_parser_input['prefix']}.genePred")
        reference_parser(*list(reference_parser_input.values()), logger=tester_logger)
        # Capture the printed output    # Check if the specific print statement is in stdout
        assert f"{genePred_annot} already exists. Using it." in caplog.text
        assert os.path.exists(genePred_annot)

    def test_reference_parser_notInGenome(self, reference_parser_input, caplog, tester_logger):
        reference_parser_input["genome_dict"] = {"chr21": "Sequence"}
        reference_parser(*list(reference_parser_input.values()), logger=tester_logger)
        # Check for the warning message in stderr
        expected_warning = "Reference annotation contains chromosomes not in genome:"
        caplog.set_level(logging.WARNING)
        assert expected_warning in caplog.text

    def test_reference_parser_correctOutput_length(self, reference_parser_input):
        refs_1exon_by_chr, refs_exons_by_chr, junctions_by_chr, junctions_by_gene, start_ends_by_gene = reference_parser(*list(reference_parser_input.values()))
        # Check if the output is correct
        assert len(refs_1exon_by_chr.keys()) == 1
        assert len(refs_1exon_by_chr["chr22"].find(0,50000000)) == 486
        assert len(refs_exons_by_chr["chr22"].find(0,50000000)) == 1824
        assert len(junctions_by_chr["chr22"]) == 4
        assert len(junctions_by_gene.keys()) == 913
        assert len(start_ends_by_gene) == 1388

    def test_reference_parser_correct1exon(self, reference_parser_input):
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

    def test_reference_parser_correctExons(self, reference_parser_input):
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

    def test_reference_parser_correctJunctionsChr(self, reference_parser_input):
        _, _, junctions_by_chr, _, _ = reference_parser(*list(reference_parser_input.values()))
        assert len(junctions_by_chr["chr22"]) == 4
        assert len(junctions_by_chr["chr22"]["donors"]) == 5550
        assert junctions_by_chr["chr22"]["donors"][0] == 10939423
        assert len(junctions_by_chr["chr22"]["acceptors"]) == 5576
        assert junctions_by_chr["chr22"]["acceptors"][0] == 10940596
        assert junctions_by_chr["chr22"]["da_pairs"]['+'][0] == (11066515, 11067984)
        assert junctions_by_chr["chr22"]["da_pairs"]['-'][0] == (10939423, 10940596)

    def test_reference_parserc_correctJunctionsGene(self, reference_parser_input):
        _, _, _, junctions_by_gene, _ = reference_parser(*list(reference_parser_input.values()))
        assert len(junctions_by_gene.keys()) == 913
        assert len(junctions_by_gene["ENSG00000206195.11"]) == 11
        assert (15791152, 15815475) in junctions_by_gene["ENSG00000206195.11"]

    def test_reference_parser_correctStartEnds(self, reference_parser_input):
        _, _, _, _, start_ends_by_gene = reference_parser(*list(reference_parser_input.values()))
        assert len(start_ends_by_gene) == 1388
        assert start_ends_by_gene["ENSG00000206195.11"]["begin"] == {15784958, 15784962, 15784976, 15784991, 15787699}
        assert start_ends_by_gene["ENSG00000206195.11"]["end"] == {15829984, 15827434, 15790573, 15791387, 15827708}


class TestIsoformsParser:
    @pytest.fixture
    def input_file(self):
        return os.path.join(main_path, "test/test_data/isoforms/test_isoforms.genePred")

    def count_lines(self, filename):
        with open(filename) as f:
            return sum(1 for line in f)

    def test_isoforms_parser(self, input_file):
        isoforms_by_chr = isoforms_parser(input_file)
        assert len(isoforms_by_chr.keys()) == 1
        assert len(isoforms_by_chr["chr22"]) == self.count_lines(input_file)

    def test_isoforms_parser_sorted(self, input_file):
        isoforms_by_chr = isoforms_parser(input_file)
        # Check if the isoforms are sorted by txStart
        txStarts = [isoform.txStart for isoform in isoforms_by_chr["chr22"]]
        assert txStarts == sorted(txStarts)


class TestParseCorrORF:
    @pytest.fixture
    def corrORF_file(self):
        return os.path.join(main_path, "test/test_data/orfs/corrected_ORF_test.fasta")

    def test_parse_corrORF(self, corrORF_file):
        corrORF = parse_corrORF(corrORF_file)
        assert len(corrORF) == 4
        assert corrORF["PB.96068.1"].cds_length == 336
        assert corrORF["PB.96068.1"].cds_start == 40
        assert corrORF["PB.96068.1"].cds_end == 375
        assert corrORF["PB.96068.1"].psauron_score == 0.546
        assert corrORF["PB.96068.1"].cds_type == "complete"
        assert corrORF["PB.96068.1"].protein_seq == "MSLRVGARAKRNPWASGDPGGPDQCPLVVGADAWAHCGRAGPEVQVPAVDPGGGWENRRGVPAVKRILEAQEQLCFQCPLGVSKSNKKRINLWVPQKSPIFKSSVYESTDS*"
        assert corrORF["PB.96068.1"].protein_length == 112
        assert corrORF["PB.96068.1"].proteinID == "PB.96068.1"
        assert len(corrORF["PB.124735.1"].protein_seq) == corrORF["PB.124735.1"].protein_length
        
    def test_parse_corrORF_empty_file(self):
        with patch('builtins.open', mock_open(read_data="")):
            result = parse_corrORF('dummy_file.fasta')
        
        assert len(result) == 0

    def test_parse_corrORF_invalid_format(self):
        invalid_fasta = """>Invalid_format
MDEGTYIHALNNGLFTLGAPHKEVDEGPSPPEQFTAVKLSDSRITLKSGYGKYLGINSDELVVGHSDAIGPREQWEPVFKNGKMAFSASNSRFIRCSAKSKTAGEEEMIKIRSCAERETKEKDDIPEEDKGNIKQCEI"""
        
        with patch('builtins.open', mock_open(read_data=invalid_fasta)):
            with pytest.raises(SystemExit):
                parse_corrORF('dummy_file.fasta')

    def test_parse_corrORF_missing_fields(self):
        missing_fields_fasta = """>PB.83093.1\tgene_1|GeneMark.hmm|138_aa|+|575
MDEGTYIHALNNGLFTLGAPHKEVDEGPSPPEQFTAVKLSDSRITLKSGYGKYLGINSDELVVGHSDAIGPREQWEPVFKNGKMAFSASNSRFIRCSAKSKTAGEEEMIKIRSCAERETKEKDDIPEEDKGNIKQCEI"""
        
        with patch('builtins.open', mock_open(read_data=missing_fields_fasta)):
            with pytest.raises(SystemExit):
                parse_corrORF('dummy_file.fasta')

    def test_parse_corrORF_non_integer_fields(self):
        non_integer_fasta = """>PB.83093.1\tgene_1|GeneMark.hmm|138_aa|+|575|not_an_integer
MDEGTYIHALNNGLFTLGAPHKEVDEGPSPPEQFTAVKLSDSRITLKSGYGKYLGINSDELVVGHSDAIGPREQWEPVFKNGKMAFSASNSRFIRCSAKSKTAGEEEMIKIRSCAERETKEKDDIPEEDKGNIKQCEI"""
        
        with patch('builtins.open', mock_open(read_data=non_integer_fasta)):
            with pytest.raises(SystemExit):
                parse_corrORF('dummy_file.fasta')

class TestParseTD2:
    @pytest.fixture
    def td2_file(self):
        return os.path.join(main_path, "test/test_data/orfs/TD2_test.faa")

    @pytest.fixture
    def corrORF_td2_file(self):
        return os.path.join(main_path, "test/test_data/orfs/corrORF_td2_test.fasta")

    @pytest.fixture
    def corrORF_file(self):
        return os.path.join(main_path, "test/test_data/orfs/corrected_ORF_test.fasta")

    def test_parse_TD2_goodORF(self, corrORF_td2_file, corrORF_file, td2_file):
        # Ensure the output file doesn't exist before running the function
        assert os.path.exists(corrORF_td2_file) == False, f"File {corrORF_td2_file} already exists. Abort!"
        
        orfDict = parse_TD2(corrORF_td2_file, td2_file)
        
        # Check if the correct number of ORFs were parsed
        assert len(orfDict) == 4, "Expected 4 ORFs, but got a different number"

        # Check properties of a good ORF (PB.124735.1)
        assert orfDict["PB.124735.1"].protein_length == 195, "ORF length mismatch for PB.124735.1"
        assert orfDict["PB.124735.1"].cds_start == 1312, "CDS start position mismatch for PB.124735.1"
        assert orfDict["PB.124735.1"].cds_end == 1896, "CDS end position mismatch for PB.124735.1"
        assert orfDict["PB.124735.1"].protein_seq == "MESCSVSQAGVQWHYLRSLQPAPSLFKRFSNLSLSRTWDHRCPPPCPANFLVFLVETVFLHVDQASLELVTSGDLPASASQTAGIAGMNHRTQPNVILFKNFFSFVFSFFSFLFFYFLFFLSLSFLSFFFLRRCLTLFPKLEHSGTISADCNLHLPSSSNSPASASQVAGTTGVCHYAQLIFVFLTEIEFYYLY*", "ORF sequence mismatch for PB.124735.1"
        assert orfDict["PB.124735.1"].proteinID == "PB.124735.1", "Protein ID mismatch for PB.124735.1"

        # Check properties of a fixed ORF (PB.96068.1)
        assert orfDict["PB.96068.1"].protein_length == 112, "ORF length mismatch for PB.96068.1"
        assert orfDict["PB.96068.1"].cds_start == 40, "CDS start position mismatch for PB.96068.1"
        assert orfDict["PB.96068.1"].cds_end == 375, "CDS end position mismatch for PB.96068.1"
        assert orfDict["PB.96068.1"].protein_seq == "MSLRVGARAKRNPWASGDPGGPDQCPLVVGADAWAHCGRAGPEVQVPAVDPGGGWENRRGVPAVKRILEAQEQLCFQCPLGVSKSNKKRINLWVPQKSPIFKSSVYESTDS*", "ORF sequence mismatch for PB.96068.1"
        assert orfDict["PB.96068.1"].proteinID == "PB.96068.1", "Protein ID mismatch for PB.96068.1"

        # check properties of non_PacBio
        assert orfDict["Non.standard12"].protein_length == 407, "ORF length mismatch for Non.standard12"
        assert orfDict["Non.standard12"].cds_start == 123, "CDS start position mismatch for Non.standard12"
        assert orfDict["Non.standard12"].cds_end == 1343, "CDS end position mismatch for Non.standard12"
        assert orfDict["Non.standard12"].protein_seq == "ASDFASFFADSFAADFFF*", "ORF sequence mismatch for Non.standard12"
        assert orfDict["Non.standard12"].proteinID == "Non.standard12", "Protein ID mismatch for Non.standard12"
        # Compare output file with expected file
        expected_records = list(SeqIO.parse(corrORF_file, 'fasta'))
        actual_records = list(SeqIO.parse(corrORF_td2_file, 'fasta'))
        
        # Check if the number of records in both files match
        assert len(expected_records) == len(actual_records), "Number of records in output file doesn't match expected file"
        
        # Compare each record in the output file with the expected file
        for expected, actual in zip(expected_records, actual_records):
            assert str(expected.seq) == str(actual.seq), f"Sequence mismatch for record {expected.id}"
            assert expected.id == actual.id, f"ID mismatch for record {expected.id}"
            # Eliminated because one has a tab and the other has a space
            #assert expected.description == actual.description, f"Description mismatch for record {expected.id}"

        # Clean up: remove the output file after the test
        os.remove(corrORF_td2_file)


class TestFLcountParser:
    @pytest.fixture
    def fl_count_single_tsv(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_single_sample.tsv")

    @pytest.fixture
    def fl_count_single_csv(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_single_sample.csv")

    @pytest.fixture
    def fl_count_multi_tsv(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_multi_sample.tsv")

    @pytest.fixture
    def fl_count_multi_csv(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_multi_sample.csv")

    @pytest.fixture
    def fl_count_with_comments(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_with_comments.tsv")

    @pytest.fixture
    def fl_count_with_na(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_with_na.tsv")

    @pytest.fixture
    def fl_count_mixed_numeric(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_mixed_numeric.tsv")

    @pytest.fixture
    def fl_count_empty(self):
        return os.path.join(main_path, "test/test_data/abundance/fl_count_empty.tsv")

    def test_FLcount_parser_single_sample_tsv(self, fl_count_single_tsv):
        samples, fl_count_dict = FLcount_parser(fl_count_single_tsv)
        
        # Single sample returns flat dictionary
        assert len(samples) == 1
        assert samples[0] == "count_fl"
        assert len(fl_count_dict) == 5
        assert fl_count_dict["PB.124830.1"] == 150
        assert fl_count_dict["PB.103714.1"] == 200
        assert fl_count_dict["PB.103724.1"] == 50
        assert fl_count_dict["PB.103781.1"] == 300
        assert fl_count_dict["PB.23068.2"] == 25

    def test_FLcount_parser_single_sample_csv(self, fl_count_single_csv):
        samples, fl_count_dict = FLcount_parser(fl_count_single_csv)
        
        assert len(samples) == 1
        assert samples[0] == "count_fl"
        assert len(fl_count_dict) == 5
        assert fl_count_dict["PB.124830.1"] == 150
        assert fl_count_dict["PB.103714.1"] == 200
        assert fl_count_dict["PB.103724.1"] == 50

    def test_FLcount_parser_multi_sample_tsv(self, fl_count_multi_tsv):
        samples, fl_count_dict = FLcount_parser(fl_count_multi_tsv)
        
        # Multi-sample returns nested dictionary
        assert len(samples) == 3
        assert set(samples) == {"sample1", "sample2", "sample3"}
        assert len(fl_count_dict) == 5
        
        # Check nested structure
        assert fl_count_dict["PB.124830.1"]["sample1"] == 150
        assert fl_count_dict["PB.124830.1"]["sample2"] == 200
        assert fl_count_dict["PB.124830.1"]["sample3"] == 180
        
        assert fl_count_dict["PB.103781.1"]["sample1"] == 300
        assert fl_count_dict["PB.103781.1"]["sample2"] == 280
        assert fl_count_dict["PB.103781.1"]["sample3"] == 310

    def test_FLcount_parser_multi_sample_csv(self, fl_count_multi_csv):
        samples, fl_count_dict = FLcount_parser(fl_count_multi_csv)
        
        assert len(samples) == 3
        assert set(samples) == {"sample1", "sample2", "sample3"}
        assert len(fl_count_dict) == 5
        
        assert fl_count_dict["PB.103714.1"]["sample1"] == 200
        assert fl_count_dict["PB.103714.1"]["sample2"] == 220
        assert fl_count_dict["PB.103714.1"]["sample3"] == 190

    def test_FLcount_parser_with_comments(self, fl_count_with_comments):
        samples, fl_count_dict = FLcount_parser(fl_count_with_comments)
        
        # Should skip comment lines and parse correctly
        assert len(samples) == 1
        assert len(fl_count_dict) == 3
        assert fl_count_dict["PB.124830.1"] == 150
        assert fl_count_dict["PB.103714.1"] == 200
        assert fl_count_dict["PB.103724.1"] == 50

    def test_FLcount_parser_with_na_values(self, fl_count_with_na):
        samples, fl_count_dict = FLcount_parser(fl_count_with_na)
        
        # NA values should be converted to 0
        assert len(samples) == 3
        assert fl_count_dict["PB.124830.1"]["sample1"] == 150
        assert fl_count_dict["PB.124830.1"]["sample2"] == 0  # NA converted to 0
        assert fl_count_dict["PB.124830.1"]["sample3"] == 180
        
        assert fl_count_dict["PB.103724.1"]["sample1"] == 0  # NA converted to 0
        assert fl_count_dict["PB.103724.1"]["sample2"] == 60
        assert fl_count_dict["PB.103724.1"]["sample3"] == 55

    def test_FLcount_parser_mixed_numeric_types(self, fl_count_mixed_numeric):
        _, fl_count_dict = FLcount_parser(fl_count_mixed_numeric)
        
        # Should handle both int and float values
        assert fl_count_dict["PB.124830.1"] == 150
        assert fl_count_dict["PB.103714.1"] == 200.5
        assert fl_count_dict["PB.103724.1"] == 50
        assert fl_count_dict["PB.103781.1"] == 300.75
        
        # Check types
        assert isinstance(fl_count_dict["PB.124830.1"], int)
        assert isinstance(fl_count_dict["PB.103714.1"], float)

    def test_FLcount_parser_empty_file(self, fl_count_empty):
        samples, fl_count_dict = FLcount_parser(fl_count_empty)
        
        # Empty file (only header) should return empty dict
        assert len(samples) == 1
        assert len(fl_count_dict) == 0

