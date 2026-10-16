import pytest, sys, os
import subprocess
import numpy as np
from unittest.mock import patch, mock_open,patch, MagicMock
import pysam
from src.utilities.short_reads import star_mapping, get_TSS_bed, get_bam_header, get_ratio_TSS, \
    get_bam_index, count_TSS_reads

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
sys.path.insert(0, main_path)

# Most of the functions in thsi module are jsut to run programs or create directories.
# There is no need for testing, as they wont be modified, and the tools they call already
# have their own tests.

### star_mapping ###
@pytest.fixture
def setup_test_environment(tmp_path):
    # Create a temporary directory structure
    index_dir = tmp_path / "index"
    output_dir = tmp_path / "output"
    mapping_dir = output_dir / "STAR_mapping"
    index_dir.mkdir()
    output_dir.mkdir()
    mapping_dir.mkdir()

    # Create a mock SR_fofn file
    sr_fofn_content = "sample1_R1.fastq.gz sample1_R2.fastq.gz\nsample2_R1.fastq\n"
    sr_fofn_file = tmp_path / "sr_fofn.txt"
    sr_fofn_file.write_text(sr_fofn_content)

    return {
        "index_dir": str(index_dir),
        "output_dir": str(output_dir),
        "sr_fofn": str(sr_fofn_file),
        "mapping_dir": str(mapping_dir)
    }

# TODO: Change this test to use the star_cmd function. Mapping per se will not be tested
# @patch('subprocess.call')
# def test_star_mapping(mock_subprocess_call, setup_test_environment):
#     env = setup_test_environment
#     # Run the function
#     star_cmd(env["index_dir"], env["sr_fofn"], env["output_dir"], 4)

#     # Assertions
#     assert mock_subprocess_call.call_count == 2

#     # Check calls for compressed files (sample1)
#     compressed_call_args = mock_subprocess_call.call_args_list[0][0][0]
#     assert compressed_call_args[0] == 'STAR'
#     assert '--runThreadN' in compressed_call_args
#     assert '--genomeDir' in compressed_call_args
#     assert '--readFilesIn' in compressed_call_args
#     assert 'sample1_R1.fastq.gz' in compressed_call_args
#     assert 'sample1_R2.fastq.gz' in compressed_call_args
#     assert '--readFilesCommand' in compressed_call_args
#     assert 'zcat' in compressed_call_args

#     # Check calls for uncompressed files (sample2)
#     uncompressed_call_args = mock_subprocess_call.call_args_list[1][0][0]
#     assert uncompressed_call_args[0] == 'STAR'
#     assert '--runThreadN' in uncompressed_call_args
#     assert '--genomeDir' in uncompressed_call_args
#     assert '--readFilesIn' in uncompressed_call_args
#     assert 'sample2_R1.fastq' in uncompressed_call_args
#     assert '--readFilesCommand' not in uncompressed_call_args

#     # Check if output directories are created
#     assert os.path.exists(env["mapping_dir"])

### get_TSS_bed ###

@pytest.fixture
def setup_test_environment_bed(tmp_path):
    # Create a temporary directory structure
    out_directory = tmp_path / "output"
    out_directory.mkdir()

    # Create a mock corrected_gtf file
    corrected_gtf = out_directory / "corrected.gtf"
    corrected_gtf.touch()

    # Create a mock chr_order file
    chr_order = out_directory / "chr_order.txt"
    chr_order.touch()

    # Create mock temporary files
    (out_directory / "coverage_inside_TSS.bed_tmp").touch()
    (out_directory / "coverage_outside_TSS.bed_tmp").touch()

    return {
        "out_directory": str(out_directory),
        "corrected_gtf": str(corrected_gtf),
        "chr_order": str(chr_order)
    }


@patch('BCBio.GFF.parse')
@patch('pybedtools.BedTool')
def test_get_TSS_bed(mock_bedtool, mock_bcbio_parse, setup_test_environment_bed):
    env = setup_test_environment_bed

    # Mock BCBio.GFF.parse to return a list with one record
    mock_feature = MagicMock()
    mock_feature.qualifiers = {"transcript_id": ["ENST00000000001"]}
    mock_feature.location = MagicMock()
    mock_feature.location.__str__.return_value = "[1:1000](+)"

    mock_record = MagicMock()
    mock_record.id = "chr1"
    mock_record.features = [mock_feature]

    mock_bcbio_parse.return_value = [mock_record]

    # Mock pybedtools.BedTool
    mock_bedtool_instance = MagicMock()
    mock_bedtool.return_value = mock_bedtool_instance

    # Run the function
    with patch('builtins.open', new_callable=mock_open) as mock_file:
        with patch('os.remove') as mock_remove:
            inside_sorted, outside_sorted = get_TSS_bed(env["corrected_gtf"], env["chr_order"])

    # Assertions
    assert inside_sorted == os.path.join(env["out_directory"], "inside_TSS.bed")
    assert outside_sorted == os.path.join(env["out_directory"], "outside_TSS.bed")

    # Check if files were opened correctly
    mock_file.assert_any_call(env["corrected_gtf"])
    mock_file.assert_any_call(os.path.join(env["out_directory"], 'coverage_inside_TSS.bed_tmp'), 'w')
    mock_file.assert_any_call(os.path.join(env["out_directory"], 'coverage_outside_TSS.bed_tmp'), 'w')

    # Check if BCBio.GFF.parse was called correctly
    mock_bcbio_parse.assert_called_once()

    # Check if pybedtools.BedTool was called correctly
    assert mock_bedtool.call_count == 2

    # Check if sort method was called on BedTool instances
    assert mock_bedtool_instance.sort.call_count == 2

    # Check if os.remove was called for temporary files
    assert mock_remove.call_count == 2
    mock_remove.assert_any_call(os.path.join(env["out_directory"], 'coverage_inside_TSS.bed_tmp'))
    mock_remove.assert_any_call(os.path.join(env["out_directory"], 'coverage_outside_TSS.bed_tmp'))

@patch('BCBio.GFF')
@patch('pybedtools.BedTool')
def test_get_TSS_bed_negative_strand(mock_bedtool, mock_bcbio_parse, setup_test_environment_bed):
    env = setup_test_environment_bed

    # Mock BCBio_GFF.parse to return a record with negative strand
    mock_feature = MagicMock()
    mock_feature.qualifiers = {"transcript_id": ["ENST00000000002"]}
    mock_feature.location = MagicMock()
    mock_feature.location.__str__.return_value = "[1000:2000](-)"

    mock_record = MagicMock()
    mock_record.id = "chr1"
    mock_record.features = [mock_feature]

    mock_bcbio_parse.return_value = [mock_record]

    # Run the function
    with patch('builtins.open', new_callable=mock_open) as mock_file:
        inside_sorted, outside_sorted = get_TSS_bed(env["corrected_gtf"], env["chr_order"])

    # Assertions (similar to the previous test)
    assert inside_sorted == os.path.join(env["out_directory"], "inside_TSS.bed")
    assert outside_sorted == os.path.join(env["out_directory"], "outside_TSS.bed")

    # Additional checks specific to negative strand processing can be added here

@patch('os.remove')
@patch('BCBio.GFF')
@patch('pybedtools.BedTool')
def test_get_TSS_bed_file_cleanup(mock_bedtool, mock_bcbio_parse, mock_remove, setup_test_environment_bed):
    env = setup_test_environment_bed

    # Mock BCBio_GFF.parse (minimal setup)
    mock_bcbio_parse.return_value = []

    # Run the function
    with patch('builtins.open', new_callable=mock_open):
        get_TSS_bed(env["corrected_gtf"], env["chr_order"])

    # Check if temporary files were removed
    assert mock_remove.call_count == 2
    mock_remove.assert_any_call(os.path.join(env["out_directory"], 'coverage_inside_TSS.bed_tmp'))
    mock_remove.assert_any_call(os.path.join(env["out_directory"], 'coverage_outside_TSS.bed_tmp'))

### get_bam_header ###

@pytest.fixture
def setup_test_environment_bam(tmp_path):
    # Input BAM and output directory live in different places, as with --SR_bam
    input_dir = tmp_path / "input"
    out_dir = tmp_path / "output"
    input_dir.mkdir()
    out_dir.mkdir()
    bam_file = input_dir / "test.bam"
    bam_file.touch()
    return str(bam_file), str(out_dir)

def test_get_bam_header_file_exists(setup_test_environment_bam):
    bam_file, out_dir = setup_test_environment_bam
    expected_output = os.path.join(out_dir, "chr_order.txt")
    
    # Create a mock chr_order.txt file
    with open(expected_output, 'w') as f:
        f.write("mock content")
    
    result = get_bam_header(bam_file, out_dir)
    assert result == expected_output
    assert os.path.isfile(result)

@patch('src.utilities.short_reads.run_command')
def test_get_bam_header_writes_inside_outdir(mock_run_command, setup_test_environment_bam):
    # Regression test for #611: nothing must be written next to the input BAM or its parent
    bam_file, out_dir = setup_test_environment_bam
    result = get_bam_header(bam_file, out_dir)
    assert result == os.path.join(out_dir, "chr_order.txt")
    log_file = mock_run_command.call_args[0][2]
    assert log_file == os.path.join(out_dir, "logs", "samtools_header.log")

@patch('subprocess.run')
def test_get_bam_header_file_not_exists(mock_subprocess, setup_test_environment_bam):
    bam_file, out_dir = setup_test_environment_bam
    with pytest.raises(SystemExit):
        get_bam_header(bam_file, out_dir)
    

@patch('subprocess.run')
def test_get_bam_header_subprocess_error(mock_subprocess, setup_test_environment_bam):
    bam_file, out_dir = setup_test_environment_bam
    mock_subprocess.side_effect = subprocess.CalledProcessError(1, 'cmd')
    
    with pytest.raises(SystemExit):
        get_bam_header(bam_file, out_dir)

def test_get_bam_header_invalid_bam(tmp_path):
    non_existent_bam = str(tmp_path / "non_existent.bam")
    
    with pytest.raises(FileNotFoundError):
        get_bam_header(non_existent_bam, str(tmp_path))


### get_ratio_TSS ###

@pytest.fixture
def test_files():
    chr_order = os.path.join(main_path, "test/test_data/other/chr_order.txt")
    inside_bed, outside_bed = get_TSS_bed(os.path.join(main_path,"test/test_data/isoforms/test_isoforms.gtf"),
                                          chr_order)
    return {
        "bam_files": [
            os.path.join(main_path, "test/test_data/bam/Rep1_test.bam"),
            os.path.join(main_path, "test/test_data/bam/Rep2_test.bam")
        ],
        "inside_bed": inside_bed,
        "outside_bed": outside_bed,
        "chr_order": chr_order
    }

# PB.103709.3 has ratios 49.51 (Rep1) and 2801.0 (Rep2); PB.137296.2 only has Rep1 (9.91, Rep2 is NaN)
@pytest.mark.parametrize("metric, expected", [
    ("max", 2801.0),
    ("mean", 1425.2574257425742),
    ("median", 1425.2574257425742),
    ("3quartile", 2113.128712871287),
])
def test_get_ratio_TSS_metrics(test_files, tmp_path, metric, expected):
    result = get_ratio_TSS(
        test_files["inside_bed"],
        test_files["outside_bed"],
        test_files["bam_files"],
        metric,
        4,
        str(tmp_path / "bam_index")
    )

    assert isinstance(result, dict)
    assert len(result) > 0
    assert result["PB.103709.3"]["return_ratio"] == pytest.approx(expected)
    # NaN replicates are skipped, and isoforms without coverage in any replicate stay NaN
    assert result["PB.137296.2"]["return_ratio"] == pytest.approx(9.910891, abs=1e-6)
    assert np.isnan(result["PB.103820.1"]["return_ratio"])


### get_bam_index and count_TSS_reads ###

def write_bam(bam_path, reads, sort_order="coordinate"):
    """Writes a BAM on chr1 from (name, flag, start, cigar) tuples. A cigar of None is an unmapped read."""
    header = {"HD": {"VN": "1.6", "SO": sort_order}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for name, flag, start, cigar in reads:
            read = pysam.AlignedSegment(bam.header)
            read.query_name = name
            read.flag = flag
            read.reference_id = 0
            read.reference_start = start
            if cigar is not None:
                read.cigarstring = cigar
                read.mapping_quality = 60
            read.query_sequence = "A" * 50
            bam.write(read)
    return str(bam_path)

@pytest.fixture
def tss_bam(tmp_path):
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    reads = [
        ("spliced", 0, 500, "20M1000N30M"),  # 500-1550: its intron covers the window
        ("primary", 0, 1020, "50M"),
        ("secondary", 256, 1050, "50M"),
        ("duplicate", 1024, 1060, "50M"),
        ("unmapped", 4, 1070, None),         # unmapped mate placed next to its pair
        ("far", 0, 2000, "50M"),
    ]
    return write_bam(input_dir / "Aligned.sortedByCoord.out.bam", reads)

class TestGetBamIndex:
    """Test suite for get_bam_index."""

    def test_uses_index_next_to_bam(self, tss_bam, tmp_path):
        pysam.index(tss_bam)
        assert get_bam_index(tss_bam, str(tmp_path / "bam_index"), 1) == f"{tss_bam}.bai"
        assert not os.path.exists(tmp_path / "bam_index")

    def test_creates_index_in_index_dir(self, tss_bam, tmp_path):
        """The input folder may be read-only, so the new index goes to index_dir (issue #611)."""
        index_dir = str(tmp_path / "bam_index")
        index = get_bam_index(tss_bam, index_dir, 2)
        assert os.path.dirname(index) == index_dir
        assert os.path.isfile(index)
        assert os.listdir(os.path.dirname(tss_bam)) == ["Aligned.sortedByCoord.out.bam"]

    def test_replicates_with_the_same_name_get_their_own_index(self, tmp_path):
        """STAR names every replicate Aligned.sortedByCoord.out.bam, in a folder per sample."""
        index_dir = str(tmp_path / "bam_index")
        indexes = []
        for sample in ["rep1", "rep2"]:
            (tmp_path / sample).mkdir()
            bam = write_bam(tmp_path / sample / "Aligned.sortedByCoord.out.bam", [("primary", 0, 1020, "50M")])
            indexes.append(get_bam_index(bam, index_dir, 1))
        assert indexes[0] != indexes[1]

    def test_index_older_than_bam_is_made_again(self, tss_bam, tmp_path):
        index_dir = str(tmp_path / "bam_index")
        index = get_bam_index(tss_bam, index_dir, 1)
        os.utime(index, (0, 0))
        assert get_bam_index(tss_bam, index_dir, 1) == index
        assert os.path.getmtime(index) >= os.path.getmtime(tss_bam)

    @patch("src.utilities.short_reads.qc_logger")
    def test_unsorted_bam(self, mock_logger, tmp_path):
        bam = write_bam(tmp_path / "unsorted.bam", [("b", 0, 2000, "50M"), ("a", 0, 100, "50M")], sort_order="unsorted")
        with pytest.raises(SystemExit):
            get_bam_index(bam, str(tmp_path / "bam_index"), 1)
        assert "must be sorted by coordinate" in mock_logger.error.call_args[0][0]
        assert not os.path.exists(tmp_path / "bam_index" / (os.path.abspath(bam).strip("/").replace("/", "_") + ".bai"))

class TestCountTSSReads:
    """Test suite for count_TSS_reads."""

    windows = [("chr1", 1000, 1100), ("chr1", 1900, 2000), ("chr1", 2000, 2010), ("chrUn", 0, 100)]

    def test_counts_as_bedtools(self, tss_bam, tmp_path):
        """Same counts as bedtools coverage -counts: reads spanning the window with an intron, secondary
        and duplicate reads count, unmapped reads don't, and a contig missing from the BAM gives 0."""
        bam_index = get_bam_index(tss_bam, str(tmp_path / "bam_index"), 1)
        assert count_TSS_reads(tss_bam, bam_index, self.windows, 1) == [4, 0, 1, 0]

    @pytest.mark.parametrize("cpus", [2, 3, 8])
    def test_same_counts_in_parallel(self, tss_bam, tmp_path, cpus):
        """Windows are dealt round-robin to the processes; the counts must come back in the input order."""
        bam_index = get_bam_index(tss_bam, str(tmp_path / "bam_index"), 1)
        windows = self.windows * 5
        assert count_TSS_reads(tss_bam, bam_index, windows, cpus) == [4, 0, 1, 0] * 5
