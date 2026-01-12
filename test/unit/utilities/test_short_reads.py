import pytest, sys, os
import subprocess
import numpy as np
from unittest.mock import patch, mock_open,patch, MagicMock
from src.utilities.short_reads import star_mapping, get_TSS_bed, get_bam_header, get_ratio_TSS

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
@patch('subprocess.call')
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
    # Create a temporary directory structure
    bam_file = tmp_path / "test.bam"
    bam_file.touch()
    return str(bam_file)

def test_get_bam_header_file_exists(setup_test_environment_bam):
    bam_file = setup_test_environment_bam
    expected_output = os.path.dirname(bam_file) + "/chr_order.txt"
    
    # Create a mock chr_order.txt file
    with open(expected_output, 'w') as f:
        f.write("mock content")
    
    result = get_bam_header(bam_file)
    assert result == expected_output
    assert os.path.isfile(result)
    

@patch('subprocess.run')
def test_get_bam_header_file_not_exists(mock_subprocess, setup_test_environment_bam):
    bam_file = setup_test_environment_bam
    with pytest.raises(SystemExit):
        get_bam_header(bam_file)
    

@patch('subprocess.run')
def test_get_bam_header_subprocess_error(mock_subprocess, setup_test_environment_bam):
    bam_file = setup_test_environment_bam
    mock_subprocess.side_effect = subprocess.CalledProcessError(1, 'cmd')
    
    with pytest.raises(SystemExit):
        get_bam_header(bam_file)

def test_get_bam_header_invalid_bam(tmp_path):
    non_existent_bam = str(tmp_path / "non_existent.bam")
    
    with pytest.raises(FileNotFoundError):
        get_bam_header(non_existent_bam)


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

def test_get_ratio_TSS_mean(test_files):
    
    result = get_ratio_TSS(
        test_files["inside_bed"],
        test_files["outside_bed"],
        test_files["bam_files"],
        test_files["chr_order"],
        "mean"
    )

    assert isinstance(result, dict)
    assert len(result) > 0  
    assert abs(result["PB.103709.3"]["return_ratio"] - 1425.2574257425742 <= 0.01) # Ensure we got some results
    assert np.isnan(result["PB.103820.1"]["return_ratio"])
