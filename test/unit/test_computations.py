import pytest
import os, sys
from unittest.mock import Mock, patch, MagicMock

# Mock problematic imports before importing qc_computations
sys.modules['pybedtools'] = MagicMock()
sys.modules['src.utilities.short_reads'] = MagicMock()

main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.append(main_path)
data_path = os.path.join(main_path, 'test','test_data')

from src.qc_computations import full_length_quantification


### Fixtures ###

@pytest.fixture
def mock_qc_logger():
    """Mock the qc_logger to verify logging calls without output"""
    with patch('src.qc_computations.qc_logger') as mock_logger:
        yield mock_logger


class MockIsoform:
    """Mock myQueryTranscripts object with minimal required attributes"""
    def __init__(self, isoform_id):
        self.id = isoform_id
        self.FL = 'NA'
        self.FL_dict = {}


@pytest.fixture
def create_isoforms_info():
    """Fixture factory to create isoforms_info dict with specified isoform IDs"""
    def _create(isoform_ids):
        return {iso_id: MockIsoform(iso_id) for iso_id in isoform_ids}
    return _create


@pytest.fixture
def fl_count_single_sample():
    """Path to single sample FL count test file"""
    return os.path.join(data_path, 'FL_count_single_sample.txt')


@pytest.fixture
def fl_count_multi_chain():
    """Path to multi-sample chain format FL count test file"""
    return os.path.join(data_path, 'FL_count_multi_chain.txt')


@pytest.fixture
def fl_count_multi_demux():
    """Path to multi-sample demux format FL count test file"""
    return os.path.join(data_path, 'FL_count_multi_demux.txt')


@pytest.fixture
def fl_count_with_na():
    """Path to FL count file with NA values"""
    return os.path.join(data_path, 'FL_count_with_NA.txt')


@pytest.fixture
def fl_count_float():
    """Path to FL count file with float values"""
    return os.path.join(data_path, 'FL_count_float_values.txt')


### Tests for full_length_quantification ###

def test_full_length_quantification_single_sample_perfect_match(
    fl_count_single_sample, create_isoforms_info, mock_qc_logger
):
    """Test single sample FL count with all isoforms matching"""
    # Create isoforms_info with the same isoforms as in the FL file
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1', 'PB.103709.3']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform', 'length', 'exons']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_single_sample, isoforms_info, fields
    )
    
    # Verify return values are same objects (in-place modification)
    assert result_info is isoforms_info
    assert result_fields == fields  # Unchanged for single sample
    
    # Verify FL values were set correctly
    assert isoforms_info['PB.103698.1'].FL == 10
    assert isoforms_info['PB.103704.2'].FL == 25
    assert isoforms_info['PB.103705.1'].FL == 5
    assert isoforms_info['PB.103707.1'].FL == 100
    assert isoforms_info['PB.103709.3'].FL == 0
    
    # Verify logger info calls
    assert mock_qc_logger.info.call_count >= 2
    mock_qc_logger.info.assert_any_call("**** Reading Full-length read abundance files.")
    mock_qc_logger.info.assert_any_call("Single-sample PacBio FL count format detected.")
    
    # Verify no warnings for perfect match
    assert mock_qc_logger.warning.call_count == 0


def test_full_length_quantification_multi_chain_perfect_match(
    fl_count_multi_chain, create_isoforms_info, mock_qc_logger
):
    """Test multi-sample chain format FL count with all isoforms matching"""
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1', 'PB.103709.3']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform', 'length']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_multi_chain, isoforms_info, fields
    )
    
    # Verify FL_dict populated correctly
    assert isoforms_info['PB.103698.1'].FL_dict == {'sample1': 10, 'sample2': 15}
    assert isoforms_info['PB.103704.2'].FL_dict == {'sample1': 25, 'sample2': 30}
    assert isoforms_info['PB.103705.1'].FL_dict == {'sample1': 5, 'sample2': 0}
    assert isoforms_info['PB.103707.1'].FL_dict == {'sample1': 100, 'sample2': 50}
    assert isoforms_info['PB.103709.3'].FL_dict == {'sample1': 0, 'sample2': 20}
    
    # Verify FL is sum of all samples
    assert isoforms_info['PB.103698.1'].FL == 25
    assert isoforms_info['PB.103704.2'].FL == 55
    assert isoforms_info['PB.103705.1'].FL == 5
    assert isoforms_info['PB.103707.1'].FL == 150
    assert isoforms_info['PB.103709.3'].FL == 20
    
    # Verify fields_class_cur extended with sample names (sorted)
    assert result_fields == ['isoform', 'length', 'FL.sample1', 'FL.sample2']
    
    # Verify logger info calls
    mock_qc_logger.info.assert_any_call("Multi-sample PacBio FL count format detected.")
    assert mock_qc_logger.warning.call_count == 0


def test_full_length_quantification_multi_demux_perfect_match(
    fl_count_multi_demux, create_isoforms_info, mock_qc_logger
):
    """Test multi-sample demux format FL count with all isoforms matching"""
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1', 'PB.103709.3']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_multi_demux, isoforms_info, fields
    )
    
    # Verify FL_dict populated correctly
    assert isoforms_info['PB.103698.1'].FL_dict == {'sample1': 10, 'sample2': 15}
    assert isoforms_info['PB.103704.2'].FL_dict == {'sample1': 25, 'sample2': 30}
    
    # Verify FL is sum of all samples
    assert isoforms_info['PB.103698.1'].FL == 25
    assert isoforms_info['PB.103704.2'].FL == 55
    
    # Verify fields extended
    assert 'FL.sample1' in result_fields
    assert 'FL.sample2' in result_fields
    
    # Verify logger
    mock_qc_logger.info.assert_any_call("Multi-sample PacBio FL count format detected.")


def test_full_length_quantification_isoform_in_fl_not_in_info(
    fl_count_single_sample, create_isoforms_info, mock_qc_logger
):
    """Test isoform present in FL file but missing from isoforms_info"""
    # Create isoforms_info missing one isoform that's in the FL file
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1']
    # Missing: PB.103709.3
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_single_sample, isoforms_info, fields
    )
    
    # Verify warning logged for missing isoform
    mock_qc_logger.warning.assert_any_call(
        "PB.103709.3 found in FL count file but not in input fasta."
    )
    
    # Verify other isoforms processed correctly
    assert isoforms_info['PB.103698.1'].FL == 10
    assert isoforms_info['PB.103704.2'].FL == 25
    assert isoforms_info['PB.103705.1'].FL == 5
    assert isoforms_info['PB.103707.1'].FL == 100


def test_full_length_quantification_isoform_in_info_not_in_fl_single(
    fl_count_single_sample, create_isoforms_info, mock_qc_logger
):
    """Test isoform present in isoforms_info but missing from FL file (single sample)"""
    # Create isoforms_info with extra isoform not in FL file
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1', 'PB.103709.3', 'PB.999999.1']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_single_sample, isoforms_info, fields
    )
    
    # Verify warning logged for missing isoform
    mock_qc_logger.warning.assert_any_call(
        "Isoform PB.999999.1 not found in FL count file. Assign count as 0."
    )
    
    # Verify missing isoform gets FL = 0
    assert isoforms_info['PB.999999.1'].FL == 0
    
    # Verify other isoforms processed correctly
    assert isoforms_info['PB.103698.1'].FL == 10
    assert isoforms_info['PB.103704.2'].FL == 25


def test_full_length_quantification_isoform_in_info_not_in_fl_multi(
    fl_count_multi_chain, create_isoforms_info, mock_qc_logger
):
    """Test isoform present in isoforms_info but missing from FL file (multi-sample)"""
    # Create isoforms_info with extra isoform not in FL file
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1', 'PB.103709.3', 'PB.888888.1']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_multi_chain, isoforms_info, fields
    )
    
    # Verify warning logged
    mock_qc_logger.warning.assert_any_call(
        "Isoform PB.888888.1 not found in FL count file. Assign count as 0."
    )
    
    # Verify missing isoform gets FL_dict with zeros and FL = 0
    assert isoforms_info['PB.888888.1'].FL_dict == {'sample1': 0, 'sample2': 0}
    assert isoforms_info['PB.888888.1'].FL == 0
    
    # Verify other isoforms processed correctly
    assert isoforms_info['PB.103698.1'].FL_dict == {'sample1': 10, 'sample2': 15}
    assert isoforms_info['PB.103698.1'].FL == 25


def test_full_length_quantification_bidirectional_mismatch(
    fl_count_single_sample, create_isoforms_info, mock_qc_logger
):
    """Test with isoforms missing in both directions"""
    # Create isoforms_info missing some from FL file, and having extras not in FL
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103707.1', 'PB.999999.1', 'PB.888888.1']
    # Missing from isoforms_info: PB.103705.1, PB.103709.3
    # Extra in isoforms_info: PB.999999.1, PB.888888.1
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_single_sample, isoforms_info, fields
    )
    
    # Verify warnings for both directions
    assert mock_qc_logger.warning.call_count >= 4
    
    # Warnings for isoforms in FL but not in info
    mock_qc_logger.warning.assert_any_call(
        "PB.103705.1 found in FL count file but not in input fasta."
    )
    mock_qc_logger.warning.assert_any_call(
        "PB.103709.3 found in FL count file but not in input fasta."
    )
    
    # Warnings for isoforms in info but not in FL
    mock_qc_logger.warning.assert_any_call(
        "Isoform PB.999999.1 not found in FL count file. Assign count as 0."
    )
    mock_qc_logger.warning.assert_any_call(
        "Isoform PB.888888.1 not found in FL count file. Assign count as 0."
    )
    
    # Verify matched isoforms processed correctly
    assert isoforms_info['PB.103698.1'].FL == 10
    assert isoforms_info['PB.103704.2'].FL == 25
    assert isoforms_info['PB.103707.1'].FL == 100
    
    # Verify missing isoforms get 0
    assert isoforms_info['PB.999999.1'].FL == 0
    assert isoforms_info['PB.888888.1'].FL == 0


def test_full_length_quantification_float_values(
    fl_count_float, create_isoforms_info, mock_qc_logger
):
    """Test that float FL count values are preserved"""
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_float, isoforms_info, fields
    )
    
    # Verify float values preserved
    assert isoforms_info['PB.103698.1'].FL == 10.5
    assert isoforms_info['PB.103704.2'].FL == 25.75
    assert isoforms_info['PB.103705.1'].FL == 5.25
    
    # Verify types
    assert isinstance(isoforms_info['PB.103698.1'].FL, float)
    assert isinstance(isoforms_info['PB.103704.2'].FL, float)


def test_full_length_quantification_zero_counts(
    fl_count_single_sample, create_isoforms_info, mock_qc_logger
):
    """Test that zero FL counts are handled correctly (not treated as missing)"""
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1', 'PB.103707.1', 'PB.103709.3']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_single_sample, isoforms_info, fields
    )
    
    # PB.103709.3 has FL count = 0 in the file
    assert isoforms_info['PB.103709.3'].FL == 0
    
    # Verify no warning for zero count (different from missing isoform)
    # Should only have info messages, no warnings
    assert mock_qc_logger.warning.call_count == 0


def test_full_length_quantification_with_na_values(
    fl_count_with_na, create_isoforms_info, mock_qc_logger
):
    """Test FL count file with NA values (converted to 0 by parser)"""
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform']
    
    # Call function
    result_info, result_fields = full_length_quantification(
        fl_count_with_na, isoforms_info, fields
    )
    
    # Verify NA values were converted to 0 by parser and handled correctly
    # PB.103698.1: sample1=10, sample2=0 (was NA), sample3=5
    assert isoforms_info['PB.103698.1'].FL_dict == {'sample1': 10, 'sample2': 0, 'sample3': 5}
    assert isoforms_info['PB.103698.1'].FL == 15
    
    # PB.103704.2: sample1=0 (was NA), sample2=30, sample3=15
    assert isoforms_info['PB.103704.2'].FL_dict == {'sample1': 0, 'sample2': 30, 'sample3': 15}
    assert isoforms_info['PB.103704.2'].FL == 45
    
    # PB.103705.1: sample1=5, sample2=0, sample3=0 (was NA)
    assert isoforms_info['PB.103705.1'].FL_dict == {'sample1': 5, 'sample2': 0, 'sample3': 0}
    assert isoforms_info['PB.103705.1'].FL == 5


def test_full_length_quantification_many_samples(
    fl_count_with_na, create_isoforms_info, mock_qc_logger
):
    """Test multi-sample with 3+ samples to verify field extension"""
    isoform_ids = ['PB.103698.1', 'PB.103704.2', 'PB.103705.1']
    isoforms_info = create_isoforms_info(isoform_ids)
    fields = ['isoform', 'length']
    
    # Call function (fl_count_with_na has 3 samples)
    _, result_fields = full_length_quantification(
        fl_count_with_na, isoforms_info, fields
    )
    
    # Verify all sample fields added in sorted order
    assert result_fields == ['isoform', 'length', 'FL.sample1', 'FL.sample2', 'FL.sample3']
    
    # Verify sorting (samples should be alphabetically sorted)
    fl_fields = [f for f in result_fields if f.startswith('FL.')]
    assert fl_fields == sorted(fl_fields)   


