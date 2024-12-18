import pytest
from unittest.mock import Mock, patch
from Bio import SeqIO
from pathlib import Path
import sys, os


main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.helpers import (
    rename_isoform_seqids, get_corr_filenames, get_isoform_hits_name, 
    get_class_junc_filenames, get_omitted_name, sequence_correction
)

### rename_isoform_seqids ### 

@pytest.fixture
def mock_fasta_file():
    """Fixture to return the path to the mock FASTA file."""
    return Path(main_path,"test","test_data","isoform_mock.fasta")

@pytest.fixture
def mock_fastq_file():
    """Fixture to return the path to the mock FASTA file."""
    return Path(main_path,"test","test_data","isoform_mock.fastq")

@pytest.fixture
def mock_gz_file():
    """Fixture to return the path to the mock FASTA file."""
    return Path(main_path,"test","test_data","isoform_mock.fasta.gz")

def test_rename_isoform_seqids_fasta(mock_fasta_file):
    # Arrange
    output_fasta = str(mock_fasta_file.with_name("isoform_mock.renamed.fasta"))
    expected_ids = ["PB.1.1", "PB.2.1", "PBfusion.3.1"]

    # Act
    result_file = rename_isoform_seqids(str(mock_fasta_file))

    # Assert
    assert result_file == output_fasta, "Output file name does not match expected."
    assert os.path.exists(output_fasta), "Output file was not created."

    # Verify the content of the output FASTA file
    with open(output_fasta, "r") as output:
        renamed_records = list(SeqIO.parse(output, "fasta"))
        assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

        for record, expected_id in zip(renamed_records, expected_ids):
            assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

    # Cleanup: Remove the output file after test
    os.remove(output_fasta)

# Reads a fastq file and returns the fasta file with IDs 
def test_rename_isoform_seqids_fastq(mock_fastq_file):
    # Arrange
    output_fastq = str(mock_fastq_file.with_name("isoform_mock.renamed.fasta"))
    expected_ids = ["PB.1.1", "PB.2.1", "PBfusion.3.1"]

    # Act
    result_file = rename_isoform_seqids(str(mock_fastq_file))

    # Assert
    assert result_file == output_fastq, "Output file name does not match expected."
    assert os.path.exists(output_fastq), "Output file was not created."

    # Verify the content of the output FASTA file
    with open(output_fastq, "r") as output:
        renamed_records = list(SeqIO.parse(output, "fasta"))
        assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

        for record, expected_id in zip(renamed_records, expected_ids):
            assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

    # Cleanup: Remove the output file after test
    os.remove(output_fastq)

# Tests if it can read a gz file
def test_rename_isoform_seqids_gz(mock_gz_file):
    # Arrange
    output_fastq = str(mock_gz_file.with_name("isoform_mock.renamed.fasta"))
    expected_ids = ["PB.1.1", "PB.2.1", "PBfusion.3.1"]

    # Act
    result_file = rename_isoform_seqids(str(mock_gz_file))

    # Assert
    assert result_file == output_fastq, "Output file name does not match expected."
    assert os.path.exists(output_fastq), "Output file was not created."

    # Verify the content of the output FASTA file
    with open(output_fastq, "r") as output:
        renamed_records = list(SeqIO.parse(output, "fasta"))
        assert len(renamed_records) == len(expected_ids), "Number of records mismatch."

        for record, expected_id in zip(renamed_records, expected_ids):
            assert record.id == expected_id, f"Expected ID {expected_id} but got {record.id}."

    # Cleanup: Remove the output file after test
    os.remove(output_fastq)

def test_rename_isoform_seqids_invalid_ids(mock_fasta_file):
    # Arrange
    invalid_fasta = str(mock_fasta_file.with_name("invalid_mock.fasta"))
    with open(invalid_fasta, "w") as invalid:
        invalid.write(">Invalid_ID\nACGTACGTACGT\n")

    # Act & Assert
    with pytest.raises(SystemExit) as excinfo:
        rename_isoform_seqids(invalid_fasta)

    assert str(excinfo.value) == "1", "Expected a SystemExit with code 1 for invalid IDs."

    # Cleanup: Remove the invalid mock file after test
    os.remove(invalid_fasta)

def test_rename_isoform_seqids_force_ignore(mock_fasta_file):
    # Arrange
    output_fasta = str(mock_fasta_file.with_name("isoform_mock.renamed.fasta"))

    # Act
    result_file = rename_isoform_seqids(str(mock_fasta_file), force_id_ignore=True)

    # Assert
    assert result_file == output_fasta, "Output file name does not match expected."
    assert os.path.exists(output_fasta), "Output file was not created."

    # Verify all records are renamed
    with open(output_fasta, "r") as output:
        renamed_records = list(SeqIO.parse(output, "fasta"))
        assert all(record.id for record in renamed_records), "Some records were not renamed."

    # Cleanup: Remove the output file after test
    os.remove(output_fasta)

### Input reading functions ###

@pytest.fixture
def sample_paths():
    return {
        'outdir': '/tmp',
        'prefix': 'test'
    }

def test_get_corr_filenames(sample_paths):
    expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
    corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF = get_corr_filenames(**sample_paths)
    
    assert corrGTF == expected_prefix + "_corrected.gtf"
    assert corrSAM == expected_prefix + "_corrected.sam"
    assert corrFASTA == expected_prefix + "_corrected.fasta"
    assert corrORF == expected_prefix + "_corrected.faa"
    assert corrCDS_GTF_GFF == expected_prefix + "_corrected.gtf.cds.gff"

def test_get_isoform_hits_name(sample_paths):
    expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
    isoform_hits_name = get_isoform_hits_name(**sample_paths)
    
    assert isoform_hits_name == expected_prefix + "_isoform_hits.txt"

def test_get_class_junc_filenames(sample_paths):
    expected_prefix = os.path.join(sample_paths['outdir'], sample_paths['prefix'])
    outputClassPath, outputJuncPath = get_class_junc_filenames(**sample_paths)
    
    assert outputClassPath == expected_prefix + "_classification.txt"
    assert outputJuncPath == expected_prefix + "_junctions.txt"

def test_get_omitted_name(sample_paths):
    expected_prefix = os.path.join(os.path.join(sample_paths['outdir'], sample_paths['prefix']))
    omitted_name = get_omitted_name(**sample_paths)
    
    assert omitted_name == expected_prefix + "_omitted_due_to_min_ref_len.txt"

@pytest.mark.parametrize("outdir,prefix", [
    ("/var/log", "app"),
    ("/home/user/documents", "report"),
    (".", "test"),
    ("/tmp/my folder", "test file"),
])
def test_all_functions_parametrized(outdir, prefix):
    expected_prefix = os.path.abspath(os.path.join(outdir, prefix))
    
    # Test get_corr_filenames
    corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF = get_corr_filenames(outdir, prefix)
    assert corrGTF == expected_prefix + "_corrected.gtf"
    assert corrSAM == expected_prefix + "_corrected.sam"
    assert corrFASTA == expected_prefix + "_corrected.fasta"
    assert corrORF == expected_prefix + "_corrected.faa"
    assert corrCDS_GTF_GFF == expected_prefix + "_corrected.gtf.cds.gff"
    
    # Test get_isoform_hits_name
    isoform_hits_name = get_isoform_hits_name(outdir, prefix)
    assert isoform_hits_name == expected_prefix + "_isoform_hits.txt"
    
    # Test get_class_junc_filenames
    outputClassPath, outputJuncPath = get_class_junc_filenames(outdir, prefix)
    assert outputClassPath == expected_prefix + "_classification.txt"
    assert outputJuncPath == expected_prefix + "_junctions.txt"
    
    # Test get_omitted_name
    omitted_name = get_omitted_name(outdir, prefix)
    assert omitted_name == expected_prefix + "_omitted_due_to_min_ref_len.txt"

def test_empty_prefix():
    outdir = "/tmp"
    prefix = ""
    expected_prefix = os.path.abspath(os.path.join(outdir, prefix))
    
    # Test all functions with empty prefix
    corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF = get_corr_filenames(outdir, prefix)
    isoform_hits_name = get_isoform_hits_name(outdir, prefix)
    outputClassPath, outputJuncPath = get_class_junc_filenames(outdir, prefix)
    omitted_name = get_omitted_name(outdir, prefix)
    
    assert all(path.startswith(expected_prefix) for path in [corrGTF, corrSAM, corrFASTA, corrORF, corrCDS_GTF_GFF, isoform_hits_name, outputClassPath, outputJuncPath, omitted_name])

### sequence_correction ###


@pytest.fixture
def mock_file_ops():
    return Mock()

@pytest.fixture
def mock_cmd_runner():
    return Mock()

@pytest.fixture
def mock_cmd_templates():
    return {
        "gmap": "gmap_cmd {cpus} {dir} {name} {sense} {i} {o}",
        "minimap2": "minimap2_cmd {cpus} {sense} {g} {i} {o}",
        "deSALT": "desalt_cmd {cpus} {dir} {i} {o}",
        "uLTRA": "ultra_cmd {cpus} {prefix} {g} {a} {i} {o_dir}",
        "GFFREAD": "gffread_cmd"
    }

@pytest.fixture
def default_args():
    return {
        "outdir": "/out",
        "output": "test",
        "cpus": 4,
        "chunks": 2,
        "fasta": True,
        "genome_dict": {"chr1": "ATCG"},
        "badstrandGTF": "bad.gtf",
        "genome": "genome.fa",
        "isoforms": "iso.fa",
        "aligner_choice": "minimap2",
        "gmap_index": None,
        "sense": False,
        "annotation": None
    }

def test_sequence_correction_fasta_exists(mock_file_ops, mock_cmd_runner, mock_cmd_templates, default_args):
    mock_file_ops.exists.return_value = True
    
    with patch('src.helpers.get_corr_filenames', return_value=('corrGTF', 'corrSAM', 'corrFASTA', None, None)):
        sequence_correction(**default_args)
    
    mock_file_ops.exists.assert_called_once()
    mock_cmd_runner.assert_not_called()

# def test_sequence_correction_fasta_not_exists(mock_file_ops, mock_cmd_runner, mock_cmd_templates, default_args):
#     mock_file_ops.exists.return_value = False
#     mock_file_ops.splitext.return_value = ('corrSAM', '')
#     mock_file_ops.basename.return_value = 'corrSAM'
    
#     with patch('src.helpers.get_corr_filenames', return_value=('corrGTF', 'corrSAM', 'corrFASTA', None, None)):
#         with patch('src.helpers.err_correct'):
#             with patch('src.helpers.convert_sam_to_gff3'):
#                 sequence_correction(**default_args, file_ops=mock_file_ops, cmd_runner=mock_cmd_runner, cmd_templates=mock_cmd_templates)
    
#     assert mock_cmd_runner.call_count == 2  # One for alignment, one for GTF conversion

# @pytest.mark.parametrize("aligner_choice", ["gmap", "minimap2", "deSALT", "uLTRA"])
# def test_sequence_correction_different_aligners(mock_file_ops, mock_cmd_runner, mock_cmd_templates, default_args, aligner_choice):
#     mock_file_ops.exists.return_value = False
#     mock_file_ops.splitext.return_value = ('corrSAM', '')
#     mock_file_ops.basename.return_value = 'corrSAM'
#     default_args['aligner_choice'] = aligner_choice
    
#     with patch('src.helpers.get_corr_filenames', return_value=('corrGTF', 'corrSAM', 'corrFASTA', None, None)):
#         with patch('src.helpers.err_correct'):
#             with patch('src.helpers.convert_sam_to_gff3'):
#                 sequence_correction(**default_args, file_ops=mock_file_ops, cmd_runner=mock_cmd_runner, cmd_templates=mock_cmd_templates)
    
#     assert mock_cmd_runner.call_count == 2
#     assert aligner_choice in mock_cmd_runner.call_args_list[0][0][0]

# def test_sequence_correction_gtf_input(mock_file_ops, mock_cmd_runner, mock_cmd_templates, default_args):
#     mock_file_ops.exists.return_value = False
#     default_args['fasta'] = False
    
#     with patch('src.helpers.get_corr_filenames', return_value=('corrGTF', 'corrSAM', 'corrFASTA', None, None)):
#         sequence_correction(**default_args, file_ops=mock_file_ops, cmd_runner=mock_cmd_runner, cmd_templates=mock_cmd_templates)
    
#     assert mock_cmd_runner.call_count == 1  # Only for GTF to FASTA conversion
#     assert "gffread_cmd" in mock_cmd_runner.call_args[0][0]

# def test_sequence_correction_invalid_aligner(mock_file_ops, mock_cmd_runner, mock_cmd_templates, default_args):
#     default_args['aligner_choice'] = 'invalid_aligner'
    
#     with pytest.raises(ValueError, match="Unsupported aligner choice: invalid_aligner"):
#         sequence_correction(**default_args, file_ops=mock_file_ops, cmd_runner=mock_cmd_runner, cmd_templates=mock_cmd_templates)

# def test_sequence_correction_gtf_chromosome_not_in_genome(mock_file_ops, mock_cmd_runner, mock_cmd_templates, default_args):
#     mock_file_ops.exists.return_value = False
#     default_args['fasta'] = False
#     mock_file_ops.open.return_value.__enter__.return_value.readlines.return_value = [
#         "chr2\tsome\texon\t1\t100\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"trans1\";\n"
#     ]
    
#     with pytest.raises(ValueError, match="ERROR: gtf chromosome 'chr2' not found in genome reference file."):
#         sequence_correction(**default_args, file_ops=mock_file_ops, cmd_runner=mock_cmd_runner, cmd_templates=mock_cmd_templates)