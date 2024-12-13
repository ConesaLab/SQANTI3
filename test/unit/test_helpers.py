import pytest
from io import StringIO
from unittest.mock import mock_open, patch
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pathlib import Path
import sys, os

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.helpers import (
    rename_isoform_seqids
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


### write_collapsed_GFF_with_CDS ###

