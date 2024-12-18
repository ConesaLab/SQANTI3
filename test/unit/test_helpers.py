
from io import StringIO
from typing import Dict
import pytest
from unittest.mock import Mock, patch, mock_open
from Bio import SeqIO
from pathlib import Path
import sys, os


main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

from src.helpers import (
    process_gtf_line, rename_isoform_seqids, get_corr_filenames, get_isoform_hits_name, 
    get_class_junc_filenames, get_omitted_name
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

## process_gtf_line ##

@pytest.fixture
def genome_dict() -> Dict[str, str]:
    return {"chr1": "ATCG", "chr2": "GCTA"}

@pytest.fixture
def mock_corrGTF_out():
    return main_path + "/test/test_data/corrected.gtf"

@pytest.fixture
def mock_discard_gtf():
    return main_path + "/test/test_data/discard.gtf"

def test_comment_line(genome_dict, mock_corrGTF_out, mock_discard_gtf):
    line = "# This is a comment\n"
    result = process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    assert result == None
    assert not os.path.exists(mock_corrGTF_out)
    assert not os.path.exists(mock_discard_gtf)


def test_malformed_line(genome_dict, mock_corrGTF_out, mock_discard_gtf, capsys):
    line = "chr1\tgene\n"
    process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    captured = capsys.readouterr()
    assert "WARNING: Skipping malformed GTF line" in captured.out
    assert not os.path.exists(mock_corrGTF_out)
    assert not os.path.exists(mock_discard_gtf)

def test_chromosome_not_in_genome(genome_dict, mock_corrGTF_out, mock_discard_gtf):
    line = "chr3\tEnsembl\texon\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
    with pytest.raises(ValueError, match="ERROR: GTF chromosome 'chr3' not found in genome reference file."):
        process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)

def test_valid_transcript_line(genome_dict, mock_corrGTF_out, mock_discard_gtf):
    line = "chr1\tEnsembl\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
    
    process_gtf_line(line, genome_dict, Path(mock_corrGTF_out), mock_discard_gtf)
    with open(mock_corrGTF_out, "r") as f:
        assert f.read() == line
    assert not os.path.exists(mock_discard_gtf)
    os.remove(mock_corrGTF_out)

def test_valid_exon_line(genome_dict, mock_corrGTF_out, mock_discard_gtf):
    line = "chr2\tEnsembl\texon\t1\t1000\t.\t-\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
    process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    with open(mock_corrGTF_out, "r") as f:
        assert f.read() == line
    assert not os.path.exists(mock_discard_gtf)
    os.remove(mock_corrGTF_out)

def test_unknown_strand(genome_dict, mock_corrGTF_out, mock_discard_gtf, capsys):
    line = "chr1\tEnsembl\texon\t1\t1000\t.\t.\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
    process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    captured = capsys.readouterr()
    assert "WARNING: Discarding unknown strand transcript" in captured.out
    assert not os.path.exists(mock_corrGTF_out)
    with open(mock_discard_gtf, "r") as f:
        assert f.read() == line
    os.remove(mock_discard_gtf)

def test_non_transcript_exon_line(genome_dict, mock_corrGTF_out, mock_discard_gtf):
    line = "chr1\tEnsembl\tgene\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\";\n"
    process_gtf_line(line, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    assert not os.path.exists(mock_corrGTF_out)
    assert not os.path.exists(mock_discard_gtf)

def test_adding_lines_to_file(genome_dict,mock_corrGTF_out,mock_discard_gtf):
    line1 = "chr1\tEnsembl\ttranscript\t1\t1000\t.\t+\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
    line2 = "chr2\tEnsembl\texon\t1\t1000\t.\t-\t.\tgene_id \"ENSG00000223972\"; transcript_id \"ENST00000456328\";\n"
    process_gtf_line(line1, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    process_gtf_line(line2, genome_dict, mock_corrGTF_out, mock_discard_gtf)
    with open(mock_corrGTF_out, "r") as f:
        assert f.read() == line1 + line2
    assert not os.path.exists(mock_discard_gtf)
    os.remove(mock_corrGTF_out)
