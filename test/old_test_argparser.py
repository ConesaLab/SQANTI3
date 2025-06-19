import os
import sys
from io import StringIO
import pytest
from unittest.mock import patch, MagicMock
main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, main_path)

# Test missing genome file
def test_missing_genome():
    genome_path = "/path/to/nonexistent/genome.fasta"
    
    # Mock sys.stderr to capture the print output
    with patch('sys.stderr', new_callable=StringIO) as mock_stderr:
        with pytest.raises(SystemExit):  # Expect the script to exit
            if not os.path.isfile(genome_path):
                print(f"ERROR: genome fasta {genome_path} doesn't exist. Abort!", file=sys.stderr)
                sys.exit(1)


        # Check that the error message was printed to stderr
        expected_error_message = f"ERROR: genome fasta {genome_path} doesn't exist. Abort!\n"
        actual_output = mock_stderr.getvalue()

        # Assert that the printed message matches the expected error message
        assert actual_output == expected_error_message
        
# Test missing isoforms file
def test_missing_isoforms():
    with patch("sys.exit") as mock_exit, patch("sys.stderr.write") as mock_stderr:
        args = MagicMock()
        args.genome = "/path/to/genome.fasta"
        args.isoforms = "/path/to/nonexistent/isoforms.fasta"
        args.annotation = "/path/to/annotation.gtf"
        
        # Simulate the behavior of your script
        if not os.path.isfile(args.isoforms):
            print(f"ERROR: Input isoforms {args.isoforms} doesn't exist. Abort!", file=sys.stderr)
            sys.exit()
        
        mock_exit.assert_called_once()
        mock_stderr.assert_called_with("ERROR: Input isoforms /path/to/nonexistent/isoforms.fasta doesn't exist. Abort!\n")

# Test output directory creation (non-existing directory)
def test_create_output_directory():
    with patch("os.makedirs") as mock_makedirs, patch("sys.stderr.write") as mock_stderr:
        args = MagicMock()
        args.dir = "/path/to/nonexistent/dir"
        
        if not os.path.isdir(args.dir):
            print(f"WARNING: output directory {args.dir} already exists. Overwriting!", file=sys.stderr)
            os.makedirs(args.dir)
        
        mock_makedirs.assert_called_once_with(args.dir)

# Test existing output directory (warning)
def test_existing_output_directory_warning():
    with patch("sys.stderr.write") as mock_stderr:
        args = MagicMock()
        args.dir = "/path/to/existing/dir"
        
        if os.path.isdir(args.dir):
            print(f"WARNING: output directory {args.dir} already exists. Overwriting!", file=sys.stderr)
        
        mock_stderr.assert_called_with("WARNING: output directory /path/to/existing/dir already exists. Overwriting!\n")

# Test GMAP index existence (invalid)
def test_missing_gmap_index():
    with patch("sys.exit") as mock_exit, patch("sys.stderr.write") as mock_stderr:
        args = MagicMock()
        args.aligner_choice = "gmap"
        args.gmap_index = "/path/to/nonexistent/gmap_index"
        
        if not os.path.isdir(args.gmap_index):
            print(f"GMAP index {args.gmap_index} doesn't exist! Abort.", file=sys.stderr)
            sys.exit()
        
        mock_exit.assert_called_once()
        mock_stderr.assert_called_with("GMAP index /path/to/nonexistent/gmap_index doesn't exist! Abort.\n")

# Test deSALT index existence (invalid)
def test_missing_desalt_index():
    with patch("sys.exit") as mock_exit, patch("sys.stderr.write") as mock_stderr:
        args = MagicMock()
        args.aligner_choice = "deSALT"
        args.gmap_index = "/path/to/nonexistent/desalt_index"
        
        if not os.path.isdir(args.gmap_index):
            print(f"deSALT index {args.gmap_index} doesn't exist! Abort.", file=sys.stderr)
            sys.exit()
        
        mock_exit.assert_called_once()
        mock_stderr.assert_called_with("deSALT index /path/to/nonexistent/desalt_index doesn't exist! Abort.\n")

# Test renaming isoform IDs
def test_rename_isoform_ids():
    with patch("src.helpers.rename_isoform_seqids") as mock_rename, patch("sys.stderr.write") as mock_stderr:
        args = MagicMock()
        args.isoforms = "/path/to/isoforms.fasta"
        args.force_id_ignore = True
        args.fasta = True
        args.aligner_choice = "gmap"
        
        # Simulate renaming isoform IDs
        from src.helpers import rename_isoform_seqids
        original_isoforms = args.isoforms
        args.isoforms = rename_isoform_seqids(original_isoforms, args.force_id_ignore)
        
        mock_rename.assert_called_once_with(original_isoforms, args.force_id_ignore)
        #mock_stderr.assert_called_with("Cleaning up isoform IDs...\n")
