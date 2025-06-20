import sys,os,pytest
import json
import pandas as pd
from Bio import SeqIO

main_path=os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, main_path)

# filter_isoforms
@pytest.fixture
def fasta_file():
    return os.path.join(main_path, "test", "test_data", "test_isoforms.fasta")

@pytest.fixture
def fastq_file():
    return os.path.join(main_path, "test", "test_data", "test_isoforms.fastq")

@pytest.fixture
def isoform_list():
    return ["PB.3853.1","PB.124830.1","PB.103763.1","PB.103724.1"]

@pytest.fixture
def prefix():
    return os.path.join(main_path, "test", "test_data")

from src.filter_output import filter_isoforms

def test_filter_isoforms(fasta_file, prefix, isoform_list):
    filter_isoforms(fasta_file, prefix, isoform_list)
    filtered_file = f"{prefix}.filtered.fasta"
    with open(filtered_file, 'r') as f:
        filtered_ids = [record.id for record in SeqIO.parse(f, 'fasta')]
    assert len(filtered_ids) == len(isoform_list)
    assert all(id in filtered_ids for id in isoform_list)
    os.remove(filtered_file)

def test_filter_isoforms_empty(fasta_file, prefix):
    isoform_list = []
    filter_isoforms(fasta_file, prefix, isoform_list)
    filtered_file = f"{prefix}.filtered.fasta"
    with open(filtered_file, 'r') as f:
        filtered_ids = [record.id for record in SeqIO.parse(f, 'fasta')]
    assert len(filtered_ids) == 0
    os.remove(filtered_file)

def test_filter_isoforms_fastq(fastq_file, prefix, isoform_list):
    filter_isoforms(fastq_file, prefix, isoform_list)
    filtered_file = f"{prefix}.filtered.fastq"
    with open(filtered_file, 'r') as f:
        filtered_ids = [record.id for record in SeqIO.parse(f, 'fastq')]
    assert len(filtered_ids) == len(isoform_list)
    assert all(id in filtered_ids for id in isoform_list)
    os.remove(filtered_file)

# filter_gtf
from src.filter_output import filter_gtf

@pytest.fixture
def gtf_file():
    return os.path.join(main_path, "test", "test_data", "test_isoforms.gtf")

def test_filter_gtf(gtf_file, prefix, isoform_list):
    filter_gtf(gtf_file, prefix, isoform_list)
    filtered_file = f"{prefix}.filtered.gtf"
    with open(filtered_file, 'r') as f:
        filtered_lines = f.readlines()
    assert len(filtered_lines) > 0  # Ensure some lines are present
    assert all(any(id in line for id in isoform_list) for line in filtered_lines)
    os.remove(filtered_file)

def test_filter_gtf_empty(gtf_file, prefix):
    isoform_list = []
    filter_gtf(gtf_file, prefix, isoform_list)
    filtered_file = f"{prefix}.filtered.gtf"
    with open(filtered_file, 'r') as f:
        filtered_lines = f.readlines()
    assert len(filtered_lines) == 0
    os.remove(filtered_file)

# filter_faa
from src.filter_output import filter_faa
@pytest.fixture
def faa_file():
    return os.path.join(main_path, "test", "test_data", "TD2_test.faa")

@pytest.fixture
def faa_isoform_list():
    return ["PB.83093.1.p1"]

def test_filter_faa(faa_file, prefix, faa_isoform_list):
    filter_faa(faa_file, prefix, faa_isoform_list)
    filtered_file = f"{prefix}.filtered.faa"
    with open(filtered_file, 'r') as f:
        filtered_lines = f.readlines()
    assert len(filtered_lines) > 0  # Ensure some lines are present
    assert all(any(id in line for id in faa_isoform_list) for line in filtered_lines if line.startswith('>'))
    os.remove(filtered_file)
