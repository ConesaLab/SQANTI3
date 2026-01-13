"""
Shared fixtures for process-level tests.
This module provides fixtures for building isoforms_info dictionaries
and other complex data structures needed for testing QC pipeline functions.
"""
import os
import sys
import pytest
from collections import defaultdict
from dataclasses import dataclass

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.qc_classes import myQueryTranscripts


@pytest.fixture
def sample_isoforms_info():
    """
    Create a sample isoforms_info dictionary with realistic test data.
    This represents the state after classification but before QC computations.
    
    Returns:
        dict: Dictionary of pbid -> myQueryTranscripts objects
    """
    isoforms_info = {}
    
    # Create test isoforms matching those in test_data
    isoforms_info["PB.124830.1"] = myQueryTranscripts(
        isoform="PB.124830.1",
        chrom="chr22",
        start=39349957,
        end=39377911,
        strand="+",
        length=2000,
        exons=2,
        structural_category="novel_in_catalog",
        genes=["GENE1"],
        transcripts=["ENST00001.1"],
        ref_length=2000,
        ref_exons=2,
        diff_to_TSS=0,
        diff_to_TTS=0
    )
    
    isoforms_info["PB.103714.1"] = myQueryTranscripts(
        isoform="PB.103714.1",
        chrom="chr22",
        start=21919422,
        end=21920084,
        strand="-",
        length=663,
        exons=1,
        structural_category="full-splice_match",
        genes=["GENE2"],
        transcripts=["ENST00002.1"],
        ref_length=663,
        ref_exons=1,
        diff_to_TSS=0,
        diff_to_TTS=0
    )
    
    isoforms_info["PB.103724.1"] = myQueryTranscripts(
        isoform="PB.103724.1",
        chrom="chr22",
        start=25448098,
        end=25459681,
        strand="+",
        length=1200,
        exons=4,
        structural_category="novel_not_in_catalog",
        genes=["GENE3"],
        transcripts=["ENST00003.1"],
        ref_length=1100,
        ref_exons=3,
        diff_to_TSS=100,
        diff_to_TTS=-50
    )
    
    isoforms_info["PB.103781.1"] = myQueryTranscripts(
        isoform="PB.103781.1",
        chrom="chr22",
        start=38741566,
        end=38742142,
        strand="+",
        length=577,
        exons=1,
        structural_category="full-splice_match",
        genes=["GENE4"],
        transcripts=["ENST00004.1"],
        ref_length=577,
        ref_exons=1,
        diff_to_TSS=0,
        diff_to_TTS=0
    )
    
    isoforms_info["PB.23068.2"] = myQueryTranscripts(
        isoform="PB.23068.2",
        chrom="chr22",
        start=37805229,
        end=37807429,
        strand="+",
        length=1500,
        exons=2,
        structural_category="full-splice_match",
        genes=["GENE5"],
        transcripts=["ENST00005.1"],
        ref_length=1500,
        ref_exons=2,
        diff_to_TSS=0,
        diff_to_TTS=0
    )
    
    return isoforms_info


@pytest.fixture
def fields_class_cur():
    """
    Return the base classification fields list.
    This list gets extended by various QC functions.
    """
    return [
        "isoform", "chrom", "strand", "length", "exons",
        "structural_category", "associated_gene", "associated_transcript",
        "ref_length", "ref_exons", "diff_to_TSS", "diff_to_TTS",
        "subcategory", "RTS_stage", "all_canonical", "min_sample_cov",
        "min_cov", "sd_cov", "FL", "n_indels", "n_indels_junc",
        "bite", "iso_exp", "gene_exp", "ratio_exp", "FSM_class",
        "coding", "ORF_length", "CDS_length", "CDS_start", "CDS_end",
        "CDS_genomic_start", "CDS_genomic_end"
    ]


@pytest.fixture
def fl_count_file_single():
    """Path to single-sample FL count test file."""
    return os.path.join(main_path, "test/test_data/abundance/fl_count_single_sample.tsv")


@pytest.fixture
def fl_count_file_multi():
    """Path to multi-sample FL count test file."""
    return os.path.join(main_path, "test/test_data/abundance/fl_count_multi_sample.tsv")


@pytest.fixture
def junction_file():
    """Path to test junctions file."""
    return os.path.join(main_path, "test/test_data/junctions/test_junctions.txt")


@pytest.fixture
def genome_file():
    """Path to test genome file."""
    return os.path.join(main_path, "test/test_data/genome/genome_test.fasta")
