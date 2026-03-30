"""
Shared fixtures for requantification tests.
This module provides reusable fixtures for building test data structures
needed for testing the requantification pipeline.
"""
import os
import sys
import pytest
import pandas as pd

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)


@pytest.fixture
def sample_counts_simple():
    """
    Simple count matrix with 3 isoforms and 2 samples.
    
    Returns:
        pd.DataFrame: Count matrix
    """
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
        'sample1': [100, 200, 50],
        'sample2': [150, 250, 60]
    })


@pytest.fixture
def sample_counts_complex():
    """
    Complex count matrix with multiple artifacts and samples.
    
    Returns:
        pd.DataFrame: Count matrix
    """
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1', 'artifact1', 'artifact2', 'artifact3'],
        'sample1': [1000, 2000, 1500, 100, 200, 150],
        'sample2': [1200, 2200, 1600, 120, 180, 140],
        'sample3': [900, 1800, 1400, 90, 210, 160]
    })


@pytest.fixture
def sample_rescue_simple():
    """
    Simple rescue mapping with one artifact mapped to one target.
    
    Returns:
        pd.DataFrame: Rescue mapping
    """
    return pd.DataFrame({
        'artifact': ['artifact1'],
        'assigned_transcript': ['PB.1.1']
    })


@pytest.fixture
def sample_rescue_multi_target():
    """
    Rescue mapping where one artifact maps to multiple targets.
    
    Returns:
        pd.DataFrame: Rescue mapping
    """
    return pd.DataFrame({
        'artifact': ['artifact1', 'artifact1'],
        'assigned_transcript': ['PB.1.1', 'PB.2.1']
    })


@pytest.fixture
def sample_rescue_complex():
    """
    Complex rescue mapping with various scenarios.
    
    Returns:
        pd.DataFrame: Rescue mapping
    """
    return pd.DataFrame({
        'artifact': ['artifact1', 'artifact2', 'artifact3', 'artifact3'],
        'assigned_transcript': ['PB.1.1', 'PB.2.1', 'PB.1.1', 'PB.3.1']
    })


@pytest.fixture
def sample_classification_simple():
    """
    Simple classification table.
    
    Returns:
        pd.DataFrame: Classification table
    """
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
        'filter_result': ['Isoform', 'Isoform', 'Artifact'],
        'associated_gene': ['GENE1', 'GENE2', 'GENE1'],
        'length': [1000, 2000, 1500]
    })


@pytest.fixture
def sample_classification_complex():
    """
    Complex classification table with multiple scenarios.
    
    Returns:
        pd.DataFrame: Classification table
    """
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'PB.3.1', 'artifact1', 'artifact2', 'artifact3'],
        'filter_result': ['Isoform', 'Isoform', 'Isoform', 'Artifact', 'Artifact', 'Artifact'],
        'associated_gene': ['GENE1', 'GENE2', 'GENE3', 'GENE1', 'GENE2', 'novel_gene_1'],
        'length': [1000, 2000, 1500, 1200, 1800, 1600]
    })


@pytest.fixture
def sample_counts_file(tmp_path, sample_counts_simple):
    """
    Create a temporary count file.
    
    Args:
        tmp_path: pytest fixture for temporary directory
        sample_counts_simple: Simple count matrix fixture
    
    Returns:
        Path: Path to the temporary count file
    """
    counts_file = tmp_path / "test_counts.tsv"
    sample_counts_simple.to_csv(counts_file, sep='\t', index=False)
    return counts_file


@pytest.fixture
def expected_reassigned_counts_simple():
    """
    Expected reassigned counts for simple test case.
    After artifact1 is redistributed to PB.1.1.
    
    Returns:
        pd.DataFrame: Expected reassigned counts
    """
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1'],
        'sample1': [150, 200],  # PB.1.1 gets artifact1's 50
        'sample2': [210, 250]   # PB.1.1 gets artifact1's 60
    })


@pytest.fixture
def sample_zero_count_artifact():
    """
    Count matrix with artifact having zero counts.
    
    Returns:
        pd.DataFrame: Count matrix
    """
    return pd.DataFrame({
        'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
        'sample1': [100, 200, 0],
        'sample2': [150, 250, 0]
    })


@pytest.fixture
def sample_all_artifacts():
    """
    Count matrix with only artifacts (no valid isoforms).
    
    Returns:
        pd.DataFrame: Count matrix
    """
    return pd.DataFrame({
        'isoform': ['artifact1', 'artifact2', 'artifact3'],
        'sample1': [100, 200, 150],
        'sample2': [120, 180, 140]
    })


@pytest.fixture
def sample_classification_all_artifacts():
    """
    Classification table where all transcripts are artifacts.
    
    Returns:
        pd.DataFrame: Classification table
    """
    return pd.DataFrame({
        'isoform': ['artifact1', 'artifact2', 'artifact3'],
        'filter_result': ['Artifact', 'Artifact', 'Artifact'],
        'associated_gene': ['GENE1', 'GENE2', 'novel_gene_1'],
        'length': [1000, 2000, 1500]
    })
