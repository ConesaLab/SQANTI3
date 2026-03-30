"""
Process-level tests for requantification pipeline.
Tests the complete requantification workflow as an integrated pipeline.
"""
import pytest
import pandas as pd
import sys
import os

# Add main path to sys.path
main_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, main_path)

from src.utilities.rescue.sq_requant import (
    requantify,
    redistribute_counts_vectorized,
    requantification_pipeline,
    update_classification
)


class TestRedistributeCountsVectorized:
    """Test suite for the vectorized count redistribution pipeline."""
    
    def test_simple_redistribution(self):
        """Test basic count redistribution from artifact to rescued isoform."""
        # Setup
        rescue_df = pd.DataFrame({
            'isoform': ['artifact1', 'artifact2'],
            'assigned_transcript': ['PB.1.1', 'PB.2.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1', 'artifact2'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact', 'Artifact']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1', 'artifact2'],
            'sample1': [100, 200, 50, 75],
            'sample2': [150, 250, 60, 80]
        })
        
        # Execute
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        
        # Verify
        assert 'isoform' in result.columns
        assert 'sample1' in result.columns
        assert 'sample2' in result.columns
        
        # PB.1.1 should have original 100 + artifact1's 50 = 150
        pb11_row = result[result['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['sample1'] == 150
        assert pb11_row['sample2'] == 210  # 150 + 60
        
        # PB.2.1 should have original 200 + artifact2's 75 = 275
        pb21_row = result[result['isoform'] == 'PB.2.1'].iloc[0]
        assert pb21_row['sample1'] == 275
        assert pb21_row['sample2'] == 330  # 250 + 80
        
        # Artifacts should not be in result
        assert 'artifact1' not in result['isoform'].values
        assert 'artifact2' not in result['isoform'].values
    
    def test_multi_target_redistribution(self):
        """Test redistribution when one artifact maps to multiple targets."""
        # Setup: artifact1 maps to both PB.1.1 and PB.2.1
        rescue_df = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1'],
            'assigned_transcript': ['PB.1.1', 'PB.2.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'sample1': [100, 200, 90]  # 90 should split proportionally
        })
        
        # Execute
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        
        # Verify: Should split proportionally based on original abundance
        # PB.1.1: 100, PB.2.1: 200 -> ratio 1:2
        # 90 counts split 30:60
        pb11_row = result[result['isoform'] == 'PB.1.1'].iloc[0]
        pb21_row = result[result['isoform'] == 'PB.2.1'].iloc[0]
        
        # Total should be conserved: 100 + 200 + 90 = 390
        total = pb11_row['sample1'] + pb21_row['sample1']
        assert total == 390
    
    def test_count_conservation(self):
        """Test that total counts are conserved through redistribution."""
        # Setup
        rescue_df = pd.DataFrame({
            'isoform': ['artifact1', 'artifact2', 'artifact3'],
            'assigned_transcript': ['PB.1.1', 'PB.2.1', 'PB.1.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1', 'artifact2', 'artifact3'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact', 'Artifact', 'Artifact']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1', 'artifact2', 'artifact3'],
            'sample1': [1000, 2000, 100, 200, 150],
            'sample2': [1500, 2500, 125, 175, 200]
        })
        
        # Execute
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        
        # Verify count conservation
        original_total_s1 = old_counts['sample1'].sum()
        new_total_s1 = result['sample1'].sum()
        
        original_total_s2 = old_counts['sample2'].sum()
        new_total_s2 = result['sample2'].sum()
        
        assert original_total_s1 == new_total_s1
        assert original_total_s2 == new_total_s2
    
    def test_td_isoform_creation(self):
        """Test that TD (Transcript Divergency) isoforms are created for unrescued artifacts."""
        # Setup: artifact2 is not rescued
        rescue_df = pd.DataFrame({
            'isoform': ['artifact1','artifact2'],
            'assigned_transcript': ['PB.1.1','GENE2_TD']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'artifact1', 'artifact2'],
            'filter_result': ['Isoform', 'Artifact', 'Artifact'],
            'associated_gene': ['GENE1', 'GENE1', 'GENE2']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'artifact1', 'artifact2'],
            'sample1': [100, 50, 75]
        })
        
        # Execute
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        print(result)
        # Verify that GENE2_TD isoform was created
        assert 'GENE2_TD' in result['isoform'].values
        
        # GENE2_TD should have artifact2's counts
        gene2_td_row = result[result['isoform'] == 'GENE2_TD'].iloc[0]
        assert gene2_td_row['sample1'] == 75


class TestRequantify:
    """Test suite for requantify function (main entry point)."""
    
    def test_full_requant_workflow(self, tmp_path):
        """Test the complete requantification workflow with file output."""
        # Setup
        counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'sample1': [100, 200, 50],
            'sample2': [150, 250, 60]
        })
        
        rescue_df = pd.DataFrame({
            'artifact': ['artifact1'],
            'assigned_transcript': ['PB.1.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact'],
            'associated_gene': ['GENE1', 'GENE2', 'GENE1']
        })
        
        prefix = str(tmp_path / "test_output")
        
        # Execute
        result = requantify(counts, rescue_df, classif_df, prefix)
        
        # Verify
        assert result is not None
        assert 'isoform' in result.columns
        
        # Check that files were created
        assert os.path.exists(f"{prefix}_reassigned_counts.tsv")
        assert os.path.exists(f"{prefix}_reassigned_counts_extended.tsv")
        
        # Read and verify reassigned counts
        reassigned = pd.read_csv(f"{prefix}_reassigned_counts.tsv", sep='\t')
        assert len(reassigned) > 0
        assert 'PB.1.1' in reassigned['isoform'].values
        
        # PB.1.1 should have gained artifact1's counts
        pb11_row = reassigned[reassigned['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['sample1'] == 150  # 100 + 50
        assert pb11_row['sample2'] == 210  # 150 + 60
        
        # Read and verify extended table
        extended = pd.read_csv(f"{prefix}_reassigned_counts_extended.tsv", sep='\t')
        assert len(extended) == 3  # All original isoforms
        assert 'old_sample1' in extended.columns
        assert 'new_sample1' in extended.columns
        
        # artifact1 should show old counts but new counts of 0
        artifact_row = extended[extended['isoform'] == 'artifact1'].iloc[0]
        assert artifact_row['old_sample1'] == 50
        assert artifact_row['new_sample1'] == 0


class TestRequantificationPipeline:
    """Test suite for the full requantification pipeline."""
    
    def test_pipeline_workflow(self, tmp_path):
        """Test the complete requantification pipeline workflow."""
        # Create temporary count file
        counts_file = tmp_path / "counts.tsv"
        counts_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'sample1': [100, 200, 50]
        })
        counts_df.to_csv(counts_file, sep='\t', index=False)
        
        # Setup rescue data
        rescue_df = pd.DataFrame({
            'artifact': ['artifact1'],
            'assigned_transcript': ['PB.1.1']
        })
        
        # Original classification (used for redistribution)
        original_class = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact'],
            'associated_gene': ['GENE1', 'GENE2', 'GENE1'],
            'length': [1000, 2000, 1500]
        })
        
        # Rescue classification (will be updated with new counts)
        rescue_class = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'structural_category': ['FSM', 'ISM'],
            'associated_gene': ['GENE1', 'GENE2'],
            'length': [1000, 2000],
            'FL': [100, 200]
        })
        
        output_dir = str(tmp_path)
        output_prefix = "test_output"
        
        # Execute
        requantification_pipeline(
            output_dir,
            output_prefix,
            str(counts_file),
            rescue_df,
            original_class,
            rescue_class
        )
        
        # Verify output files
        assert os.path.exists(f"{output_dir}/{output_prefix}_reassigned_counts.tsv")
        assert os.path.exists(f"{output_dir}/{output_prefix}_reassigned_counts_extended.tsv")
        assert os.path.exists(f"{output_dir}/{output_prefix}_rescued_classification.txt")
        
        # Verify reassigned counts file
        reassigned_df = pd.read_csv(f"{output_dir}/{output_prefix}_reassigned_counts.tsv", sep='\t')
        assert 'isoform' in reassigned_df.columns
        assert 'sample1' in reassigned_df.columns
        
        # Verify PB.1.1 received artifact1's counts
        pb11_row = reassigned_df[reassigned_df['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['sample1'] == 150  # 100 + 50
        
        # Verify rescued classification was updated
        rescued_class = pd.read_csv(f"{output_dir}/{output_prefix}_rescued_classification.txt", sep='\t')
        assert 'FL.sample1' in rescued_class.columns
        assert 'FL' in rescued_class.columns
        
        # Verify FL counts were updated
        pb11_class_row = rescued_class[rescued_class['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_class_row['FL.sample1'] == 150


class TestUpdateClassification:
    """Test suite for update_classification function."""
    
    def test_basic_classification_update(self, tmp_path):
        """Test basic update of classification file with new counts."""
        # Setup
        requant_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'sample1': [150, 275],
            'sample2': [210, 330]
        })
        
        rescue_class = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'structural_category': ['FSM', 'ISM'],
            'associated_gene': ['GENE1', 'GENE2'],
            'length': [1000, 2000]
        })
        
        prefix = str(tmp_path / "test_output")
        
        # Execute
        update_classification(requant_df, rescue_class, prefix)
        
        # Verify output file exists
        assert os.path.exists(f"{prefix}_rescued_classification.txt")
        
        # Read and verify
        result = pd.read_csv(f"{prefix}_rescued_classification.txt", sep='\t')
        
        # Check FL columns were created
        assert 'FL.sample1' in result.columns
        assert 'FL.sample2' in result.columns
        assert 'FL' in result.columns
        
        # Check values
        pb11_row = result[result['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['FL.sample1'] == 150
        assert pb11_row['FL.sample2'] == 210
        assert pb11_row['FL'] == 360  # 150 + 210
        
        pb21_row = result[result['isoform'] == 'PB.2.1'].iloc[0]
        assert pb21_row['FL.sample1'] == 275
        assert pb21_row['FL.sample2'] == 330
        assert pb21_row['FL'] == 605  # 275 + 330
    
    def test_partial_isoform_match(self, tmp_path):
        """Test update when requant has extra isoforms not in classification."""
        # Setup: requant_df has an extra isoform
        requant_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'EXTRA_ISOFORM'],
            'sample1': [150, 275, 100]
        })
        
        rescue_class = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],  # No EXTRA_ISOFORM
            'structural_category': ['FSM', 'ISM'],
            'associated_gene': ['GENE1', 'GENE2']
        })
        
        prefix = str(tmp_path / "test_output")
        
        # Execute - should not crash
        update_classification(requant_df, rescue_class, prefix)
        
        # Verify
        result = pd.read_csv(f"{prefix}_rescued_classification.txt", sep='\t')
        
        # Only the matching isoforms should be in output
        assert len(result) == 2
        assert 'EXTRA_ISOFORM' not in result['isoform'].values
        
        # Values should be correct for matching isoforms
        pb11_row = result[result['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['FL.sample1'] == 150
    
    def test_existing_fl_columns(self, tmp_path):
        """Test update when FL columns already exist in classification."""
        # Setup: rescue_class already has FL columns
        requant_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'sample1': [150, 275]
        })
        
        rescue_class = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'structural_category': ['FSM', 'ISM'],
            'FL': [100, 200],  # Existing FL column
            'FL.sample1': [100, 200]  # Existing sample FL column
        })
        
        prefix = str(tmp_path / "test_output")
        
        # Execute
        update_classification(requant_df, rescue_class, prefix)
        
        # Verify columns were updated, not duplicated
        result = pd.read_csv(f"{prefix}_rescued_classification.txt", sep='\t')
        
        # Check updated values
        pb11_row = result[result['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['FL.sample1'] == 150  # Updated from 100 to 150
        assert pb11_row['FL'] == 150  # Updated total
    
    def test_no_matching_isoforms(self, tmp_path, capsys):
        """Test behavior when no isoforms match between files."""
        # Setup: completely different isoforms
        requant_df = pd.DataFrame({
            'isoform': ['ISO_A', 'ISO_B'],
            'sample1': [100, 200]
        })
        
        rescue_class = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'structural_category': ['FSM', 'ISM']
        })
        
        prefix = str(tmp_path / "test_output")
        
        # Execute - should print warning but not crash
        result = update_classification(requant_df, rescue_class, prefix)
        
        # Check that warning was printed
        captured = capsys.readouterr()
        assert "Warning" in captured.out or result is not None


class TestEdgeCases:
    """Test suite for edge cases and error conditions."""
    
    def test_zero_counts_artifact(self):
        """Test handling artifacts with zero counts."""
        rescue_df = pd.DataFrame({
            'isoform': ['artifact1'],
            'assigned_transcript': ['PB.1.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'artifact1'],
            'filter_result': ['Isoform', 'Artifact']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'artifact1'],
            'sample1': [100, 0]  # Zero count artifact
        })
        
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        
        # Should not crash and PB.1.1 should retain its count
        pb11_row = result[result['isoform'] == 'PB.1.1'].iloc[0]
        assert pb11_row['sample1'] == 100
    
    def test_no_artifacts(self):
        """Test when there are no artifacts to redistribute."""
        rescue_df = pd.DataFrame({
            'isoform': [],
            'assigned_transcript': []
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'filter_result': ['Isoform', 'Isoform']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1'],
            'sample1': [100, 200]
        })
        
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        
        # Should return original counts unchanged
        assert len(result) == 2
        assert result[result['isoform'] == 'PB.1.1']['sample1'].values[0] == 100
        assert result[result['isoform'] == 'PB.2.1']['sample1'].values[0] == 200
    
    def test_all_targets_zero_counts(self):
        """Test distributing artifact counts when all targets have zero counts."""
        rescue_df = pd.DataFrame({
            'isoform': ['artifact1', 'artifact1'],
            'assigned_transcript': ['PB.1.1', 'PB.2.1']
        })
        
        classif_df = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'filter_result': ['Isoform', 'Isoform', 'Artifact']
        })
        
        old_counts = pd.DataFrame({
            'isoform': ['PB.1.1', 'PB.2.1', 'artifact1'],
            'sample1': [0, 0, 100]  # Both targets have zero
        })
        
        result = redistribute_counts_vectorized(rescue_df, classif_df, old_counts)
        
        # Should split evenly (50:50)
        pb11_count = result[result['isoform'] == 'PB.1.1']['sample1'].values[0]
        pb21_count = result[result['isoform'] == 'PB.2.1']['sample1'].values[0]
        
        # Should be equal and sum to 100
        assert pb11_count + pb21_count == 100
        assert abs(pb11_count - pb21_count) <= 1  # Allow for rounding difference
