"""
Integration tests for RNA-seq Analysis Platform.

Tests the complete end-to-end workflows integrating all modules without Streamlit UI.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile


# =============================================================================
# TEST 1: RAW_COUNTS Full Workflow
# =============================================================================
def test_integration_raw_counts_full_workflow(test_data_dir, mock_gseapy):
    """
    Verify complete pipeline: upload RAW_COUNTS → metadata → DE → viz → export.

    This test runs the FULL analysis workflow programmatically (no Streamlit).
    """
    from rnaseq_parser import RNASeqParser, DataType
    from de_analysis import DEAnalysisEngine
    from visualizations import (
        create_volcano_plot,
        create_clustered_heatmap,
        create_pca_plot,
    )
    from export_engine import ExportEngine, ExportData

    # Step 1: Parse RAW_COUNTS file
    parser = RNASeqParser()
    result = parser.parse(test_data_dir / "sample_counts.csv")

    assert result.data_type == DataType.RAW_COUNTS
    assert result.can_run_de == True
    assert result.expression_df is not None
    assert result.expression_df.shape == (10, 100)  # 10 samples × 100 genes

    # Step 2: Run DE analysis
    # NOTE: Both parser output AND DE engine input are samples × genes (NO TRANSPOSE)
    metadata = pd.DataFrame(
        {"condition": ["Control"] * 5 + ["Treatment"] * 5},
        index=[f"Sample_{i}" for i in range(1, 11)],
    )

    engine = DEAnalysisEngine()
    de_results = engine.run_all_comparisons(
        result.expression_df,  # samples × genes - same as parser output (NO TRANSPOSE)
        metadata,
        comparisons=[("Treatment", "Control")],
        design_factor="condition",
    )

    assert ("Treatment", "Control") in de_results
    de_result = de_results[("Treatment", "Control")]
    assert not de_result.results_df.empty
    assert "log2FoldChange" in de_result.results_df.columns
    assert "padj" in de_result.results_df.columns
    assert de_result.normalized_counts is not None  # RAW_COUNTS provides normalized

    # Step 3: Generate visualizations
    volcano = create_volcano_plot(de_result.results_df)
    assert volcano is not None  # Plotly Figure object

    # NOTE: Heatmap expects genes × samples, so we transpose log_normalized_counts
    sample_conditions = {
        f"Sample_{i}": "Control" if i <= 5 else "Treatment" for i in range(1, 11)
    }
    heatmap = create_clustered_heatmap(
        de_result.log_normalized_counts.T,  # Transpose: samples × genes → genes × samples for heatmap
        sample_conditions,
        de_results_df=de_result.results_df,
    )
    assert heatmap is not None

    # NOTE: PCA expects samples × genes (no transpose needed)
    pca = create_pca_plot(
        de_result.log_normalized_counts,  # samples × genes - no transpose
        sample_conditions,
    )
    assert pca is not None

    # Step 4: Export (verify no errors)
    export_data = ExportData(
        data_type=DataType.RAW_COUNTS,
        de_results=de_results,
        expression_matrix=de_result.log_normalized_counts,
        enrichment_results={},
        figures={"volcano": volcano, "heatmap": heatmap, "pca": pca},
        settings={"padj_threshold": 0.05, "lfc_threshold": 1.0},
        sample_conditions={
            f"Sample_{i}": "Control" if i <= 5 else "Treatment" for i in range(1, 11)
        },
        active_comparison=("Treatment", "Control"),
    )

    export_engine = ExportEngine()
    # Just verify construction works - actual file export tested in test_export.py
    assert export_data.de_results is not None


# =============================================================================
# TEST 2: PRE_ANALYZED Mode (Skip Metadata)
# =============================================================================
def test_integration_preanalyzed_mode_disables_metadata():
    """
    PRE_ANALYZED input should populate de_results_df, NOT expression_df.
    Metadata assignment and DE analysis should be skipped.
    """
    from rnaseq_parser import RNASeqParser, DataType

    # Create PRE_ANALYZED format file (DE results table)
    preanalyzed_data = pd.DataFrame(
        {
            "gene": ["GENE_1", "GENE_2", "GENE_3"],
            "log2FoldChange": [2.5, -1.8, 0.3],
            "pvalue": [0.001, 0.01, 0.5],
            "padj": [0.01, 0.05, 0.8],
        }
    )

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        preanalyzed_data.to_csv(f, index=False)
        temp_path = f.name

    parser = RNASeqParser()
    result = parser.parse(temp_path)

    # Verify PRE_ANALYZED detection
    assert result.data_type == DataType.PRE_ANALYZED
    assert result.can_run_de == False  # Cannot re-run DE
    assert result.expression_df is None  # No expression matrix
    assert result.de_results_df is not None  # Has DE results
    assert "log2FoldChange" in result.de_results_df.columns

    # Cleanup
    Path(temp_path).unlink()


# =============================================================================
# TEST 3: NORMALIZED Mode (Skip DE)
# =============================================================================
def test_integration_normalized_mode_disables_de(test_data_dir):
    """
    NORMALIZED input (floats) should disable DE analysis but enable viz.
    """
    from rnaseq_parser import RNASeqParser, DataType
    from visualizations import create_clustered_heatmap, create_pca_plot

    parser = RNASeqParser()
    result = parser.parse(test_data_dir / "normalized_data.csv")

    # Verify NORMALIZED detection
    assert result.data_type == DataType.NORMALIZED
    assert result.can_run_de == False  # DE disabled for normalized
    assert (
        result.expression_df is not None
    )  # Expression matrix present (samples × genes)
    assert result.de_results_df is None  # No pre-analyzed DE

    # Verify viz still works
    conditions = {
        f"Sample_{i}": "Control" if i <= 5 else "Treatment" for i in range(1, 11)
    }
    log_expr = np.log2(result.expression_df + 1)  # samples × genes

    # Heatmap: transpose to genes × samples
    heatmap = create_clustered_heatmap(log_expr.T, conditions)  # genes × samples
    assert heatmap is not None

    # PCA: samples × genes (no transpose)
    pca = create_pca_plot(log_expr, conditions)  # samples × genes
    assert pca is not None


# =============================================================================
# TEST 4: Multi-Comparison Storage
# =============================================================================
def test_integration_multi_comparison_storage(test_data_dir):
    """
    Multiple comparisons should be stored separately in de_results dict.
    """
    from rnaseq_parser import RNASeqParser
    from de_analysis import DEAnalysisEngine

    parser = RNASeqParser()
    result = parser.parse(test_data_dir / "sample_counts.csv")

    # Create 3-condition metadata
    metadata = pd.DataFrame(
        {"condition": ["Control"] * 3 + ["TreatmentA"] * 3 + ["TreatmentB"] * 4},
        index=[f"Sample_{i}" for i in range(1, 11)],
    )

    engine = DEAnalysisEngine()
    de_results = engine.run_all_comparisons(
        result.expression_df,  # samples × genes - NO TRANSPOSE (same as parser output)
        metadata,
        comparisons=[("TreatmentA", "Control"), ("TreatmentB", "Control")],
        design_factor="condition",
    )

    # Verify both comparisons stored
    assert len(de_results) == 2
    assert ("TreatmentA", "Control") in de_results
    assert ("TreatmentB", "Control") in de_results

    # Verify each has independent results
    result_a = de_results[("TreatmentA", "Control")]
    result_b = de_results[("TreatmentB", "Control")]

    assert not result_a.results_df.empty
    assert not result_b.results_df.empty
    # Results should differ (different comparisons)
    assert not result_a.results_df["log2FoldChange"].equals(
        result_b.results_df["log2FoldChange"]
    )


# =============================================================================
# TEST 5: Error Handling - Invalid Input
# =============================================================================
def test_integration_invalid_input_raises_error():
    """
    Invalid input (negative counts) should raise ParserValidationError.
    """
    from rnaseq_parser import RNASeqParser, ParserValidationError

    # Create invalid data with negative values
    invalid_data = pd.DataFrame(
        {
            "gene": ["GENE_1", "GENE_2"],
            "Sample_1": [100, -50],  # Negative value
            "Sample_2": [200, 150],
        }
    )

    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        invalid_data.to_csv(f, index=False)
        temp_path = f.name

    parser = RNASeqParser()

    with pytest.raises(ParserValidationError) as exc_info:
        parser.parse(temp_path)

    assert "negative" in exc_info.value.message.lower()

    # Cleanup
    Path(temp_path).unlink()
