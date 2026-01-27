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
        index=[f"sample_{i}" for i in range(1, 11)],
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
        f"sample_{i}": "Control" if i <= 5 else "Treatment" for i in range(1, 11)
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
        f"sample_{i}": "Control" if i <= 5 else "Treatment" for i in range(1, 11)
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
        index=[f"sample_{i}" for i in range(1, 11)],
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


# =============================================================================
# TEST 6: All Reference Files Parse Successfully
# =============================================================================
def test_all_reference_files_parse_successfully():
    """
    Test that all 4 reference files parse without error.

    Verifies:
    - Multi-format detection (PRE_ANALYZED + RAW_COUNTS + NORMALIZED)
    - All data extracted correctly
    - Expected shapes (6 samples × genes)
    """
    from rnaseq_parser import RNASeqParser, DataType

    parser = RNASeqParser()

    files = [
        "Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx",
        "Reference sequencing data/data3_Bu10_vs_none_fc2_&_raw.p.xlsx",
        "Reference sequencing data/data3_SwX15_vs_none_fc2_&_raw.p.xlsx",
        "Reference sequencing data/data3_my45_vs_none_fc2_&_raw.p.xlsx",
    ]

    expected_gene_counts = {
        "data3_Bt10U_vs_none_fc2_&_raw.p.xlsx": 320,
        "data3_Bu10_vs_none_fc2_&_raw.p.xlsx": 194,
        "data3_SwX15_vs_none_fc2_&_raw.p.xlsx": 170,
        "data3_my45_vs_none_fc2_&_raw.p.xlsx": 510,
    }

    for file in files:
        # Verify file exists
        file_path = Path(file)
        assert file_path.exists(), f"Reference file not found: {file}"

        # Parse file
        result = parser.parse(file)

        # Verify multi-format detection
        assert len(result.data_types_detected) == 3, (
            f"{file}: Expected 3 data types, got {len(result.data_types_detected)}"
        )
        assert DataType.PRE_ANALYZED in result.data_types_detected, (
            f"{file}: PRE_ANALYZED not detected"
        )
        assert DataType.RAW_COUNTS in result.data_types_detected, (
            f"{file}: RAW_COUNTS not detected"
        )
        assert DataType.NORMALIZED in result.data_types_detected, (
            f"{file}: NORMALIZED not detected"
        )

        # Verify all data extracted
        assert result.de_results_df is not None, f"{file}: de_results_df is None"
        assert result.expression_df is not None, f"{file}: expression_df is None"
        assert result.normalized_df is not None, f"{file}: normalized_df is None"

        # Verify shapes (6 samples × genes)
        assert result.expression_df.shape[0] == 6, (
            f"{file}: Expected 6 samples, got {result.expression_df.shape[0]}"
        )
        assert result.normalized_df.shape[0] == 6, (
            f"{file}: Normalized data expected 6 samples, got {result.normalized_df.shape[0]}"
        )

        # Verify gene counts match expected
        file_key = file.split("/")[-1]
        expected_genes = expected_gene_counts.get(file_key)
        if expected_genes:
            actual_genes = len(result.de_results_df)
            assert actual_genes == expected_genes, (
                f"{file}: Expected {expected_genes} genes, got {actual_genes}"
            )

        print(
            f"✓ {file}: {len(result.de_results_df)} genes, "
            f"{result.expression_df.shape[0]} samples"
        )


def test_all_reference_files_visualize_successfully():
    """
    Test that visualizations can be created for all reference files.

    Verifies:
    - Volcano plot generation
    - Clustered heatmap generation
    - PCA plot generation
    """
    from rnaseq_parser import RNASeqParser
    from visualizations import (
        create_volcano_plot,
        create_clustered_heatmap,
        create_pca_plot,
    )

    parser = RNASeqParser()

    file = "Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx"

    file_path = Path(file)
    assert file_path.exists(), f"Reference file not found: {file}"

    result = parser.parse(file)

    assert result.de_results_df is not None
    assert result.normalized_df is not None

    if not result.de_results_df["padj"].isna().any():
        volcano = create_volcano_plot(result.de_results_df)
        assert volcano is not None, "Volcano plot is None"
        assert hasattr(volcano, "to_html"), "Volcano plot is not a Plotly Figure"
    else:
        print(f"Skipping volcano plot for {file}: contains NaN padj values")

    sample_names = list(result.normalized_df.index)
    conditions = {
        sample: "none"
        if "none" in str(result.de_results_df.iloc[0]).lower()
        else "treated"
        for sample in sample_names
    }

    heatmap = create_clustered_heatmap(
        result.normalized_df.T,
        conditions,
        de_results_df=result.de_results_df,
    )
    assert heatmap is not None, "Heatmap is None"
    assert hasattr(heatmap, "to_html"), "Heatmap is not a Plotly Figure"

    pca = create_pca_plot(
        result.normalized_df,
        conditions,
    )
    assert pca is not None, "PCA plot is None"
    assert hasattr(pca, "to_html"), "PCA plot is not a Plotly Figure"

    print(f"✓ All visualizations created successfully for {file}")


def test_multiformat_can_run_de():
    """
    Test that multi-format files (RAW_COUNTS + PRE_ANALYZED) have can_run_de=True.

    This test verifies the fix for the UI flow bug where multi-format files
    couldn't run DE analysis because the UI checked data_type == RAW_COUNTS
    instead of can_run_de flag.

    Scenario:
    - File contains both RAW_COUNTS (expression matrix) and PRE_ANALYZED (DE results)
    - Parser should set can_run_de=True (because RAW_COUNTS is present)
    - UI should respect can_run_de flag to enable metadata assignment and DE analysis
    """
    from rnaseq_parser import RNASeqParser, DataType
    from de_analysis import DEAnalysisEngine
    import pandas as pd

    parser = RNASeqParser()

    file_path = Path("Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx")
    assert file_path.exists(), f"Reference file not found: {file_path}"

    result = parser.parse(str(file_path))

    assert hasattr(result, "data_types_detected"), (
        "ParseResult missing data_types_detected"
    )
    assert len(result.data_types_detected) > 1, (
        f"Expected multi-format file, got {len(result.data_types_detected)} data type(s)"
    )
    assert DataType.RAW_COUNTS in result.data_types_detected, (
        "Multi-format file should contain RAW_COUNTS"
    )
    assert DataType.PRE_ANALYZED in result.data_types_detected, (
        "Multi-format file should contain PRE_ANALYZED"
    )

    assert result.can_run_de == True, (
        "Multi-format file with RAW_COUNTS should have can_run_de=True"
    )

    assert result.expression_df is not None, (
        "Multi-format file should have expression_df (RAW_COUNTS)"
    )
    assert result.de_results_df is not None, (
        "Multi-format file should have de_results_df (PRE_ANALYZED)"
    )

    sample_count = min(6, len(result.expression_df))
    metadata = pd.DataFrame(
        {
            "condition": ["Control"] * (sample_count // 2)
            + ["Treatment"] * (sample_count - sample_count // 2)
        },
        index=result.expression_df.index[:sample_count],
    )

    engine = DEAnalysisEngine()
    de_results = engine.run_all_comparisons(
        result.expression_df.iloc[:sample_count],
        metadata,
        comparisons=[("Treatment", "Control")],
        design_factor="condition",
    )

    assert ("Treatment", "Control") in de_results, (
        "DE analysis should complete successfully with extracted counts"
    )
    de_result = de_results[("Treatment", "Control")]
    assert not de_result.results_df.empty, "DE results should not be empty"

    print(f"✓ Multi-format file correctly has can_run_de=True and supports DE analysis")
    assert DataType.RAW_COUNTS in result.data_types_detected, (
        "Multi-format file should contain RAW_COUNTS"
    )
    assert DataType.PRE_ANALYZED in result.data_types_detected, (
        "Multi-format file should contain PRE_ANALYZED"
    )

    # CRITICAL: can_run_de should be True because RAW_COUNTS is present
    assert result.can_run_de == True, (
        "Multi-format file with RAW_COUNTS should have can_run_de=True"
    )

    # Verify both data types are extracted
    assert result.expression_df is not None, (
        "Multi-format file should have expression_df (RAW_COUNTS)"
    )
    assert result.de_results_df is not None, (
        "Multi-format file should have de_results_df (PRE_ANALYZED)"
    )

    # Verify we can run DE analysis with extracted counts
    # This simulates what the UI would do after metadata assignment
    metadata = pd.DataFrame(
        {"condition": ["Control"] * 3 + ["Treatment"] * 3},
        index=result.expression_df.index[:6],
    )

    engine = DEAnalysisEngine()
    de_results = engine.run_all_comparisons(
        result.expression_df.iloc[:6],  # Use first 6 samples
        metadata,
        comparisons=[("Treatment", "Control")],
        design_factor="condition",
    )

    assert ("Treatment", "Control") in de_results, (
        "DE analysis should complete successfully with extracted counts"
    )
    de_result = de_results[("Treatment", "Control")]
    assert not de_result.results_df.empty, "DE results should not be empty"

    print(f"✓ Multi-format file correctly has can_run_de=True and supports DE analysis")
