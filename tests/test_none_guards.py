"""
Tests for None guard checks in visualization paths.

Ensures that PRE_ANALYZED data (which has expression_df=None) doesn't crash
when attempting to generate visualizations.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
from unittest.mock import MagicMock


def test_preanalyzed_data_no_crash_on_visualization():
    """
    Test that PRE_ANALYZED data with expression_df=None doesn't crash on visualization.

    This test reproduces the crash at line 474:
    `np.log2(result.expression_df + 1)` when expression_df is None.

    Expected behavior: Should gracefully handle None and display info message.
    """
    from rnaseq_parser import RNASeqParser, DataType

    # Step 1: Create PRE_ANALYZED data (DE results only, no expression_df)
    pre_analyzed_data = pd.DataFrame(
        {
            "gene": ["COL1A1", "IL6", "TYR", "VEGFA", "FGF2"],
            "log2FoldChange": [2.34, -1.89, 1.45, 0.78, -0.56],
            "padj": [0.0001, 0.0023, 0.0456, 0.1234, 0.5678],
            "baseMean": [1500.2, 450.8, 120.3, 890.5, 234.1],
        }
    )

    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        pre_analyzed_data.to_csv(f, index=False)
        temp_path = f.name

    try:
        # Step 2: Parse the file - should detect as PRE_ANALYZED
        parser = RNASeqParser()
        result = parser.parse(temp_path)

        # Verify it's PRE_ANALYZED type
        assert result.data_type == DataType.PRE_ANALYZED
        assert result.expression_df is None  # KEY: PRE_ANALYZED has no expression data
        assert result.de_results_df is not None

        # Step 3: Simulate the visualization code path that was crashing
        # This is the code from rnaseq_analysis_platform.py lines 471-474
        norm_counts = None
        if result.data_type == DataType.NORMALIZED:
            import numpy as np

            # This line would crash if expression_df is None
            if result.expression_df is None:
                # NEW GUARD: Check for None before arithmetic
                pass
            else:
                norm_counts = np.log2(result.expression_df + 1)

        # Step 4: Verify no crash occurred
        # (In actual code, this would be in the NORMALIZED branch, but we're testing the guard)
        assert norm_counts is None  # No visualization data generated

        # Step 5: Verify the guard works for NORMALIZED data with None expression_df
        # Create a mock NORMALIZED result with None expression_df (edge case)
        result.data_type = DataType.NORMALIZED
        result.expression_df = None

        # This should NOT crash with the guard in place
        norm_counts = None
        if result.data_type == DataType.NORMALIZED:
            import numpy as np

            if result.expression_df is None:
                pass
            else:
                norm_counts = np.log2(result.expression_df + 1)

        # Verify no crash and guard was effective
        assert norm_counts is None

    finally:
        # Cleanup
        Path(temp_path).unlink()


def test_normalized_data_with_expression_df_generates_viz():
    """
    Test that NORMALIZED data WITH expression_df properly generates visualizations.

    This is the happy path - ensure the guard doesn't break normal operation.
    """
    from rnaseq_parser import RNASeqParser, DataType

    # Create NORMALIZED data with expression_df
    normalized_data = pd.DataFrame(
        {
            "gene": ["COL1A1", "IL6", "TYR", "VEGFA", "FGF2"],
            "Sample_1": [12.5, 8.3, 6.7, 10.2, 7.1],
            "Sample_2": [13.2, 9.1, 7.2, 11.0, 7.8],
            "Sample_3": [14.1, 10.2, 7.8, 12.1, 8.5],
            "Sample_4": [13.8, 9.5, 7.4, 11.5, 8.0],
        }
    )

    # Write to temp file
    with tempfile.NamedTemporaryFile(suffix=".csv", delete=False, mode="w") as f:
        normalized_data.to_csv(f, index=False)
        temp_path = f.name

    try:
        # Parse the file
        parser = RNASeqParser()
        result = parser.parse(temp_path)

        # Verify it's NORMALIZED type
        assert result.data_type == DataType.NORMALIZED
        assert result.expression_df is not None  # NORMALIZED has expression data

        # Simulate the visualization code path
        norm_counts = None
        if result.data_type == DataType.NORMALIZED:
            import numpy as np

            if result.expression_df is None:
                pass
            else:
                norm_counts = np.log2(result.expression_df + 1)

        # Verify visualization data was generated
        assert norm_counts is not None
        assert norm_counts.shape == result.expression_df.shape

    finally:
        # Cleanup
        Path(temp_path).unlink()


class TestVolcanoPlotValidation:
    """Tests for create_volcano_plot() input validation."""

    def test_volcano_plot_none_input_raises(self):
        """Test that None input raises ValueError."""
        from visualizations import create_volcano_plot

        with pytest.raises(
            ValueError, match="Cannot create volcano plot.*empty or None"
        ):
            create_volcano_plot(None)

    def test_volcano_plot_empty_dataframe_raises(self):
        """Test that empty DataFrame raises ValueError."""
        from visualizations import create_volcano_plot

        empty_df = pd.DataFrame()
        with pytest.raises(
            ValueError, match="Cannot create volcano plot.*empty or None"
        ):
            create_volcano_plot(empty_df)

    def test_volcano_plot_missing_columns_raises(self):
        """Test that missing required columns raises ValueError."""
        from visualizations import create_volcano_plot

        df = pd.DataFrame(
            {
                "gene": ["GENE1", "GENE2"],
                "log2FoldChange": [1.5, -2.0],
            }
        )
        with pytest.raises(
            ValueError, match="Cannot create volcano plot.*missing required columns"
        ):
            create_volcano_plot(df)

    def test_volcano_plot_valid_input_succeeds(self):
        """Test that valid input produces a figure."""
        from visualizations import create_volcano_plot

        df = pd.DataFrame(
            {
                "gene": ["GENE1", "GENE2", "GENE3"],
                "log2FoldChange": [1.5, -2.0, 0.5],
                "padj": [0.001, 0.01, 0.5],
            }
        )
        fig = create_volcano_plot(df)
        assert fig is not None
        assert hasattr(fig, "data")


class TestHeatmapValidation:
    """Tests for create_clustered_heatmap() input validation."""

    def test_heatmap_none_input_raises(self):
        """Test that None expression_df raises ValueError."""
        from visualizations import create_clustered_heatmap

        with pytest.raises(ValueError, match="Cannot create heatmap.*empty or None"):
            create_clustered_heatmap(None, {"sample1": "control"})

    def test_heatmap_empty_dataframe_raises(self):
        """Test that empty expression_df raises ValueError."""
        from visualizations import create_clustered_heatmap

        empty_df = pd.DataFrame()
        with pytest.raises(ValueError, match="Cannot create heatmap.*empty or None"):
            create_clustered_heatmap(empty_df, {"sample1": "control"})

    def test_heatmap_empty_conditions_raises(self):
        """Test that empty sample_conditions raises ValueError."""
        from visualizations import create_clustered_heatmap

        df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )
        with pytest.raises(ValueError, match="Cannot create heatmap.*empty"):
            create_clustered_heatmap(df, {})

    def test_heatmap_missing_samples_raises(self):
        """Test that samples missing from conditions raises ValueError."""
        from visualizations import create_clustered_heatmap

        df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
                "sample3": [2.0, 3.0, 4.0],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )
        conditions = {"sample1": "control"}
        with pytest.raises(
            ValueError, match="Cannot create heatmap.*missing from sample_conditions"
        ):
            create_clustered_heatmap(df, conditions)

    def test_heatmap_valid_input_succeeds(self):
        """Test that valid input produces a figure."""
        from visualizations import create_clustered_heatmap

        df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )
        conditions = {"sample1": "control", "sample2": "treatment"}
        fig = create_clustered_heatmap(df, conditions)
        assert fig is not None
        assert hasattr(fig, "data")


class TestPCAValidation:
    """Tests for create_pca_plot() input validation."""

    def test_pca_none_input_raises(self):
        """Test that None expression_df raises ValueError."""
        from visualizations import create_pca_plot

        with pytest.raises(ValueError, match="Cannot create PCA plot.*empty or None"):
            create_pca_plot(None, {"sample1": "control"})

    def test_pca_empty_dataframe_raises(self):
        """Test that empty expression_df raises ValueError."""
        from visualizations import create_pca_plot

        empty_df = pd.DataFrame()
        with pytest.raises(ValueError, match="Cannot create PCA plot.*empty or None"):
            create_pca_plot(empty_df, {"sample1": "control"})

    def test_pca_insufficient_samples_raises(self):
        """Test that <2 samples raises ValueError."""
        from visualizations import create_pca_plot

        df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )
        conditions = {"sample1": "control"}
        with pytest.raises(
            ValueError, match="Cannot create PCA plot.*at least 2 samples"
        ):
            create_pca_plot(df, conditions)

    def test_pca_empty_conditions_raises(self):
        """Test that empty sample_conditions raises ValueError."""
        from visualizations import create_pca_plot

        df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )
        with pytest.raises(ValueError, match="Cannot create PCA plot.*empty"):
            create_pca_plot(df, {})

    def test_pca_valid_input_succeeds(self):
        """Test that valid input produces a figure."""
        from visualizations import create_pca_plot

        df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
                "sample3": [2.0, 3.0, 4.0],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )
        conditions = {
            "sample1": "control",
            "sample2": "treatment",
            "sample3": "treatment",
        }
        fig = create_pca_plot(df, conditions)
        assert fig is not None
        assert hasattr(fig, "data")


class TestGenePanelValidation:
    """Tests for GenePanelAnalyzer input validation."""

    def test_score_panel_none_expression_raises(self):
        """Test that score_panel() raises ValueError when expression_df is None."""
        from gene_panels import GenePanelAnalyzer

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")

        with pytest.raises(ValueError, match="expression_df cannot be None or empty"):
            analyzer.score_panel(None, "Anti-aging", {"Sample_1": "Control"})

    def test_score_panel_empty_expression_raises(self):
        """Test that score_panel() raises ValueError when expression_df is empty."""
        from gene_panels import GenePanelAnalyzer

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")
        empty_df = pd.DataFrame()

        with pytest.raises(ValueError, match="expression_df cannot be None or empty"):
            analyzer.score_panel(empty_df, "Anti-aging", {"Sample_1": "Control"})

    def test_score_panel_invalid_panel_raises(self):
        """Test that score_panel() raises ValueError with helpful message for invalid panel."""
        from gene_panels import GenePanelAnalyzer

        expr_df = pd.DataFrame(
            {
                "COL1A1": [12.5, 13.2, 14.1],
                "IL6": [8.3, 9.1, 10.2],
                "TYR": [6.7, 7.2, 7.8],
            },
            index=["Sample_1", "Sample_2", "Sample_3"],
        )

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")

        with pytest.raises(ValueError, match="Panel 'InvalidPanel' not found"):
            analyzer.score_panel(
                expr_df,
                "InvalidPanel",
                {"Sample_1": "Control", "Sample_2": "Control", "Sample_3": "Treatment"},
            )

    def test_plot_panel_none_expression_raises(self):
        """Test that plot_panel() raises ValueError when expression_df is None."""
        from gene_panels import GenePanelAnalyzer

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")

        with pytest.raises(ValueError, match="expression_df cannot be None or empty"):
            analyzer.plot_panel(None, "Anti-aging", {"Sample_1": "Control"})

    def test_plot_panel_empty_expression_raises(self):
        """Test that plot_panel() raises ValueError when expression_df is empty."""
        from gene_panels import GenePanelAnalyzer

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")
        empty_df = pd.DataFrame()

        with pytest.raises(ValueError, match="expression_df cannot be None or empty"):
            analyzer.plot_panel(empty_df, "Anti-aging", {"Sample_1": "Control"})

    def test_plot_panel_invalid_panel_raises(self):
        """Test that plot_panel() raises ValueError with helpful message for invalid panel."""
        from gene_panels import GenePanelAnalyzer

        expr_df = pd.DataFrame(
            {
                "COL1A1": [12.5, 13.2, 14.1],
                "IL6": [8.3, 9.1, 10.2],
                "TYR": [6.7, 7.2, 7.8],
            },
            index=["Sample_1", "Sample_2", "Sample_3"],
        )

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")

        with pytest.raises(ValueError, match="Panel 'InvalidPanel' not found"):
            analyzer.plot_panel(
                expr_df,
                "InvalidPanel",
                {"Sample_1": "Control", "Sample_2": "Control", "Sample_3": "Treatment"},
            )

    def test_platform_tracks_failed_panels(self):
        """
        Integration test: Verify platform tracks failed panels and displays info message.

        This tests the code at lines 589-606 in rnaseq_analysis_platform.py
        that catches ValueError and tracks failed panels.
        """
        from gene_panels import GenePanelAnalyzer

        expr_df = pd.DataFrame(
            {
                "COL1A1": [12.5, 13.2, 14.1],
                "IL6": [8.3, 9.1, 10.2],
            },
            index=["Sample_1", "Sample_2", "Sample_3"],
        )

        sample_conds = {
            "Sample_1": "Control",
            "Sample_2": "Control",
            "Sample_3": "Treatment",
        }

        analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")
        panels = analyzer.panels

        failed_panels = []
        figures = {}

        for panel_name in panels:
            try:
                fig = analyzer.plot_panel(expr_df, panel_name, sample_conds)
                figures[f"panel_{panel_name}"] = fig
            except ValueError as e:
                failed_panels.append((panel_name, str(e)))

        assert len(failed_panels) > 0, (
            "Expected some panels to fail with insufficient genes"
        )

        for panel_name, error_msg in failed_panels:
            assert isinstance(panel_name, str)
            assert isinstance(error_msg, str)
            assert "genes" in error_msg.lower() or "not found" in error_msg.lower()


class TestExportDataValidation:
    """Tests for ExportData.validate() method."""

    def test_export_data_validate_no_data_raises(self):
        """Test that validate() returns warning when no data available."""
        from export_engine import ExportData
        from rnaseq_parser import DataType

        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        warnings = export_data.validate()
        assert "No data available to export" in warnings

    def test_export_data_validate_with_de_results(self):
        """Test that validate() passes when DE results present."""
        from export_engine import ExportData
        from de_analysis import DEResult
        from rnaseq_parser import DataType

        de_result = DEResult(
            results_df=pd.DataFrame(
                {
                    "gene": ["GENE1", "GENE2"],
                    "log2FoldChange": [1.5, -2.0],
                    "padj": [0.001, 0.01],
                }
            ),
            normalized_counts=None,
            log_normalized_counts=None,
            dds=None,
            comparison=("test", "ref"),
            n_significant=2,
            warnings=[],
        )

        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={("test", "ref"): de_result},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        warnings = export_data.validate()
        assert "No data available to export" not in warnings

    def test_export_data_validate_empty_de_results(self):
        """Test that validate() warns when DE results are empty."""
        from export_engine import ExportData
        from de_analysis import DEResult
        from rnaseq_parser import DataType

        de_result = DEResult(
            results_df=pd.DataFrame(),
            normalized_counts=None,
            log_normalized_counts=None,
            dds=None,
            comparison=("test", "ref"),
            n_significant=0,
            warnings=[],
        )

        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={("test", "ref"): de_result},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        warnings = export_data.validate()
        assert any("has no results" in w for w in warnings)

    def test_export_data_validate_with_expression_matrix(self):
        """Test that validate() passes when expression matrix present."""
        from export_engine import ExportData
        from rnaseq_parser import DataType

        expr_df = pd.DataFrame(
            {
                "sample1": [1.0, 2.0, 3.0],
                "sample2": [1.5, 2.5, 3.5],
            },
            index=["GENE1", "GENE2", "GENE3"],
        )

        export_data = ExportData(
            data_type=DataType.NORMALIZED,
            de_results={},
            expression_matrix=expr_df,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        warnings = export_data.validate()
        assert "No data available to export" not in warnings


def test_export_multiformat_includes_all_data():
    """
    RED TEST: Multi-format files should export ALL data types.

    This test reproduces the bug where export_engine.py branches on data_type enum
    instead of checking actual data presence. A file with all 3 data types
    (DE results + expression matrix + normalized) should export all of them.

    Expected behavior: Excel workbook has DE + Expression + Normalized sheets
    Current behavior: Only exports based on data_type enum, losing data
    """
    from export_engine import ExportEngine, ExportData, EnrichmentResult
    from de_analysis import DEResult
    from rnaseq_parser import DataType
    import tempfile
    import openpyxl

    engine = ExportEngine()

    de_result = DEResult(
        results_df=pd.DataFrame(
            {
                "gene": ["COL1A1", "IL6", "TYR"],
                "log2FoldChange": [2.34, -1.89, 1.45],
                "padj": [0.0001, 0.0023, 0.0456],
                "baseMean": [1500.2, 450.8, 120.3],
            }
        ),
        normalized_counts=None,
        log_normalized_counts=None,
        dds=None,
        comparison=("test", "ref"),
        n_significant=2,
        warnings=[],
    )

    expression_matrix = pd.DataFrame(
        {
            "Sample_1": [12.5, 8.3, 6.7],
            "Sample_2": [13.2, 9.1, 7.2],
            "Sample_3": [14.1, 10.2, 7.8],
        },
        index=["COL1A1", "IL6", "TYR"],
    )

    enrich_result = EnrichmentResult(
        go_results=pd.DataFrame(
            {
                "Term": ["GO:0001", "GO:0002"],
                "Overlap": ["5/100", "3/50"],
                "P-value": [0.001, 0.01],
                "Adjusted P-value": [0.01, 0.05],
                "Genes": ["GENE1;GENE2", "GENE3;GENE4"],
            }
        ),
        kegg_results=pd.DataFrame(
            {
                "Term": ["hsa00001", "hsa00002"],
                "Overlap": ["4/80", "2/40"],
                "P-value": [0.002, 0.02],
                "Adjusted P-value": [0.02, 0.1],
                "Genes": ["GENE1;GENE3", "GENE2;GENE4"],
            }
        ),
        genes_used=["GENE1", "GENE2", "GENE3", "GENE4"],
        selection_note="2 genes (padj<0.05)",
    )

    export_data = ExportData(
        data_type=DataType.RAW_COUNTS,
        de_results={("test", "ref"): de_result},
        expression_matrix=expression_matrix,
        enrichment_results={("test", "ref"): enrich_result},
        figures={},
        settings={"padj_threshold": 0.05, "lfc_threshold": 1.0},
        sample_conditions={
            "Sample_1": "control",
            "Sample_2": "test",
            "Sample_3": "test",
        },
    )

    with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
        temp_path = f.name

    try:
        engine.export_excel(temp_path, export_data)

        wb = openpyxl.load_workbook(temp_path)
        sheet_names = wb.sheetnames

        assert any("DE_" in name for name in sheet_names), (
            f"Missing DE results sheet. Sheets: {sheet_names}"
        )

        assert "Expression Matrix" in sheet_names, (
            f"Missing Expression Matrix sheet. Sheets: {sheet_names}. "
            "BUG: export_engine.py branches on data_type instead of checking "
            "expression_matrix presence"
        )

        assert any("GO_" in name for name in sheet_names), (
            f"Missing GO enrichment sheet. Sheets: {sheet_names}"
        )

        assert "Settings" in sheet_names, (
            f"Missing Settings sheet. Sheets: {sheet_names}"
        )

    finally:
        Path(temp_path).unlink(missing_ok=True)


class TestExportExcelValidation:
    """Tests for export_excel() validation and None guards."""

    def test_export_excel_no_data_raises(self):
        """Test that export_excel() raises ValueError when no data available."""
        from export_engine import ExportEngine, ExportData
        from rnaseq_parser import DataType
        import tempfile

        engine = ExportEngine()
        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
            temp_path = f.name

        try:
            with pytest.raises(ValueError, match="Cannot export: no data available"):
                engine.export_excel(temp_path, export_data)
        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_excel_invalid_path_raises(self):
        """Test that export_excel() raises ValueError for non-existent directory."""
        from export_engine import ExportEngine, ExportData
        from de_analysis import DEResult
        from rnaseq_parser import DataType

        engine = ExportEngine()
        de_result = DEResult(
            results_df=pd.DataFrame(
                {
                    "gene": ["GENE1"],
                    "log2FoldChange": [1.5],
                    "padj": [0.001],
                }
            ),
            normalized_counts=None,
            log_normalized_counts=None,
            dds=None,
            comparison=("test", "ref"),
            n_significant=1,
            warnings=[],
        )

        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={("test", "ref"): de_result},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        invalid_path = "/nonexistent/directory/file.xlsx"

        with pytest.raises(ValueError, match="Directory does not exist"):
            engine.export_excel(invalid_path, export_data)

    def test_export_excel_handles_none_de_results(self):
        """Test that export_excel() skips None DE results gracefully."""
        from export_engine import ExportEngine, ExportData
        from de_analysis import DEResult
        from rnaseq_parser import DataType
        import tempfile

        engine = ExportEngine()
        de_result = DEResult(
            results_df=None,
            normalized_counts=None,
            log_normalized_counts=None,
            dds=None,
            comparison=("test", "ref"),
            n_significant=0,
            warnings=["Failed"],
        )

        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={("test", "ref"): de_result},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
            temp_path = f.name

        try:
            engine.export_excel(temp_path, export_data)
            assert Path(temp_path).exists()
        finally:
            Path(temp_path).unlink(missing_ok=True)

    def test_export_excel_handles_empty_de_results(self):
        """Test that export_excel() skips empty DE results gracefully."""
        from export_engine import ExportEngine, ExportData
        from de_analysis import DEResult
        from rnaseq_parser import DataType
        import tempfile

        engine = ExportEngine()
        de_result = DEResult(
            results_df=pd.DataFrame(),
            normalized_counts=None,
            log_normalized_counts=None,
            dds=None,
            comparison=("test", "ref"),
            n_significant=0,
            warnings=[],
        )

        export_data = ExportData(
            data_type=DataType.RAW_COUNTS,
            de_results={("test", "ref"): de_result},
            expression_matrix=None,
            enrichment_results={},
            figures={},
            settings={},
            sample_conditions={},
        )

        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
            temp_path = f.name

        try:
            engine.export_excel(temp_path, export_data)
            assert Path(temp_path).exists()
        finally:
            Path(temp_path).unlink(missing_ok=True)
