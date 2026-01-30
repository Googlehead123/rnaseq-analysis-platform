"""
Tests for new PhD-level features:
- PathwayEnrichment extensions (GSEA, additional databases)
- InterpretationEngine
- CellTypeDeconvolution
- AdvancedQC
- Visualization extensions (enrichment dotplot, deconvolution plot)
- EnrichmentResult extension
"""

import pytest
import pandas as pd
import numpy as np
import plotly.graph_objects as go


# =============================================================================
# PathwayEnrichment Extensions
# =============================================================================

class TestPathwayEnrichmentExtensions:
    """Test extended pathway enrichment module."""

    def test_available_databases(self):
        from pathway_enrichment import PathwayEnrichment
        pe = PathwayEnrichment()
        dbs = pe.AVAILABLE_DATABASES
        assert len(dbs) == 7
        assert 'GO Biological Process' in dbs
        assert 'GO Molecular Function' in dbs
        assert 'GO Cellular Component' in dbs
        assert 'KEGG' in dbs
        assert 'Reactome' in dbs
        assert 'WikiPathway' in dbs
        assert 'MSigDB Hallmark' in dbs

    def test_create_ranking_metric(self):
        from pathway_enrichment import PathwayEnrichment
        pe = PathwayEnrichment()
        mock_de = pd.DataFrame({
            'gene': ['A', 'B', 'C'],
            'pvalue': [0.001, 0.05, 0.5],
            'log2FoldChange': [2.0, -1.5, 0.1],
            'padj': [0.01, 0.1, 0.8]
        })
        ranks = pe.create_ranking_metric(mock_de)
        assert len(ranks) == 3
        assert ranks.iloc[0] > 0  # Most positive (upregulated) first

    def test_create_ranking_metric_handles_nan(self):
        from pathway_enrichment import PathwayEnrichment
        pe = PathwayEnrichment()
        mock_de = pd.DataFrame({
            'gene': ['A', 'B', 'C'],
            'pvalue': [0.001, np.nan, 0.5],
            'log2FoldChange': [2.0, -1.5, np.nan],
            'padj': [0.01, 0.1, 0.8]
        })
        ranks = pe.create_ranking_metric(mock_de)
        assert len(ranks) == 1  # Only gene A has both valid pvalue and lfc

    def test_gsea_method_exists(self):
        from pathway_enrichment import PathwayEnrichment
        pe = PathwayEnrichment()
        assert hasattr(pe, 'run_gsea')
        assert callable(pe.run_gsea)

    def test_new_enrichment_methods_exist(self):
        from pathway_enrichment import PathwayEnrichment
        pe = PathwayEnrichment()
        for method in ['get_go_mf_enrichment', 'get_go_cc_enrichment',
                       'get_reactome_enrichment', 'get_wikipathway_enrichment',
                       'get_hallmark_enrichment']:
            assert hasattr(pe, method), f"Missing method: {method}"
            assert callable(getattr(pe, method))

    def test_gsea_empty_input(self):
        from pathway_enrichment import PathwayEnrichment
        pe = PathwayEnrichment()
        result, error = pe.run_gsea(pd.Series(dtype=float))
        assert error is not None
        assert result.empty


# =============================================================================
# InterpretationEngine
# =============================================================================

class TestInterpretationEngine:
    """Test interpretation engine module."""

    def test_interpret_de_results(self):
        from interpretation_engine import InterpretationEngine
        engine = InterpretationEngine()
        mock_de = pd.DataFrame({
            'gene': ['A', 'B', 'C', 'D'],
            'log2FoldChange': [3.0, -2.5, 0.1, 1.5],
            'padj': [0.001, 0.01, 0.8, 0.03],
            'pvalue': [0.0001, 0.001, 0.5, 0.01],
            'baseMean': [100, 200, 50, 150]
        })
        results = engine.interpret_de_results(mock_de, 'Treatment_vs_Control')
        assert len(results) >= 1
        assert results[0].category == 'de'
        assert results[0].level in ('critical', 'important', 'informative')
        assert len(results[0].description) > 0

    def test_interpret_de_no_significant(self):
        from interpretation_engine import InterpretationEngine
        engine = InterpretationEngine()
        mock_de = pd.DataFrame({
            'gene': ['A', 'B'],
            'log2FoldChange': [0.1, -0.2],
            'padj': [0.8, 0.9],
            'pvalue': [0.5, 0.6],
            'baseMean': [100, 200]
        })
        results = engine.interpret_de_results(mock_de, 'T_vs_C')
        # Should have "no significant genes" interpretation
        titles = [r.title for r in results]
        assert any('No Significant' in t for t in titles)

    def test_interpret_enrichment_empty(self):
        from interpretation_engine import InterpretationEngine
        engine = InterpretationEngine()
        results = engine.interpret_enrichment(pd.DataFrame())
        assert len(results) >= 1
        assert results[0].title == 'No Significant Enrichment'

    def test_interpret_enrichment_immune(self):
        from interpretation_engine import InterpretationEngine
        engine = InterpretationEngine()
        mock_enr = pd.DataFrame({
            'Term': ['immune response', 'T cell activation', 'metabolism'],
            'Adjusted P-value': [0.001, 0.01, 0.05]
        })
        results = engine.interpret_enrichment(mock_enr)
        titles = [r.title for r in results]
        assert any('Immune' in t or 'immune' in t.lower() for t in titles)

    def test_generate_methods_text(self):
        from interpretation_engine import InterpretationEngine
        engine = InterpretationEngine()
        text = engine.generate_methods_text({
            'de_tool': 'PyDESeq2',
            'padj': 0.05,
            'lfc': 1.0,
            'databases': ['GO Biological Process', 'KEGG']
        })
        assert 'PyDESeq2' in text
        assert '0.05' in text
        assert 'Methods' in text

    def test_modest_effect_size_detection(self):
        from interpretation_engine import InterpretationEngine
        engine = InterpretationEngine()
        # All significant genes have small LFC
        mock_de = pd.DataFrame({
            'gene': [f'G{i}' for i in range(20)],
            'log2FoldChange': [0.3] * 10 + [-0.3] * 10,
            'padj': [0.01] * 20,
            'pvalue': [0.001] * 20,
            'baseMean': [100] * 20
        })
        results = engine.interpret_de_results(mock_de, 'T_vs_C', lfc_threshold=0.1)
        titles = [r.title for r in results]
        assert any('Modest' in t for t in titles)


# =============================================================================
# CellTypeDeconvolution
# =============================================================================

class TestCellTypeDeconvolution:
    """Test cell type deconvolution module."""

    def test_default_skin_signature(self):
        from cell_deconvolution import CellTypeDeconvolution
        deconv = CellTypeDeconvolution()
        sig = deconv.get_default_skin_signature()
        assert sig.shape[0] >= 20
        assert sig.shape[1] >= 4
        assert 'Keratinocyte' in sig.columns
        assert 'Fibroblast' in sig.columns

    def test_nnls_basic(self):
        from cell_deconvolution import CellTypeDeconvolution, DeconvolutionResult
        deconv = CellTypeDeconvolution()
        sig = deconv.get_default_skin_signature()
        
        mock_expr = pd.DataFrame(
            np.random.rand(len(sig), 3) * 100,
            index=sig.index,
            columns=['S1', 'S2', 'S3']
        )
        result = deconv.run_nnls(mock_expr, sig)
        assert isinstance(result, DeconvolutionResult)
        assert result.proportions.shape == (3, sig.shape[1])
        assert result.method == 'NNLS'

    def test_nnls_proportions_sum_to_one(self):
        from cell_deconvolution import CellTypeDeconvolution
        deconv = CellTypeDeconvolution()
        sig = deconv.get_default_skin_signature()
        
        mock_expr = pd.DataFrame(
            np.random.rand(len(sig), 5) * 100,
            index=sig.index,
            columns=[f'S{i}' for i in range(5)]
        )
        result = deconv.run_nnls(mock_expr, sig)
        row_sums = result.proportions.sum(axis=1)
        assert np.allclose(row_sums, 1.0, atol=0.01)

    def test_nnls_insufficient_overlap(self):
        from cell_deconvolution import CellTypeDeconvolution
        deconv = CellTypeDeconvolution()
        sig = deconv.get_default_skin_signature()
        
        # Expression with completely different genes
        mock_expr = pd.DataFrame(
            np.random.rand(10, 3) * 100,
            index=[f'FAKE_GENE_{i}' for i in range(10)],
            columns=['S1', 'S2', 'S3']
        )
        with pytest.raises(ValueError, match="Insufficient gene overlap"):
            deconv.run_nnls(mock_expr, sig)


# =============================================================================
# AdvancedQC
# =============================================================================

class TestAdvancedQC:
    """Test advanced QC module."""

    def _make_pca_df(self):
        return pd.DataFrame({
            'PC1': [1.0, 1.2, 0.8, 1.1, 10.0],
            'PC2': [0.5, 0.6, 0.4, 0.5, 10.0]
        }, index=['S1', 'S2', 'S3', 'S4', 'Outlier'])

    def test_detect_outliers_basic(self):
        from advanced_qc import AdvancedQC
        qc = AdvancedQC()
        pca_df = self._make_pca_df()
        outliers, distances = qc.detect_outliers(pca_df, threshold_sd=2.0)
        assert 'Outlier' in outliers
        assert len(distances) == 5

    def test_detect_outliers_none(self):
        from advanced_qc import AdvancedQC
        qc = AdvancedQC()
        pca_df = pd.DataFrame({
            'PC1': [1.0, 1.1, 0.9, 1.0],
            'PC2': [0.5, 0.5, 0.5, 0.5]
        }, index=['S1', 'S2', 'S3', 'S4'])
        outliers, _ = qc.detect_outliers(pca_df)
        assert len(outliers) == 0

    def test_batch_effects(self):
        from advanced_qc import AdvancedQC
        qc = AdvancedQC()
        pca_df = self._make_pca_df()
        batch = pd.Series(['A', 'A', 'B', 'B', 'A'], index=pca_df.index)
        results = qc.assess_batch_effects(pca_df, batch)
        assert 'PC1' in results
        assert 'PC2' in results
        assert 0 <= results['PC1'] <= 1
        assert 0 <= results['PC2'] <= 1

    def test_outlier_plot_returns_plotly(self):
        from advanced_qc import AdvancedQC
        qc = AdvancedQC()
        pca_df = self._make_pca_df()
        outliers, distances = qc.detect_outliers(pca_df)
        conds = {s: 'ctrl' for s in pca_df.index}
        fig = qc.create_outlier_plot(pca_df, conds, outliers, distances)
        assert isinstance(fig, go.Figure)

    def test_batch_plot_returns_plotly(self):
        from advanced_qc import AdvancedQC
        qc = AdvancedQC()
        results = {'PC1': 0.25, 'PC2': 0.10}
        fig = qc.create_batch_effect_plot(results)
        assert isinstance(fig, go.Figure)


# =============================================================================
# Visualization Extensions
# =============================================================================

class TestVisualizationExtensions:
    """Test new visualization functions."""

    def test_enrichment_dotplot(self):
        from visualizations import create_enrichment_dotplot
        mock_enr = pd.DataFrame({
            'Term': ['Pathway A', 'Pathway B', 'Pathway C'],
            'Overlap': ['5/100', '3/80', '10/200'],
            'P-value': [0.001, 0.01, 0.05],
            'Adjusted P-value': [0.01, 0.05, 0.1],
            'Genes': ['A;B;C;D;E', 'F;G;H', 'I;J;K']
        })
        fig = create_enrichment_dotplot(mock_enr)
        assert isinstance(fig, go.Figure)

    def test_enrichment_dotplot_empty(self):
        from visualizations import create_enrichment_dotplot
        fig = create_enrichment_dotplot(pd.DataFrame())
        assert isinstance(fig, go.Figure)

    def test_deconvolution_plot(self):
        from visualizations import create_deconvolution_plot
        mock_props = pd.DataFrame({
            'Keratinocyte': [0.5, 0.6, 0.4],
            'Fibroblast': [0.3, 0.2, 0.35],
            'Immune': [0.2, 0.2, 0.25]
        }, index=['S1', 'S2', 'S3'])
        fig = create_deconvolution_plot(mock_props)
        assert isinstance(fig, go.Figure)

    def test_deconvolution_plot_with_conditions(self):
        from visualizations import create_deconvolution_plot
        mock_props = pd.DataFrame({
            'CellA': [0.6, 0.4],
            'CellB': [0.4, 0.6]
        }, index=['S1', 'S2'])
        conds = {'S1': 'Ctrl', 'S2': 'Treat'}
        fig = create_deconvolution_plot(mock_props, conds)
        assert isinstance(fig, go.Figure)


# =============================================================================
# EnrichmentResult Extension
# =============================================================================

class TestEnrichmentResultExtension:
    """Test backward-compatible EnrichmentResult extension."""

    def test_backward_compatible(self):
        from export_engine import EnrichmentResult
        old = EnrichmentResult(
            go_results=pd.DataFrame(),
            kegg_results=pd.DataFrame(),
            genes_used=['A', 'B'],
            selection_note='test'
        )
        assert old.error is None
        assert old.go_results.empty
        assert old.kegg_results.empty

    def test_new_fields_defaults(self):
        from export_engine import EnrichmentResult
        obj = EnrichmentResult(
            go_results=pd.DataFrame(),
            kegg_results=pd.DataFrame(),
            genes_used=['A'],
            selection_note='test'
        )
        assert isinstance(obj.gsea_results, dict)
        assert obj.reactome_results.empty
        assert obj.hallmark_results.empty
        assert obj.go_mf_results.empty
        assert obj.go_cc_results.empty
        assert obj.wikipathway_results.empty

    def test_new_fields_populated(self):
        from export_engine import EnrichmentResult
        obj = EnrichmentResult(
            go_results=pd.DataFrame({'Term': ['A']}),
            kegg_results=pd.DataFrame(),
            genes_used=['X'],
            selection_note='test',
            reactome_results=pd.DataFrame({'Term': ['B']}),
            gsea_results={'Hallmark': pd.DataFrame({'Term': ['C']})}
        )
        assert len(obj.reactome_results) == 1
        assert 'Hallmark' in obj.gsea_results
