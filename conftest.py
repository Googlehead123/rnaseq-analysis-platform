"""
Pytest configuration and fixtures for RNA-seq Analysis Platform tests.
"""

import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, Mock
import pytest
import pandas as pd
import numpy as np


# ============================================================================
# Streamlit Mocking Fixtures
# ============================================================================


@pytest.fixture
def mock_streamlit(monkeypatch):
    """Mock Streamlit module to prevent actual UI rendering during tests."""
    mock_st = MagicMock()

    # Mock common Streamlit functions
    mock_st.title = MagicMock()
    mock_st.header = MagicMock()
    mock_st.subheader = MagicMock()
    mock_st.write = MagicMock()
    mock_st.markdown = MagicMock()
    mock_st.error = MagicMock()
    mock_st.warning = MagicMock()
    mock_st.success = MagicMock()
    mock_st.info = MagicMock()
    mock_st.text = MagicMock()
    mock_st.code = MagicMock()
    mock_st.dataframe = MagicMock()
    mock_st.table = MagicMock()
    mock_st.metric = MagicMock()
    mock_st.json = MagicMock()

    # Mock input widgets
    mock_st.button = MagicMock(return_value=False)
    mock_st.checkbox = MagicMock(return_value=False)
    mock_st.radio = MagicMock(return_value=None)
    mock_st.selectbox = MagicMock(return_value=None)
    mock_st.multiselect = MagicMock(return_value=[])
    mock_st.slider = MagicMock(return_value=0)
    mock_st.select_slider = MagicMock(return_value=0)
    mock_st.text_input = MagicMock(return_value="")
    mock_st.number_input = MagicMock(return_value=0)
    mock_st.text_area = MagicMock(return_value="")
    mock_st.date_input = MagicMock()
    mock_st.time_input = MagicMock()
    mock_st.file_uploader = MagicMock(return_value=None)
    mock_st.color_picker = MagicMock(return_value="#000000")

    # Mock layout
    mock_st.sidebar = MagicMock()
    mock_st.columns = MagicMock(return_value=[MagicMock(), MagicMock()])
    mock_st.container = MagicMock()
    mock_st.expander = MagicMock()
    mock_st.tabs = MagicMock(return_value=[MagicMock()])

    # Mock media
    mock_st.image = MagicMock()
    mock_st.audio = MagicMock()
    mock_st.video = MagicMock()

    # Mock charts
    mock_st.plotly_chart = MagicMock()
    mock_st.pyplot = MagicMock()
    mock_st.altair_chart = MagicMock()
    mock_st.vega_lite_chart = MagicMock()
    mock_st.line_chart = MagicMock()
    mock_st.area_chart = MagicMock()
    mock_st.bar_chart = MagicMock()
    mock_st.map = MagicMock()

    # Mock status
    mock_st.progress = MagicMock()
    mock_st.spinner = MagicMock()
    mock_st.balloons = MagicMock()
    mock_st.snow = MagicMock()

    # Mock control flow
    mock_st.stop = MagicMock()
    mock_st.form = MagicMock()
    mock_st.form_submit_button = MagicMock(return_value=False)

    # Mock session state
    mock_st.session_state = {}

    # Mock set_page_config
    mock_st.set_page_config = MagicMock()

    # Mock download button
    mock_st.download_button = MagicMock(return_value=False)

    monkeypatch.setattr("streamlit", mock_st)
    return mock_st


@pytest.fixture
def mock_session_state(mock_streamlit):
    """Provide a clean session state dictionary for each test."""
    session_state = {}
    mock_streamlit.session_state = session_state
    return session_state


# ============================================================================
# Test Data Directory Fixtures
# ============================================================================


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return Path(__file__).parent / "tests" / "data"


# ============================================================================
# Sample Data Fixtures
# ============================================================================


@pytest.fixture
def sample_counts_df():
    """
    Sample RNA-seq count matrix for testing.
    Shape: (10 samples, 100 genes)
    """
    np.random.seed(42)
    data = np.random.negative_binomial(n=10, p=0.1, size=(10, 100))
    samples = [f"sample_{i + 1}" for i in range(10)]
    genes = [f"gene_{i + 1}" for i in range(100)]
    df = pd.DataFrame(data, index=samples, columns=genes)
    df.index.name = "sample_id"
    return df


@pytest.fixture
def sample_metadata_df():
    """
    Sample metadata for testing.
    Shape: (10 samples, 2 columns: sample_id, condition)
    """
    samples = [f"sample_{i + 1}" for i in range(10)]
    conditions = ["control"] * 5 + ["treatment"] * 5
    df = pd.DataFrame({"sample_id": samples, "condition": conditions})
    return df


@pytest.fixture
def sample_conditions_dict():
    """Sample conditions dictionary for testing."""
    return {f"sample_{i + 1}": "control" if i < 5 else "treatment" for i in range(10)}


@pytest.fixture
def sample_de_results_df():
    """
    Sample differential expression results for testing.
    Contains typical DESeq2 output columns.
    """
    np.random.seed(42)
    n_genes = 100
    genes = [f"gene_{i + 1}" for i in range(n_genes)]

    df = pd.DataFrame(
        {
            "gene_id": genes,
            "baseMean": np.random.uniform(10, 1000, n_genes),
            "log2FoldChange": np.random.normal(0, 2, n_genes),
            "lfcSE": np.random.uniform(0.1, 0.5, n_genes),
            "stat": np.random.normal(0, 3, n_genes),
            "pvalue": np.random.uniform(0, 1, n_genes),
            "padj": np.random.uniform(0, 1, n_genes),
        }
    )

    # Ensure some significant genes
    df.loc[:10, "padj"] = np.random.uniform(0, 0.05, 11)
    df.loc[:10, "log2FoldChange"] = np.random.uniform(1.5, 3, 11)

    return df


@pytest.fixture
def sample_log_normalized_df():
    """
    Sample log-normalized expression data for testing.
    Shape: (10 samples, 100 genes)
    """
    np.random.seed(42)
    data = np.random.uniform(0, 100, size=(10, 100))
    samples = [f"sample_{i + 1}" for i in range(10)]
    genes = [f"gene_{i + 1}" for i in range(100)]
    df = pd.DataFrame(data, index=samples, columns=genes)
    df.index.name = "sample_id"
    return df


# ============================================================================
# External API Mocking Fixtures
# ============================================================================


@pytest.fixture
def mock_enrichr_response():
    """Mock Enrichr API response for pathway enrichment testing."""
    return {
        "GO_Biological_Process_2021": [
            ["immune response", 0.001, 0.01, 0.05, ["gene_1", "gene_2", "gene_3"]],
            ["cell cycle", 0.005, 0.02, 0.08, ["gene_4", "gene_5"]],
            ["apoptosis", 0.01, 0.03, 0.10, ["gene_6", "gene_7", "gene_8"]],
        ]
    }


@pytest.fixture
def mock_gseapy(monkeypatch):
    """Mock gseapy module for enrichment analysis testing."""
    mock_gp = MagicMock()

    # Mock enrichr function
    mock_enrichr_result = MagicMock()
    mock_enrichr_result.results = pd.DataFrame(
        {
            "Term": ["immune response", "cell cycle", "apoptosis"],
            "P-value": [0.001, 0.005, 0.01],
            "Adjusted P-value": [0.01, 0.02, 0.03],
            "Odds Ratio": [2.5, 2.0, 1.8],
            "Combined Score": [50, 40, 35],
            "Genes": ["gene_1;gene_2;gene_3", "gene_4;gene_5", "gene_6;gene_7;gene_8"],
        }
    )
    mock_gp.enrichr = MagicMock(return_value=mock_enrichr_result)

    # Mock GSEA function
    mock_gsea_result = MagicMock()
    mock_gsea_result.res2d = pd.DataFrame(
        {
            "Term": ["pathway_1", "pathway_2"],
            "ES": [0.5, -0.4],
            "NES": [2.0, -1.8],
            "NOM p-val": [0.001, 0.005],
            "FDR q-val": [0.01, 0.02],
            "FWER p-val": [0.01, 0.03],
            "Lead_genes": ["gene_1;gene_2", "gene_3;gene_4"],
        }
    )
    mock_gp.prerank = MagicMock(return_value=mock_gsea_result)

    monkeypatch.setattr("pathway_enrichment.gp", mock_gp)
    return mock_gp


@pytest.fixture
def mock_de_result_fixture(monkeypatch):
    """
    Mock PyDESeq2 DESeqDataSet and DeseqStats for differential expression testing.
    This fixture imports DEResult at runtime to avoid import errors during collection.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats

        # Create mock DESeqDataSet
        mock_dds = MagicMock(spec=DeseqDataSet)
        mock_dds.varm = {"log2FoldChange": np.random.normal(0, 2, 100)}
        mock_dds.var_names = [f"gene_{i + 1}" for i in range(100)]

        # Create mock DeseqStats
        mock_ds = MagicMock(spec=DeseqStats)
        mock_ds.results_df = pd.DataFrame(
            {
                "baseMean": np.random.uniform(10, 1000, 100),
                "log2FoldChange": np.random.normal(0, 2, 100),
                "lfcSE": np.random.uniform(0.1, 0.5, 100),
                "stat": np.random.normal(0, 3, 100),
                "pvalue": np.random.uniform(0, 1, 100),
                "padj": np.random.uniform(0, 1, 100),
            }
        )

        return {"dds": mock_dds, "ds": mock_ds}
    except ImportError:
        # If PyDESeq2 not installed, return None
        return None


# ============================================================================
# Evidence Directory Fixture
# ============================================================================


@pytest.fixture(autouse=True)
def ensure_evidence_dir():
    """Ensure .sisyphus/evidence/ directory exists for test artifacts."""
    evidence_dir = Path(__file__).parent / ".sisyphus" / "evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    return evidence_dir
