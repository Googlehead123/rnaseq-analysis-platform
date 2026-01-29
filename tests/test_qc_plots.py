"""Tests for QC plot functions."""
import pytest
import pandas as pd
import numpy as np
from qc_plots import (
    create_library_size_barplot,
    create_count_distribution_boxplot,
    create_gene_detection_plot,
    create_sample_similarity_heatmap,
)


@pytest.fixture
def sample_counts():
    np.random.seed(42)
    return pd.DataFrame(
        np.random.poisson(100, (6, 200)),
        index=[f"Sample_{i}" for i in range(6)],
        columns=[f"Gene_{i}" for i in range(200)],
    )


def test_library_size_barplot(sample_counts):
    fig = create_library_size_barplot(sample_counts)
    assert fig is not None
    assert hasattr(fig, "data")


def test_count_distribution_boxplot(sample_counts):
    fig = create_count_distribution_boxplot(sample_counts)
    assert fig is not None
    assert hasattr(fig, "data")


def test_count_distribution_boxplot_no_log(sample_counts):
    fig = create_count_distribution_boxplot(sample_counts, log_transform=False)
    assert fig is not None


def test_gene_detection_plot(sample_counts):
    fig = create_gene_detection_plot(sample_counts)
    assert fig is not None
    assert hasattr(fig, "data")


def test_sample_similarity_heatmap(sample_counts):
    fig = create_sample_similarity_heatmap(sample_counts)
    assert fig is not None
    assert hasattr(fig, "data")
