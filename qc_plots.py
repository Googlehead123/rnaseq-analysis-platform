"""Pre-analysis QC visualizations for RNA-seq data."""

import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px


def create_library_size_barplot(counts_df: pd.DataFrame) -> go.Figure:
    """
    Bar plot of total counts (library size) per sample, sorted descending.

    Args:
        counts_df: samples × genes DataFrame of raw counts

    Returns:
        Plotly Figure object
    """
    lib_sizes = counts_df.sum(axis=1).sort_values(ascending=False)
    mean_size = lib_sizes.mean()

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=lib_sizes.index.tolist(),
            y=lib_sizes.values,
            marker_color="steelblue",
            name="Library Size",
        )
    )
    fig.add_hline(
        y=mean_size,
        line_dash="dash",
        line_color="red",
        annotation_text=f"Mean: {mean_size:,.0f}",
        annotation_position="top right",
    )
    fig.update_layout(
        title="Library Size per Sample",
        xaxis_title="Sample",
        yaxis_title="Total Counts",
        showlegend=False,
    )
    return fig


def create_count_distribution_boxplot(
    counts_df: pd.DataFrame, log_transform: bool = True
) -> go.Figure:
    """
    Box plot of expression distribution per sample.

    Args:
        counts_df: samples × genes DataFrame
        log_transform: Apply log2(x+1) transformation (default: True)

    Returns:
        Plotly Figure object
    """
    data = counts_df.copy()
    if log_transform:
        data = np.log2(data + 1)

    fig = go.Figure()
    for sample in data.index:
        fig.add_trace(
            go.Box(
                y=data.loc[sample].values,
                name=str(sample),
                showlegend=False,
            )
        )

    ylabel = "log₂(count + 1)" if log_transform else "Count"
    fig.update_layout(
        title="Expression Distribution per Sample",
        xaxis_title="Sample",
        yaxis_title=ylabel,
    )
    return fig


def create_gene_detection_plot(
    counts_df: pd.DataFrame, threshold: int = 0
) -> go.Figure:
    """
    Bar plot of number of detected genes per sample.

    Args:
        counts_df: samples × genes DataFrame
        threshold: Minimum count to consider a gene detected (default: 0)

    Returns:
        Plotly Figure object
    """
    detected = (counts_df > threshold).sum(axis=1).sort_values(ascending=False)

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=detected.index.tolist(),
            y=detected.values,
            marker_color="darkorange",
            name="Detected Genes",
        )
    )
    fig.update_layout(
        title=f"Genes Detected per Sample (count > {threshold})",
        xaxis_title="Sample",
        yaxis_title="Number of Genes",
        showlegend=False,
    )
    return fig


def create_sample_similarity_heatmap(
    counts_df: pd.DataFrame, method: str = "spearman"
) -> go.Figure:
    """
    Pairwise sample correlation heatmap.

    Args:
        counts_df: samples × genes DataFrame
        method: Correlation method ('spearman' or 'pearson')

    Returns:
        Plotly Figure object
    """
    corr_matrix = counts_df.T.corr(method=method)

    text_vals = [[f"{v:.3f}" for v in row] for row in corr_matrix.values]

    fig = go.Figure(
        data=go.Heatmap(
            z=corr_matrix.values,
            x=corr_matrix.columns.tolist(),
            y=corr_matrix.index.tolist(),
            colorscale="RdBu_r",
            zmid=0,
            text=text_vals,
            texttemplate="%{text}",
            hovertemplate="Sample X: %{x}<br>Sample Y: %{y}<br>Correlation: %{z:.3f}<extra></extra>",
        )
    )
    fig.update_layout(
        title=f"Sample Similarity ({method.capitalize()} Correlation)",
        width=600,
        height=600,
    )
    return fig
