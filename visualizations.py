"""
Interactive visualizations for RNA-seq analysis using Plotly.

Provides volcano plots, clustered heatmaps, and PCA plots.
"""

from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist
from sklearn.decomposition import PCA
from de_analysis import ensure_gene_column


def create_volcano_plot(
    results_df: pd.DataFrame, lfc_threshold: float = 1.0, padj_threshold: float = 0.05,
    top_n_labels: int = 10,
) -> go.Figure:
    """
    Create interactive volcano plot from DE results.

    Args:
        results_df: DataFrame with columns: gene, log2FoldChange, padj
        lfc_threshold: Log2 fold change threshold for significance (default: 1.0)
        padj_threshold: Adjusted p-value threshold (default: 0.05)

    Returns:
        Plotly Figure object
    """
    # Validate input
    if results_df is None or results_df.empty:
        raise ValueError(
            "Cannot create volcano plot: results_df is empty or None. "
            "Ensure your differential expression analysis produced results. "
            "Check that your input data has samples and genes."
        )

    # Ensure gene column is properly named (defensive)
    results_df = ensure_gene_column(results_df)

    required_cols = ["gene", "log2FoldChange", "padj"]
    missing = [col for col in required_cols if col not in results_df.columns]
    if missing:
        available = ", ".join(results_df.columns.tolist()[:5])
        extra = (
            f"... ({len(results_df.columns) - 5} more)"
            if len(results_df.columns) > 5
            else ""
        )
        raise ValueError(
            f"Cannot create volcano plot: missing required columns {missing}. "
            f"Found columns: {available}{extra}. "
            f"Suggestion: Ensure your DE results contain 'gene', 'log2FoldChange' (or 'log2FC'), "
            f"and 'padj' (or 'FDR', 'adjusted_pvalue'). Check column names for typos."
        )

    # Add -log10(padj) column
    df = results_df.copy()

    if df["padj"].isna().any():
        n_nan = df["padj"].isna().sum()
        raise ValueError(
            f"Cannot create volcano plot: {n_nan} NaN values in padj column. "
            f"Ensure differential expression analysis completed successfully."
        )

    df["-log10_padj"] = -np.log10(df["padj"].clip(lower=1e-300))  # Clip to avoid inf

    # Classify significance
    def classify(row):
        if row["padj"] >= padj_threshold:
            return "NS"  # Not significant
        elif row["log2FoldChange"] > lfc_threshold:
            return "Up"
        elif row["log2FoldChange"] < -lfc_threshold:
            return "Down"
        else:
            return "NS"

    df["significance"] = df.apply(classify, axis=1)

    fig = px.scatter(
        df,
        x="log2FoldChange",
        y="-log10_padj",
        color="significance",
        hover_name="gene",
        hover_data={
            "log2FoldChange": ":.2f",
            "padj": ":.2e",
            "-log10_padj": False,
            "significance": False,
        },
        color_discrete_map={"Up": "red", "Down": "blue", "NS": "lightgray"},
        labels={"log2FoldChange": "log₂(Fold Change)", "-log10_padj": "-log₁₀(padj)"},
    )

    # Add threshold lines
    fig.add_hline(y=-np.log10(padj_threshold), line_dash="dash", line_color="gray")
    fig.add_vline(x=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_vline(x=-lfc_threshold, line_dash="dash", line_color="gray")

    # Add text annotations for top N significant genes
    if top_n_labels > 0:
        top_genes = (
            df[df["padj"] < padj_threshold]
            .nsmallest(top_n_labels, "padj")
        )
        if not top_genes.empty:
            fig.add_trace(
                go.Scatter(
                    x=top_genes["log2FoldChange"],
                    y=top_genes["-log10_padj"],
                    mode="text",
                    text=top_genes["gene"],
                    textposition="top center",
                    textfont=dict(size=9),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    fig.update_layout(title="Volcano Plot", showlegend=True)

    return fig


def create_clustered_heatmap(
    expression_df: pd.DataFrame,
    sample_conditions: Dict[str, str],
    de_results_df: Optional[pd.DataFrame] = None,
    top_n_genes: int = 50,
    z_score: bool = True,
    gene_selection: str = "de",
) -> go.Figure:
    """
    Create heatmap with row (gene) clustering dendrogram.

    Gene selection logic:
    - If gene_selection="de" AND de_results_df provided: use top N by padj (most significant)
    - If gene_selection="variance" OR no de_results_df: use top N by variance (fallback)

    Args:
        expression_df: genes × samples (log2 normalized)
        sample_conditions: Dict[sample_name, condition]
        de_results_df: Optional DE results with 'gene', 'padj' columns
        top_n_genes: Number of genes to display (default: 50)
        z_score: Apply z-score normalization (default: True)
        gene_selection: "de" (default) or "variance"

    Returns:
        Plotly Figure object

    Note: We only cluster ROWS (genes), not columns (samples).
    Samples are grouped by condition for clearer visualization.
    """
    # Validate input
    if expression_df is None or expression_df.empty:
        raise ValueError(
            "Cannot create heatmap: expression_df is empty or None. "
            "Ensure your expression data contains samples and genes."
        )

    if not sample_conditions:
        raise ValueError(
            "Cannot create heatmap: sample_conditions is empty. "
            "Provide a mapping of sample names to experimental conditions."
        )

    missing_samples = [s for s in expression_df.columns if s not in sample_conditions]
    if missing_samples:
        provided = ", ".join(list(sample_conditions.keys())[:3])
        extra = (
            f"... ({len(sample_conditions) - 3} more)"
            if len(sample_conditions) > 3
            else ""
        )
        raise ValueError(
            f"Cannot create heatmap: {len(missing_samples)} samples missing from sample_conditions. "
            f"Missing: {', '.join(missing_samples[:3])}{'...' if len(missing_samples) > 3 else ''}. "
            f"Provided conditions for: {provided}{extra}. "
            f"Suggestion: Ensure all sample names in expression_df are mapped in sample_conditions."
        )

    # Gene selection: prioritize DE genes, fallback to variance
    if gene_selection == "de" and de_results_df is not None:
        # Top N genes by significance (lowest padj)
        sig_genes = (
            de_results_df.dropna(subset=["padj"])
            .nsmallest(top_n_genes, "padj")["gene"]
            .tolist()
        )
        # Filter to genes present in expression matrix
        top_genes = [g for g in sig_genes if g in expression_df.index]
        if len(top_genes) < 10:  # Fallback if too few significant genes
            gene_vars = expression_df.var(axis=1)
            top_genes = gene_vars.nlargest(top_n_genes).index.tolist()
    else:
        # Fallback: top N genes by variance
        gene_vars = expression_df.var(axis=1)
        top_genes = gene_vars.nlargest(top_n_genes).index.tolist()

    plot_data = expression_df.loc[top_genes]

    # Z-score normalize if requested (per gene, across samples)
    if z_score:
        plot_data = plot_data.sub(plot_data.mean(axis=1), axis=0).div(
            plot_data.std(axis=1), axis=0
        )

    # Sort samples by condition (group conditions together)
    sample_order = sorted(plot_data.columns, key=lambda s: sample_conditions.get(s, ""))
    plot_data = plot_data[sample_order]

    # Hierarchical clustering on genes (rows)
    if len(plot_data) > 1:
        linkage_matrix = linkage(
            pdist(plot_data.values, metric="euclidean"), method="average"
        )
        gene_order = leaves_list(linkage_matrix)
        plot_data = plot_data.iloc[gene_order]

    # Create heatmap
    fig = go.Figure(
        data=go.Heatmap(
            z=plot_data.values,
            x=plot_data.columns,
            y=plot_data.index,
            colorscale="RdBu_r",
            zmid=0,
            hovertemplate="Gene: %{y}<br>Sample: %{x}<br>Expression: %{z:.2f}<extra></extra>",
        )
    )

    fig.update_layout(
        title=f"Clustered Heatmap (Top {len(plot_data)} Genes)",
        xaxis_title="Samples",
        yaxis_title="Genes",
        height=max(400, len(plot_data) * 10),  # Scale height with gene count
    )

    return fig


def create_pca_plot(
    expression_df: pd.DataFrame, sample_conditions: Dict[str, str],
    show_ellipses: bool = True,
) -> go.Figure:
    """
    Create PCA plot for sample clustering visualization.

    Args:
        expression_df: samples × genes (log2 normalized)
        sample_conditions: Dict[sample_name, condition]

    Returns:
        Plotly Figure object
    """
    # Validate input
    if expression_df is None or expression_df.empty:
        raise ValueError(
            "Cannot create PCA plot: expression_df is empty or None. "
            "Ensure your expression data contains samples and genes."
        )

    if expression_df.shape[1] < 2:
        raise ValueError(
            f"Cannot create PCA plot: requires at least 2 samples, but got {expression_df.shape[1]}. "
            f"Suggestion: PCA needs multiple samples to compute principal components. "
            f"Ensure your expression data has at least 2 samples."
        )

    if not sample_conditions:
        raise ValueError(
            "Cannot create PCA plot: sample_conditions is empty. "
            "Provide a mapping of sample names to experimental conditions for coloring."
        )

    # Run PCA (samples × genes → PC space)
    pca = PCA(n_components=min(3, expression_df.shape[0]))
    pca_result = pca.fit_transform(expression_df.values)

    # Create DataFrame for plotting
    pca_df = pd.DataFrame(
        pca_result[:, :2],  # First 2 PCs
        columns=["PC1", "PC2"],
        index=expression_df.index,
    )
    pca_df["condition"] = [sample_conditions.get(s, "Unknown") for s in pca_df.index]
    pca_df["sample"] = pca_df.index

    # Create scatter plot
    fig = px.scatter(
        pca_df,
        x="PC1",
        y="PC2",
        color="condition",
        hover_name="sample",
        labels={
            "PC1": f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)",
            "PC2": f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)",
        },
    )

    # Add 95% confidence ellipses per condition group
    if show_ellipses:
        colors = px.colors.qualitative.Plotly
        conditions = sorted(pca_df["condition"].unique())
        for i, cond in enumerate(conditions):
            group = pca_df[pca_df["condition"] == cond]
            if len(group) < 3:
                continue
            mean_x, mean_y = group["PC1"].mean(), group["PC2"].mean()
            cov = np.cov(group["PC1"].values, group["PC2"].values)
            eigenvalues, eigenvectors = np.linalg.eigh(cov)
            # Sort descending
            order = eigenvalues.argsort()[::-1]
            eigenvalues = eigenvalues[order]
            eigenvectors = eigenvectors[:, order]
            # 95% confidence: chi2(df=2) = 5.991
            scale = np.sqrt(5.991)
            theta = np.linspace(0, 2 * np.pi, 100)
            ellipse = np.array([np.cos(theta), np.sin(theta)])
            # Transform
            transform = eigenvectors @ np.diag(np.sqrt(np.maximum(eigenvalues, 0)) * scale)
            ellipse_pts = (transform @ ellipse).T + np.array([mean_x, mean_y])
            fig.add_trace(
                go.Scatter(
                    x=ellipse_pts[:, 0],
                    y=ellipse_pts[:, 1],
                    mode="lines",
                    line=dict(color=colors[i % len(colors)], dash="dash", width=1.5),
                    showlegend=False,
                    hoverinfo="skip",
                )
            )

    fig.update_layout(title="PCA Plot", showlegend=True)

    return fig


def create_ma_plot(
    results_df: pd.DataFrame, padj_threshold: float = 0.05, lfc_threshold: float = 1.0
) -> go.Figure:
    """
    Create MA plot (log mean expression vs log2 fold change).

    Args:
        results_df: DataFrame with columns: gene, log2FoldChange, padj, baseMean
        padj_threshold: Adjusted p-value threshold (default: 0.05)
        lfc_threshold: Log2 fold change threshold (default: 1.0)

    Returns:
        Plotly Figure object
    """
    if results_df is None or results_df.empty:
        raise ValueError("Cannot create MA plot: results_df is empty or None.")

    results_df = ensure_gene_column(results_df)

    required_cols = ["gene", "log2FoldChange", "padj", "baseMean"]
    missing = [c for c in required_cols if c not in results_df.columns]
    if missing:
        raise ValueError(f"Cannot create MA plot: missing columns {missing}.")

    df = results_df.copy()
    df["log10_baseMean"] = np.log10(df["baseMean"] + 1)

    def classify(row):
        if row["padj"] >= padj_threshold:
            return "NS"
        elif row["log2FoldChange"] > lfc_threshold:
            return "Up"
        elif row["log2FoldChange"] < -lfc_threshold:
            return "Down"
        else:
            return "NS"

    df["significance"] = df.apply(classify, axis=1)

    fig = px.scatter(
        df,
        x="log10_baseMean",
        y="log2FoldChange",
        color="significance",
        hover_name="gene",
        hover_data={
            "log2FoldChange": ":.2f",
            "padj": ":.2e",
            "log10_baseMean": False,
            "significance": False,
        },
        color_discrete_map={"Up": "red", "Down": "blue", "NS": "lightgray"},
        labels={
            "log10_baseMean": "log₁₀(baseMean + 1)",
            "log2FoldChange": "log₂(Fold Change)",
        },
    )

    fig.add_hline(y=lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_hline(y=-lfc_threshold, line_dash="dash", line_color="gray")
    fig.add_hline(y=0, line_color="black", line_width=0.5)

    fig.update_layout(title="MA Plot", showlegend=True)

    return fig


def create_correlation_heatmap(
    expression_df: pd.DataFrame, sample_conditions: dict,
    method: str = "spearman",
) -> go.Figure:
    """
    Create sample-to-sample correlation heatmap.

    Args:
        expression_df: samples × genes DataFrame
        sample_conditions: Dict mapping sample names to conditions
        method: Correlation method ('spearman' or 'pearson')

    Returns:
        Plotly Figure object
    """
    if expression_df is None or expression_df.empty:
        raise ValueError("Cannot create correlation heatmap: expression_df is empty or None.")

    # Compute sample-to-sample correlations
    corr_matrix = expression_df.T.corr(method=method)

    # Sort samples by condition
    sorted_samples = sorted(
        corr_matrix.columns,
        key=lambda s: sample_conditions.get(s, ""),
    )
    corr_matrix = corr_matrix.loc[sorted_samples, sorted_samples]

    # Format text annotations
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
        title=f"Sample Correlation Heatmap ({method.capitalize()})",
        width=600,
        height=600,
    )

    return fig


def compute_de_summary(
    results_df: pd.DataFrame, padj_threshold: float = 0.05, lfc_threshold: float = 1.0
) -> dict:
    """
    Compute summary statistics from DE results.

    Args:
        results_df: DataFrame with columns: gene, log2FoldChange, padj
        padj_threshold: Adjusted p-value threshold (default: 0.05)
        lfc_threshold: Log2 fold change threshold (default: 1.0)

    Returns:
        Dict with total_genes, significant_genes, upregulated, downregulated,
        top_up_genes, top_down_genes, top_significant
    """
    if results_df is None or results_df.empty:
        raise ValueError("Cannot compute DE summary: results_df is empty or None.")

    results_df = ensure_gene_column(results_df)
    df = results_df.copy()

    sig = df[df["padj"] < padj_threshold]
    up = sig[sig["log2FoldChange"] > lfc_threshold]
    down = sig[sig["log2FoldChange"] < -lfc_threshold]

    top_up = (
        up.nlargest(10, "log2FoldChange")[["gene", "log2FoldChange", "padj"]]
        .apply(lambda r: (r["gene"], r["log2FoldChange"], r["padj"]), axis=1)
        .tolist()
    )
    top_down = (
        down.nsmallest(10, "log2FoldChange")[["gene", "log2FoldChange", "padj"]]
        .apply(lambda r: (r["gene"], r["log2FoldChange"], r["padj"]), axis=1)
        .tolist()
    )
    top_significant = (
        sig.nsmallest(10, "padj")[["gene", "log2FoldChange", "padj"]]
        .apply(lambda r: (r["gene"], r["log2FoldChange"], r["padj"]), axis=1)
        .tolist()
    )

    return {
        "total_genes": len(df),
        "significant_genes": len(sig),
        "upregulated": len(up),
        "downregulated": len(down),
        "top_up_genes": top_up,
        "top_down_genes": top_down,
        "top_significant": top_significant,
    }
