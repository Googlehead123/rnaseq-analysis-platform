"""
Interactive visualizations for RNA-seq analysis using Plotly.

Provides volcano plots, clustered heatmaps, and PCA plots.
"""

from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
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

    # Drop rows with NaN padj (normal for low-count genes filtered by PyDESeq2)
    n_nan = df["padj"].isna().sum()
    if n_nan > 0:
        df = df.dropna(subset=["padj"])
    if df.empty:
        raise ValueError(
            "Cannot create volcano plot: all padj values are NaN. "
            "Ensure differential expression analysis completed successfully."
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
    # Drop rows with NaN padj or log2FoldChange (normal for filtered genes)
    df = df.dropna(subset=["padj", "log2FoldChange", "baseMean"])
    if df.empty:
        raise ValueError("Cannot create MA plot: no valid data after removing NaN values.")

    df["log10_baseMean"] = np.log10(df["baseMean"] + 1)

    def classify(row):
        if pd.isna(row["padj"]) or row["padj"] >= padj_threshold:
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
    df = results_df.dropna(subset=["padj", "log2FoldChange"]).copy()

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


def create_enrichment_dotplot(
    enrichment_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "Enrichment Results",
) -> go.Figure:
    """
    Create enrichment dot plot using Plotly.

    Args:
        enrichment_df: DataFrame with columns: Term, Overlap, P-value, Adjusted P-value, Genes
        top_n: Number of top terms to display
        title: Plot title

    Returns:
        Plotly Figure object
    """
    if enrichment_df is None or enrichment_df.empty:
        fig = go.Figure()
        fig.update_layout(
            title=title,
            annotations=[dict(
                text="No enrichment results to display",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False, font=dict(size=16)
            )]
        )
        return fig

    df = enrichment_df.head(top_n).copy()

    # Parse gene count from Overlap column (e.g., "5/100" → 5)
    if 'Overlap' in df.columns:
        df['gene_count'] = df['Overlap'].astype(str).apply(
            lambda x: int(x.split('/')[0]) if '/' in str(x) else 0
        )
    else:
        df['gene_count'] = 10  # Default size

    # Compute -log10(adjusted p-value)
    padj_col = 'Adjusted P-value' if 'Adjusted P-value' in df.columns else 'adjusted_p_value'
    if padj_col not in df.columns:
        # Fallback: use first column that looks like a p-value
        for col in df.columns:
            if 'adj' in col.lower() and 'p' in col.lower():
                padj_col = col
                break
        else:
            padj_col = 'P-value' if 'P-value' in df.columns else df.columns[2]

    df['-log10_padj'] = -np.log10(df[padj_col].astype(float).clip(lower=1e-300))

    # Term column
    term_col = 'Term' if 'Term' in df.columns else df.columns[0]

    # Truncate long term names
    df['term_display'] = df[term_col].astype(str).apply(
        lambda x: x[:60] + '...' if len(x) > 60 else x
    )

    # Sort by significance (most significant at top)
    df = df.sort_values('-log10_padj', ascending=True)

    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df['-log10_padj'],
        y=df['term_display'],
        mode='markers',
        marker=dict(
            size=df['gene_count'].clip(lower=5, upper=40),
            color=df['-log10_padj'],
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title='-log₁₀(Adj. P)'),
            line=dict(width=1, color='DarkSlateGrey')
        ),
        text=[f"Genes: {g}" for g in df.get('Genes', [''] * len(df))],
        hovertemplate=(
            '<b>%{y}</b><br>'
            '-log₁₀(padj): %{x:.2f}<br>'
            'Gene count: %{marker.size}<br>'
            '%{text}<extra></extra>'
        )
    ))

    fig.update_layout(
        title=title,
        xaxis_title='-log₁₀(Adjusted P-value)',
        yaxis_title='',
        height=max(400, len(df) * 25 + 100),
        margin=dict(l=300),
        showlegend=False,
    )

    return fig


def create_deconvolution_plot(
    proportions_df: pd.DataFrame,
    sample_conditions: Optional[Dict[str, str]] = None,
) -> go.Figure:
    """
    Create stacked bar chart of cell type proportions.

    Args:
        proportions_df: samples × cell types DataFrame (values should sum to ~1 per row)
        sample_conditions: Optional dict mapping sample names to conditions for grouping

    Returns:
        Plotly Figure object
    """
    if proportions_df is None or proportions_df.empty:
        fig = go.Figure()
        fig.update_layout(
            title="Cell Type Deconvolution",
            annotations=[dict(
                text="No deconvolution results to display",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False, font=dict(size=16)
            )]
        )
        return fig

    df = proportions_df.copy()

    # Sort samples by condition if provided
    if sample_conditions:
        df['_condition'] = [sample_conditions.get(s, 'Unknown') for s in df.index]
        df = df.sort_values('_condition')
        df = df.drop(columns=['_condition'])

    colors = px.colors.qualitative.Plotly + px.colors.qualitative.Set2

    fig = go.Figure()

    for i, cell_type in enumerate(df.columns):
        fig.add_trace(go.Bar(
            name=cell_type,
            x=df.index,
            y=df[cell_type],
            marker_color=colors[i % len(colors)],
            hovertemplate=(
                f'<b>{cell_type}</b><br>'
                'Sample: %{x}<br>'
                'Proportion: %{y:.3f}<extra></extra>'
            )
        ))

    fig.update_layout(
        barmode='stack',
        title='Cell Type Deconvolution - Estimated Proportions',
        xaxis_title='Sample',
        yaxis_title='Proportion',
        yaxis=dict(range=[0, 1.05]),
        legend_title='Cell Type',
        height=500,
    )

    return fig


def create_gene_expression_plot(
    expression_df: pd.DataFrame,
    gene: str,
    sample_conditions: Dict[str, str],
    plot_type: str = "box",
) -> go.Figure:
    """
    Create box or violin plot for a single gene's expression across conditions.

    Args:
        expression_df: samples (rows) × genes (columns)
        gene: Gene name (must be a column in expression_df)
        sample_conditions: Dict mapping sample_name → condition
        plot_type: "box" or "violin"

    Returns:
        Plotly Figure object
    """
    if gene not in expression_df.columns:
        raise ValueError(
            f"Gene '{gene}' not found in expression data. "
            f"Available genes (first 5): {list(expression_df.columns[:5])}"
        )

    gene_values = expression_df[gene]
    conditions = [sample_conditions.get(s, "Unknown") for s in expression_df.index]
    samples = expression_df.index.tolist()

    colors = px.colors.qualitative.Set2
    unique_conditions = sorted(set(conditions))

    fig = go.Figure()

    for i, cond in enumerate(unique_conditions):
        mask = [c == cond for c in conditions]
        vals = gene_values[mask].values
        sample_names = [s for s, m in zip(samples, mask) if m]
        color = colors[i % len(colors)]

        if plot_type == "violin":
            fig.add_trace(go.Violin(
                y=vals,
                name=cond,
                points="all",
                marker_color=color,
                line_color=color,
                text=sample_names,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "Expression: %{y:.2f}<extra></extra>"
                ),
            ))
        else:
            fig.add_trace(go.Box(
                y=vals,
                name=cond,
                boxpoints="all",
                jitter=0.3,
                marker_color=color,
                line_color=color,
                text=sample_names,
                hovertemplate=(
                    "<b>%{text}</b><br>"
                    "Expression: %{y:.2f}<extra></extra>"
                ),
            ))

    fig.update_layout(
        title=f"Expression of {gene}",
        xaxis_title="Condition",
        yaxis_title="Expression (log2)",
        showlegend=True,
    )

    return fig


def create_venn_diagram(
    gene_sets: Dict[str, set],
    title: str = "DE Gene Overlap",
) -> go.Figure:
    """
    Create Venn diagram for 2-3 gene sets, or UpSet-style plot for 4+ sets.

    Args:
        gene_sets: Dict mapping comparison_name → set of gene names
        title: Plot title

    Returns:
        Plotly Figure object
    """
    set_names = list(gene_sets.keys())
    sets = [gene_sets[k] for k in set_names]
    n = len(set_names)

    if n < 2:
        raise ValueError("At least 2 gene sets are required for a Venn diagram.")

    if n == 2:
        a, b = sets
        only_a = len(a - b)
        only_b = len(b - a)
        inter = len(a & b)

        fig = go.Figure()
        fig.update_layout(
            title=title,
            shapes=[
                dict(type="circle", x0=0.1, y0=0.15, x1=0.55, y1=0.85,
                     line=dict(color="royalblue", width=2),
                     fillcolor="rgba(65,105,225,0.2)"),
                dict(type="circle", x0=0.45, y0=0.15, x1=0.9, y1=0.85,
                     line=dict(color="crimson", width=2),
                     fillcolor="rgba(220,20,60,0.2)"),
            ],
            annotations=[
                dict(x=0.25, y=0.5, text=str(only_a), showarrow=False,
                     font=dict(size=18, color="royalblue")),
                dict(x=0.5, y=0.5, text=str(inter), showarrow=False,
                     font=dict(size=18, color="purple")),
                dict(x=0.75, y=0.5, text=str(only_b), showarrow=False,
                     font=dict(size=18, color="crimson")),
                dict(x=0.25, y=0.9, text=set_names[0], showarrow=False,
                     font=dict(size=13)),
                dict(x=0.75, y=0.9, text=set_names[1], showarrow=False,
                     font=dict(size=13)),
            ],
            xaxis=dict(visible=False, range=[0, 1]),
            yaxis=dict(visible=False, range=[0, 1], scaleanchor="x"),
            width=550, height=450,
        )
        return fig

    if n == 3:
        a, b, c = sets
        abc = a & b & c
        ab_only = (a & b) - c
        ac_only = (a & c) - b
        bc_only = (b & c) - a
        a_only = a - b - c
        b_only = b - a - c
        c_only = c - a - b

        fig = go.Figure()
        fig.update_layout(
            title=title,
            shapes=[
                dict(type="circle", x0=0.15, y0=0.25, x1=0.6, y1=0.9,
                     line=dict(color="royalblue", width=2),
                     fillcolor="rgba(65,105,225,0.15)"),
                dict(type="circle", x0=0.4, y0=0.25, x1=0.85, y1=0.9,
                     line=dict(color="crimson", width=2),
                     fillcolor="rgba(220,20,60,0.15)"),
                dict(type="circle", x0=0.27, y0=0.1, x1=0.73, y1=0.65,
                     line=dict(color="forestgreen", width=2),
                     fillcolor="rgba(34,139,34,0.15)"),
            ],
            annotations=[
                dict(x=0.28, y=0.72, text=str(len(a_only)), showarrow=False,
                     font=dict(size=15, color="royalblue")),
                dict(x=0.72, y=0.72, text=str(len(b_only)), showarrow=False,
                     font=dict(size=15, color="crimson")),
                dict(x=0.5, y=0.22, text=str(len(c_only)), showarrow=False,
                     font=dict(size=15, color="forestgreen")),
                dict(x=0.5, y=0.75, text=str(len(ab_only)), showarrow=False,
                     font=dict(size=14, color="purple")),
                dict(x=0.35, y=0.42, text=str(len(ac_only)), showarrow=False,
                     font=dict(size=14, color="teal")),
                dict(x=0.65, y=0.42, text=str(len(bc_only)), showarrow=False,
                     font=dict(size=14, color="orangered")),
                dict(x=0.5, y=0.52, text=str(len(abc)), showarrow=False,
                     font=dict(size=14, color="black")),
                dict(x=0.2, y=0.93, text=set_names[0], showarrow=False,
                     font=dict(size=12)),
                dict(x=0.8, y=0.93, text=set_names[1], showarrow=False,
                     font=dict(size=12)),
                dict(x=0.5, y=0.07, text=set_names[2], showarrow=False,
                     font=dict(size=12)),
            ],
            xaxis=dict(visible=False, range=[0, 1]),
            yaxis=dict(visible=False, range=[0, 1], scaleanchor="x"),
            width=600, height=550,
        )
        return fig

    # 4+ sets: UpSet-style plot
    from itertools import combinations

    all_elements = set()
    for s in sets:
        all_elements |= s

    intersections = []
    for r in range(1, n + 1):
        for combo in combinations(range(n), r):
            combo_set = set(all_elements)
            for i in range(n):
                if i in combo:
                    combo_set &= sets[i]
                else:
                    combo_set -= sets[i]
            if combo_set:
                intersections.append((combo, len(combo_set)))

    intersections.sort(key=lambda x: x[1], reverse=True)

    combo_labels = [" & ".join(set_names[i] for i in combo) for combo, _ in intersections]
    sizes = [size for _, size in intersections]

    fig = make_subplots(
        rows=2, cols=1, row_heights=[0.65, 0.35],
        shared_xaxes=True, vertical_spacing=0.02,
    )

    fig.add_trace(
        go.Bar(
            x=list(range(len(sizes))),
            y=sizes,
            marker_color="steelblue",
            hovertemplate="<b>%{customdata}</b><br>Count: %{y}<extra></extra>",
            customdata=combo_labels,
            showlegend=False,
        ),
        row=1, col=1,
    )

    for row_idx, name in enumerate(set_names):
        member_xs = [col_idx for col_idx, (combo, _) in enumerate(intersections) if row_idx in combo]
        fig.add_trace(
            go.Scatter(
                x=member_xs,
                y=[row_idx] * len(member_xs),
                mode="markers",
                marker=dict(size=10, color="steelblue"),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=2, col=1,
        )
        non_member_xs = [col_idx for col_idx, (combo, _) in enumerate(intersections) if row_idx not in combo]
        fig.add_trace(
            go.Scatter(
                x=non_member_xs,
                y=[row_idx] * len(non_member_xs),
                mode="markers",
                marker=dict(size=10, color="lightgray"),
                showlegend=False,
                hoverinfo="skip",
            ),
            row=2, col=1,
        )

    for col_idx, (combo, _) in enumerate(intersections):
        if len(combo) > 1:
            members = sorted(combo)
            fig.add_trace(
                go.Scatter(
                    x=[col_idx, col_idx],
                    y=[members[0], members[-1]],
                    mode="lines",
                    line=dict(color="steelblue", width=2),
                    showlegend=False,
                    hoverinfo="skip",
                ),
                row=2, col=1,
            )

    fig.update_layout(
        title=title,
        height=400 + n * 30,
        width=max(500, len(intersections) * 45 + 150),
    )
    fig.update_yaxes(title_text="Intersection Size", row=1, col=1)
    fig.update_yaxes(
        tickvals=list(range(n)),
        ticktext=set_names,
        row=2, col=1,
    )
    fig.update_xaxes(visible=False, row=1, col=1)
    fig.update_xaxes(visible=False, row=2, col=1)

    return fig


def create_normalization_comparison_plot(
    raw_counts: pd.DataFrame,
    normalized_counts: pd.DataFrame,
    sample_conditions: Dict[str, str],
) -> go.Figure:
    """
    Create side-by-side box plots comparing raw and normalized count distributions.

    Args:
        raw_counts: genes × samples DataFrame of raw counts
        normalized_counts: genes × samples DataFrame of normalized counts
        sample_conditions: Dict mapping sample_name → condition

    Returns:
        Plotly Figure object
    """
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Raw (log2)", "Normalized"],
        shared_yaxes=True,
    )

    colors = px.colors.qualitative.Set2
    conditions = sorted(set(sample_conditions.values()))
    cond_color = {c: colors[i % len(colors)] for i, c in enumerate(conditions)}

    log2_raw = np.log2(raw_counts + 1)
    for sample in log2_raw.columns:
        cond = sample_conditions.get(sample, "Unknown")
        fig.add_trace(
            go.Box(
                y=log2_raw[sample].values,
                name=sample,
                marker_color=cond_color.get(cond, "gray"),
                legendgroup=cond,
                legendgrouptitle_text=cond,
                showlegend=True,
                hovertemplate="Sample: " + sample + "<br>Value: %{y:.2f}<extra></extra>",
            ),
            row=1, col=1,
        )

    for sample in normalized_counts.columns:
        cond = sample_conditions.get(sample, "Unknown")
        fig.add_trace(
            go.Box(
                y=normalized_counts[sample].values,
                name=sample,
                marker_color=cond_color.get(cond, "gray"),
                legendgroup=cond,
                legendgrouptitle_text=cond,
                showlegend=False,
                hovertemplate="Sample: " + sample + "<br>Value: %{y:.2f}<extra></extra>",
            ),
            row=1, col=2,
        )

    fig.update_layout(
        title="Normalization Comparison",
        height=500,
        showlegend=True,
    )

    return fig
