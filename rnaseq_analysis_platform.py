"""
RNA-seq Analysis Platform
A Streamlit-based application for RNA-seq data analysis and visualization.
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from typing import Dict, List, Tuple, Optional
import io
import os
import tempfile

# Import backend modules
from rnaseq_parser import RNASeqParser, DataType, ParseResult, ParserValidationError
from de_analysis import DEAnalysisEngine, DEResult
from pathway_enrichment import PathwayEnrichment
from visualizations import (
    create_volcano_plot,
    create_clustered_heatmap,
    create_pca_plot,
    create_ma_plot,
    create_correlation_heatmap,
    compute_de_summary
)
from gene_panels import GenePanelAnalyzer
from export_engine import ExportEngine, ExportData, EnrichmentResult
from demo_data import load_demo_dataset, get_demo_description
from gene_search import search_genes, get_gene_summary
from qc_plots import (
    create_library_size_barplot,
    create_count_distribution_boxplot,
    create_gene_detection_plot,
    create_sample_similarity_heatmap
)
from interpretation_engine import InterpretationEngine
from cell_deconvolution import CellTypeDeconvolution
from advanced_qc import AdvancedQC
from visualizations import create_enrichment_dotplot, create_deconvolution_plot
from visualizations import create_gene_expression_plot, create_venn_diagram, create_normalization_comparison_plot
from sklearn.decomposition import PCA
import numpy as np
from multi_file_merger import merge_parse_results, validate_merge_compatibility

PLATFORM_CSS = """
/* === RNA-seq Analysis Platform Professional Theme === */

/* Global Settings */
:root {
    --primary-color: #0f3460;
    --accent-color: #16a085;
    --background-color: #f8f9fa;
    --text-color: #2c3e50;
}

/* App Background */
.stApp {
    background-color: #f8f9fa;
    color: #2c3e50;
}

/* Sidebar Styling */
[data-testid="stSidebar"] {
    background-color: #0f3460;
    background-image: linear-gradient(180deg, #0f3460 0%, #1a1a2e 100%);
}
[data-testid="stSidebar"] [data-testid="stMarkdown"] {
    color: #ffffff;
}
[data-testid="stSidebar"] .stSlider label,
[data-testid="stSidebar"] .stSelectbox label,
[data-testid="stSidebar"] .stMultiSelect label,
[data-testid="stSidebar"] .stTextInput label {
    color: #e0e0e0 !important;
}
[data-testid="stSidebar"] hr {
    border-color: rgba(255,255,255,0.2);
}

/* Headers */
h1, h2, h3 {
    color: #0f3460 !important;
    font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
    font-weight: 700;
    letter-spacing: -0.5px;
}

/* Step Indicator Styling */
.step-indicator {
    display: flex;
    flex-direction: column;
    gap: 8px;
    padding: 10px 0;
}
.step {
    padding: 8px 16px;
    border-radius: 20px;
    font-size: 0.9rem;
    font-weight: 500;
    transition: all 0.3s ease;
    display: flex;
    align-items: center;
    gap: 10px;
}
.step.completed {
    background-color: rgba(22, 160, 133, 0.15);
    color: #16a085;
    border: 1px solid rgba(22, 160, 133, 0.2);
}
.step.active {
    background-color: #16a085;
    color: white;
    box-shadow: 0 4px 6px rgba(22, 160, 133, 0.2);
    font-weight: 600;
}
.step.pending {
    background-color: rgba(255,255,255,0.05);
    color: rgba(255,255,255,0.5);
}

/* Metric Cards */
[data-testid="stMetric"] {
    background-color: white;
    border: 1px solid #e0e0e0;
    border-radius: 8px;
    padding: 16px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.05);
    transition: transform 0.2s;
}
[data-testid="stMetric"]:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 8px rgba(0,0,0,0.1);
}
[data-testid="stMetric"] label {
    color: #7f8c8d;
    font-size: 0.85rem;
    text-transform: uppercase;
    letter-spacing: 0.5px;
}
[data-testid="stMetricValue"] {
    color: #0f3460;
    font-weight: 700;
}

/* Tabs */
.stTabs [data-baseweb="tab-list"] {
    gap: 8px;
    background-color: transparent;
    padding-bottom: 8px;
}
.stTabs [data-baseweb="tab"] {
    background-color: white;
    border-radius: 6px;
    padding: 8px 16px;
    border: 1px solid #e0e0e0;
    color: #7f8c8d;
    transition: all 0.2s;
}
.stTabs [aria-selected="true"] {
    background-color: #16a085;
    color: white !important;
    border-color: #16a085;
    font-weight: 600;
}

/* Buttons */
.stButton > button[kind="primary"] {
    background-color: #16a085;
    color: white;
    border: none;
    border-radius: 6px;
    padding: 0.6rem 1.2rem;
    font-weight: 600;
    letter-spacing: 0.3px;
    transition: all 0.2s;
}
.stButton > button[kind="primary"]:hover {
    background-color: #138d75;
    box-shadow: 0 4px 12px rgba(22, 160, 133, 0.3);
    transform: translateY(-1px);
}

/* Expanders */
.streamlit-expanderHeader {
    background-color: #f1f3f5;
    border-radius: 6px;
    color: #2c3e50;
    font-weight: 500;
}

/* DataFrames & Tables */
[data-testid="stDataFrame"] {
    border: 1px solid #e0e0e0;
    border-radius: 6px;
    overflow: hidden;
}

/* Alerts/Notifications */
.stAlert {
    border-radius: 8px;
    border: none;
    box-shadow: 0 2px 4px rgba(0,0,0,0.05);
}
"""

# --- Page Config ---
st.set_page_config(
    page_title="RNA-seq Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)
st.markdown(f"<style>{PLATFORM_CSS}</style>", unsafe_allow_html=True)

# --- Session State Initialization ---
if "parsed_result" not in st.session_state:
    st.session_state["parsed_result"] = None
if "expression_matrix" not in st.session_state:
    st.session_state["expression_matrix"] = None
if "sample_metadata" not in st.session_state:
    st.session_state["sample_metadata"] = {}
if "conditions" not in st.session_state:
    st.session_state["conditions"] = ["Control", "Treatment"]
if "comparisons" not in st.session_state:
    st.session_state["comparisons"] = []
if "de_results" not in st.session_state:
    st.session_state["de_results"] = {}
if "enrichment_results" not in st.session_state:
    st.session_state["enrichment_results"] = {}
if "figures" not in st.session_state:
    st.session_state["figures"] = {}
if "active_comparison" not in st.session_state:
    st.session_state["active_comparison"] = None
if "analysis_complete" not in st.session_state:
    st.session_state["analysis_complete"] = False
if "current_step" not in st.session_state:
    st.session_state["current_step"] = 1
if "gsea_results" not in st.session_state:
    st.session_state["gsea_results"] = {}
if "interpretations" not in st.session_state:
    st.session_state["interpretations"] = {}
if "deconvolution_results" not in st.session_state:
    st.session_state["deconvolution_results"] = None
if "outlier_info" not in st.session_state:
    st.session_state["outlier_info"] = None
if "batch_effects" not in st.session_state:
    st.session_state["batch_effects"] = None
if "pca_coordinates" not in st.session_state:
    st.session_state["pca_coordinates"] = None

# --- Help Text Constants ---
HELP_TEXTS = {
    "padj_threshold": "Adjusted p-value cutoff. Genes with padj below this are considered significant. Default 0.05 controls false discovery rate at 5%.",
    "lfc_threshold": "Log2 fold change cutoff. Genes must exceed this magnitude to be called up/down regulated. 1.0 = 2-fold change.",
    "summary": """**Purpose**: Provides a quick "health check" of your analysis results and a high-level overview of the biological response.

**How to Read This**: Review the total number of significant genes and the balance between up- and down-regulated genes. The top genes tables show the most extreme changes.

**What to Look For**:
- A reasonable number of DE genes (too few might mean low power; too many might mean a very strong effect or batch issues).
- Top genes should make biological sense for your experimental model.
- Balance between upregulated and downregulated genes.

**Common Pitfalls**:
- ⚠️ Don't focus only on the top 10 genes; they may be highly variable or have low expression despite high fold changes.
""",
    "de_results": """**Purpose**: The master table containing every gene's statistical performance and expression metrics.

**How to Read This**:
- **gene**: The identifier (Symbol or Ensembl ID).
- **log2FoldChange**: Magnitude of change (positive = up, negative = down).
- **padj**: Adjusted p-value (significance).
- **baseMean**: Average expression level across all samples.

**What to Look For**:
- Sort by `padj` to find the most statistically robust changes.
- Sort by `log2FoldChange` to find the genes with the largest biological shifts.
- Use the search box to find specific genes of interest.

**Common Pitfalls**:
- ⚠️ High fold changes in genes with very low `baseMean` are often unreliable and may be technical noise.
""",
    "volcano": """**Purpose**: Shows which genes are both statistically significant AND biologically meaningful.

**How to Read This**: Each dot is one gene. The x-axis shows the log2 fold change (magnitude of change), and the y-axis shows the statistical significance (-log10 adjusted p-value).

| |Log2FC| | Fold Change |
| :--- | :--- | :--- |
| |LFC| = 1 | 2-fold |
| |LFC| = 2 | 4-fold |
| |LFC| = 3 | 8-fold |

**What to Look For**:
- Genes in the **upper-right corner** are significantly upregulated.
- Genes in the **upper-left corner** are significantly downregulated.
- The "upper corners" pattern indicates a strong biological response.

**Common Pitfalls**:
- ⚠️ A very low p-value does NOT mean a gene is biologically important; it just means the change is consistent.
- Significance ≠ Importance.
""",
    "ma_plot": """**Purpose**: Visualizes the relationship between average expression level and fold change to detect systematic biases.

**How to Read This**: The x-axis represents the average expression (A), and the y-axis represents the log2 fold change (M).

**What to Look For**:
- Points should be centered around y=0 across the entire expression range.
- A "banana" shape or vertical shift indicates normalization bias or technical artifacts.
- Higher variance is typically expected at lower expression levels (left side).

**Common Pitfalls**:
- ⚠️ Be cautious of large fold changes in genes with very low expression, as these are often noisy and expression-dependent.
""",
    "pca": """**Purpose**: Reduces the complexity of thousands of genes into a 2D map to see how samples relate to each other.

**How to Read This**: Samples that are closer together have more similar overall gene expression profiles. The axes (PC1, PC2) represent the directions of greatest variation.

**What to Look For**:
- **Replicate tightness**: Biological replicates should cluster closely together.
- **Condition separation**: Different experimental groups should form distinct clusters.
- **Variance explained**: The % on each axis shows how much data variation is captured.

**Common Pitfalls**:
- ⚠️ If PC1 is dominated by a batch effect (e.g., processing date) rather than your treatment, the biological signal may be masked.
""",
    "heatmap": """**Purpose**: Provides a high-level view of expression patterns across samples for the most significant genes.

**How to Read This**: Colors represent **Z-scores** (standard deviations from the gene's mean). Red indicates higher expression, and blue indicates lower. Dendrograms show hierarchical clustering.

**What to Look For**:
- **Gene blocks**: Groups of genes that are consistently up or down in specific conditions.
- **Sample clustering**: Samples should ideally cluster by their experimental condition.
- **Outliers**: Individual samples that don't match their group's pattern.

**Common Pitfalls**:
- ⚠️ If samples cluster by batch or processing date instead of biology, your results may be confounded by technical variables.
""",
    "enrichment": """**Purpose**: Identifies biological "themes" or pathways that are overrepresented in your list of differentially expressed genes.

**How to Read This**: This uses **Over-Representation Analysis (ORA)**. In the dot plot, size represents the number of genes in that pathway, and color represents statistical significance.

**What to Look For**:
- Pathways relevant to your biological system (e.g., "Skin Development" or "Inflammatory Response").
- High "Gene Ratio" combined with low p-values.
- Consistency across different databases (GO, KEGG, Reactome).

**Common Pitfalls**:
- ⚠️ If your DE gene list is very large (thousands of genes), enrichment results can become non-specific and lose biological resolution.
""",
    "gsea": """**Purpose**: Detects subtle but consistent changes in entire pathways without relying on arbitrary significance cutoffs.

**How to Read This**: Unlike ORA, **GSEA** uses a ranked list of ALL genes. The **Normalized Enrichment Score (NES)** indicates if a pathway is generally up (positive) or down (negative).

**What to Look For**:
- **NES magnitude**: High absolute NES values (>1.5 or 2.0) indicate strong pathway-level shifts.
- **Leading Edge**: The subset of genes that actually drive the enrichment score.
- **FDR**: Focus on pathways with a False Discovery Rate < 0.25.

**Common Pitfalls**:
- ⚠️ GSEA results depend heavily on the quality and completeness of the gene set annotations in databases like MSigDB.
""",
    "interpretation": """**Purpose**: Automated insights to help bridge the gap between statistical results and biological meaning.

**How to Read This**: These insights are generated by analyzing your DE and enrichment results to highlight key findings and potential issues.

**What to Look For**:
- **Critical alerts**: Red flags about data quality or unexpected results.
- **Pathway clusters**: Groups of related pathways that show a coordinated response.
- **Recommendations**: Suggested follow-up experiments or validation targets.

**Common Pitfalls**:
- ⚠️ These are automated suggestions and should always be verified by your own biological expertise and literature review.
""",
    "deconvolution": """**Purpose**: Estimates the proportions of different cell types within your "bulk" tissue sample.

**How to Read This**: Uses the **NNLS (Non-Negative Least Squares)** method to compare your data against "fingerprints" of known cell types. Results show estimated percentages.

**What to Look For**:
- **Population shifts**: Changes in cell types (e.g., increased immune cells in inflamed skin).
- **Biological consistency**: Does the cell type mix match the tissue you sampled?
- **Proportions**: Ensure the estimated proportions align with your experimental expectations.

**Common Pitfalls**:
- ⚠️ These are computational estimates, not ground truth. They are highly dependent on the accuracy of the reference signature used.
""",
    "qc_advanced": """**Purpose**: Deep-dive quality control to ensure your data is reliable for publication and downstream analysis.

**How to Read This**: Focuses on identifying outlier samples and detecting potential technical artifacts like batch effects.

**What to Look For**:
- **Outliers**: Samples that fall far outside the expected clusters in PCA space (e.g., >3 standard deviations).
- **Batch effects**: Patterns that correlate with technical variables (e.g., plate ID) rather than biology.
- **Sample distances**: Quantitative metrics of how similar replicates are to each other.

**Common Pitfalls**:
- ⚠️ Removing outliers should be done with caution and only when there is a clear technical justification.
""",
    "gene_panels": """**Purpose**: Focused analysis on curated sets of genes known to be important in specific biological contexts (e.g., Dermatology).

**How to Read This**: These panels allow you to bypass the "noise" of the whole transcriptome and look directly at your genes of interest.

**What to Look For**:
- **Coordinated changes**: Do genes in a functional panel (e.g., "Skin Barrier") move together?
- **Marker validation**: Check for expected changes in known markers for your treatment.
- **Panel-specific patterns**: Unique signatures that define your experimental condition.

**Common Pitfalls**:
- ⚠️ A panel only shows a subset of genes; don't ignore significant changes that might fall outside these curated lists.
""",
    "gene_expression": "View normalized expression of individual genes across conditions. Useful for validating DE results.",
    "venn": "Overlap between DE gene sets from different comparisons. Shared genes respond to multiple conditions.",
    "upload": "Upload a gene count matrix (CSV/TSV/Excel). Rows = genes, columns = samples. Raw counts preferred for DE analysis.",
    "metadata": "Assign each sample to an experimental condition (e.g., Control, Treatment). Need ≥2 samples per condition.",
    "export": """**Purpose**: Save your results for further analysis, presentations, or publication.

**How to Read This**: Choose the format that best fits your next steps.

**What to Look For**:
- **Excel**: Best for sharing raw data and tables with collaborators.
- **PDF/Images**: Best for high-quality, publication-ready figures (PNG/SVG).
- **Session File**: Save your entire analysis state to resume later.

**Common Pitfalls**:
- ⚠️ Ensure you export the results for the correct comparison if you have multiple treatment groups.
""",
}

# --- Helper Functions ---

def set_step(step: int):
    st.session_state["current_step"] = step
    st.rerun()

def get_unique_conditions(sample_metadata: Dict[str, Dict]) -> List[str]:
    """Get list of unique conditions from metadata."""
    if not sample_metadata:
        return []
    conditions = set(m["condition"] for m in sample_metadata.values())
    return sorted(list(conditions))

@st.cache_data(show_spinner=False)
def cached_de_analysis(_counts_df, _meta_df, comparisons_tuple, design_factor):
    """Cache DE analysis results. Uses _ prefix for unhashable DataFrames."""
    engine = DEAnalysisEngine()
    return engine.run_all_comparisons(_counts_df, _meta_df, list(comparisons_tuple), design_factor)

@st.cache_data(show_spinner=False)
def cached_enrichment(genes_tuple: tuple, database: str = "go_bp"):
    """Cache pathway enrichment API calls."""
    pe = PathwayEnrichment()
    if database == "go_bp":
        return pe.get_go_enrichment(list(genes_tuple))
    elif database == "kegg":
        return pe.get_kegg_enrichment(list(genes_tuple))
    elif database == "go_mf":
        return pe.get_go_mf_enrichment(list(genes_tuple))
    elif database == "go_cc":
        return pe.get_go_cc_enrichment(list(genes_tuple))
    elif database == "reactome":
        return pe.get_reactome_enrichment(list(genes_tuple))
    elif database == "wikipathway":
        return pe.get_wikipathway_enrichment(list(genes_tuple))
    elif database == "hallmark":
        return pe.get_hallmark_enrichment(list(genes_tuple))
    return None

def count_samples_per_condition(sample_metadata: Dict[str, Dict]) -> Dict[str, int]:
    """Count how many samples are assigned to each condition."""
    counts = {}
    for meta in sample_metadata.values():
        condition = meta["condition"]
        counts[condition] = counts.get(condition, 0) + 1
    return counts

def plot_with_download(fig: go.Figure, key: str, title: str = ""):
    """Display Plotly chart with PNG/SVG download buttons."""
    st.plotly_chart(fig, use_container_width=True, key=f"chart_{key}")
    col1, col2, _ = st.columns([1, 1, 6])
    with col1:
        try:
            png_bytes = fig.to_image(format="png", scale=2, width=1200, height=800)
            st.download_button("📷 PNG", png_bytes, f"{key}.png", "image/png", key=f"dl_png_{key}")
        except Exception:
            pass
    with col2:
        try:
            svg_bytes = fig.to_image(format="svg", width=1200, height=800)
            st.download_button("📐 SVG", svg_bytes, f"{key}.svg", "image/svg+xml", key=f"dl_svg_{key}")
        except Exception:
            pass

def validate_metadata_assignment(
    sample_metadata: Dict[str, Dict], counts_df: pd.DataFrame
) -> Tuple[bool, List[str]]:
    """
    Validate metadata assignment.
    Checks:
    - All samples assigned
    - Max 10 conditions (guardrail)
    - Min 2 samples per condition
    """
    errors = []

    # Check all samples assigned
    unassigned = [s for s in counts_df.index if s not in sample_metadata]
    if unassigned:
        errors.append(f"Unassigned samples: {', '.join(unassigned)}")

    # Get unique conditions
    conditions = get_unique_conditions(sample_metadata)

    # Guardrail: Max 10 conditions
    if len(conditions) > 10:
        errors.append(
            "Maximum 10 conditions supported. Please reduce the number of unique conditions."
        )

    # Check minimum 2 samples per condition
    sample_counts = count_samples_per_condition(sample_metadata)
    for condition, count in sample_counts.items():
        if count < 2:
            errors.append(
                f"Condition '{condition}' has only {count} sample(s), need ≥2 for DE analysis"
            )

    return len(errors) == 0, errors

def validate_comparison(
    comp: Tuple[str, str], available_conditions: List[str]
) -> Tuple[bool, str]:
    """Validate a comparison pair."""
    test, ref = comp
    if test == ref:
        return False, "Test and Reference conditions must be different."
    if test not in available_conditions or ref not in available_conditions:
        return False, "Selected conditions are not valid."
    return True, ""

def add_comparison(test: str, ref: str):
    """Add a comparison to session state."""
    comp = (test, ref)
    if comp not in st.session_state["comparisons"]:
        st.session_state["comparisons"].append(comp)

def remove_comparison(comp: Tuple[str, str]):
    """Remove a comparison from session state."""
    if comp in st.session_state["comparisons"]:
        st.session_state["comparisons"].remove(comp)

# --- Sidebar ---
with st.sidebar:
    st.title("🧬 RNA-seq Analysis")
    st.markdown("---")
    
    # Styled step indicator
    current_step = st.session_state["current_step"]
    steps = ["Upload Data", "Configure", "Analyze", "Results"]
    steps_html = '<div class="step-indicator">'
    for i, step_name in enumerate(steps, 1):
        if i < current_step:
            steps_html += f'<div class="step completed">✅ {step_name}</div>'
        elif i == current_step:
            steps_html += f'<div class="step active">● {step_name}</div>'
        else:
            steps_html += f'<div class="step pending">○ {step_name}</div>'
    steps_html += '</div>'
    st.markdown(steps_html, unsafe_allow_html=True)
    
    st.markdown("---")
    # Global settings only
    st.subheader("Settings")
    padj_threshold = st.slider("P-adj threshold", 0.001, 0.1, 0.05, key="padj_thresh", help=HELP_TEXTS["padj_threshold"])
    lfc_threshold = st.slider("Log2FC threshold", 0.5, 3.0, 1.0, key="lfc_thresh", help=HELP_TEXTS["lfc_threshold"])

    st.markdown("---")
    st.subheader("💾 Session")
    if st.session_state.get("analysis_complete"):
        from session_manager import SessionManager
        from datetime import datetime
        json_str = SessionManager.save_session(dict(st.session_state))
        st.download_button(
            "Save Session",
            json_str.encode("utf-8"),
            f"rnaseq_session_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
            "application/json",
            key="save_session"
        )
    session_file = st.file_uploader("Load Session", type=["json"], key="load_session")
    if session_file:
        from session_manager import SessionManager
        loaded = SessionManager.load_session(session_file.read().decode("utf-8"))
        for key, val in loaded.items():
            if key in st.session_state:
                st.session_state[key] = val
        st.success("Session loaded! Re-upload data file to restore full analysis.")
        st.rerun()

# --- Main Area ---

# STEP 1: Upload & Preview
if current_step == 1:
    st.header("Step 1: Upload Data")
    st.markdown("Upload your RNA-seq count matrix or load demo data to explore the platform.")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        uploaded_files = st.file_uploader(
            "Upload Count Files (CSV/TSV/Excel) — upload multiple files to merge",
            type=["csv", "tsv", "txt", "xlsx"],
            accept_multiple_files=True,
            help=HELP_TEXTS["upload"]
        )
    with col2:
        st.markdown("**Or try the demo:**")
        if st.button("🧪 Load Demo Data", use_container_width=True):
            counts_df, metadata_df = load_demo_dataset()
            # Create a ParseResult-like object for compatibility
            demo_result = ParseResult(
                expression_df=counts_df,
                normalized_df=None,
                de_results_df=None,
                data_type=DataType.RAW_COUNTS,
                can_run_de=True,
                warnings=["Demo dataset loaded"],
                dropped_columns=[],
                gene_column_source="demo",
                needs_user_input=False,
                gene_column_candidates=[],
                data_types_detected=[DataType.RAW_COUNTS]
            )
            st.session_state["parsed_result"] = demo_result
            st.session_state["expression_matrix"] = counts_df
            st.session_state["sample_metadata"] = {}
            # Pre-populate metadata from demo
            for sample in counts_df.index:
                cond = metadata_df.loc[sample, "condition"]
                st.session_state["sample_metadata"][sample] = {"condition": cond}
            st.session_state["conditions"] = list(metadata_df["condition"].unique())
            st.session_state["comparisons"] = []
            st.session_state["analysis_complete"] = False
            st.rerun()

    # Handle file upload
    if uploaded_files:
        # Handle single or multiple files
        if len(uploaded_files) == 1:
            # Single file — existing behavior
            uploaded_file = uploaded_files[0]
            tmp_path = None
            try:
                suffix = f".{uploaded_file.name.split('.')[-1]}"
                with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                    tmp.write(uploaded_file.getvalue())
                    tmp_path = tmp.name
                parser = RNASeqParser()
                result = parser.parse(tmp_path)

                is_new_file = False
                if st.session_state["parsed_result"] is None:
                    is_new_file = True
                elif hasattr(st.session_state["parsed_result"], "expression_df"):
                    old_df = st.session_state["parsed_result"].expression_df
                    new_df = result.expression_df
                    if old_df is None and new_df is None:
                        is_new_file = False
                        old_de = getattr(st.session_state["parsed_result"], "de_results_df", None)
                        new_de = result.de_results_df
                        if old_de is not None and new_de is not None:
                            is_new_file = not new_de.equals(old_de)
                    elif old_df is None or new_df is None:
                        is_new_file = True
                    else:
                        is_new_file = not new_df.equals(old_df)

                if is_new_file:
                    st.session_state["parsed_result"] = result
                    st.session_state["expression_matrix"] = result.expression_df
                    st.session_state["analysis_complete"] = False
                    st.session_state["sample_metadata"] = {}
                    st.session_state["comparisons"] = []
                    st.session_state["de_results"] = {}
                    st.session_state["enrichment_results"] = {}
                    st.session_state["gsea_results"] = {}
                    st.session_state["interpretations"] = {}
                    st.session_state["deconvolution_results"] = None
                    st.success(f"✅ Loaded {result.data_type.value} data successfully!")

                    if hasattr(result, "data_types_detected") and len(result.data_types_detected) > 1 and result.expression_df is not None:
                        sample_names = result.expression_df.index.tolist()
                        detected_metadata = {}
                        detected_conditions = set()
                        import re
                        for sample in sample_names:
                            match = re.match(r"^.+?_(.+?)_", sample)
                            if match:
                                condition = match.group(1)
                                detected_metadata[sample] = {"condition": condition}
                                detected_conditions.add(condition)
                        if detected_metadata:
                            st.session_state["sample_metadata"] = detected_metadata
                            st.session_state["conditions"] = sorted(list(detected_conditions))
                            st.success(f"✓ Auto-detected conditions: {', '.join(sorted(detected_conditions))}")

            except Exception as e:
                st.error(f"Error parsing file: {str(e)}")
            finally:
                if tmp_path and os.path.exists(tmp_path):
                    os.unlink(tmp_path)

        else:
            # Multiple files — parse each and merge
            st.info(f"📂 Merging {len(uploaded_files)} files...")
            parse_results = []
            filenames = []
            tmp_paths = []
            progress = st.progress(0, text="Parsing files...")

            try:
                for i, uf in enumerate(uploaded_files):
                    progress.progress((i + 1) / len(uploaded_files), text=f"Parsing {uf.name}...")
                    suffix = f".{uf.name.split('.')[-1]}"
                    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                        tmp.write(uf.getvalue())
                        tmp_paths.append(tmp.name)
                    parser = RNASeqParser()
                    r = parser.parse(tmp_paths[-1])
                    parse_results.append(r)
                    filenames.append(uf.name)

                progress.progress(1.0, text="Merging...")

                # Validate compatibility
                validation = validate_merge_compatibility(parse_results)
                report = validation["report"]

                # Show merge report
                with st.expander("📊 Merge Report", expanded=True):
                    rc1, rc2, rc3 = st.columns(3)
                    rc1.metric("Files", report["n_files"])
                    rc2.metric("Total Samples", sum(report["samples_per_file"]))
                    rc3.metric("Gene Overlap", f"{report.get('gene_overlap_pct', 0):.1f}%")

                    if validation["issues"]:
                        for issue in validation["issues"]:
                            st.warning(f"⚠️ {issue}")

                    # Per-file breakdown
                    file_info = []
                    for fname, scount, gcount in zip(filenames, report["samples_per_file"], report["genes_per_file"]):
                        file_info.append({"File": fname, "Samples": scount, "Genes": gcount})
                    st.dataframe(pd.DataFrame(file_info), hide_index=True)

                # Merge
                merged_result = merge_parse_results(parse_results, filenames)

                # Store
                st.session_state["parsed_result"] = merged_result
                st.session_state["expression_matrix"] = merged_result.expression_df
                st.session_state["analysis_complete"] = False
                st.session_state["sample_metadata"] = {}
                st.session_state["comparisons"] = []
                st.session_state["de_results"] = {}
                st.session_state["enrichment_results"] = {}
                st.session_state["gsea_results"] = {}
                st.session_state["interpretations"] = {}
                st.session_state["deconvolution_results"] = None

                if merged_result.expression_df is not None:
                    st.success(
                        f"✅ Merged {len(uploaded_files)} files: "
                        f"{merged_result.expression_df.shape[0]} samples × "
                        f"{merged_result.expression_df.shape[1]} genes"
                    )

                for w in merged_result.warnings:
                    if "dropped" in w.lower() or "overlap" in w.lower():
                        st.warning(f"ℹ️ {w}")

            except ValueError as e:
                st.error(f"❌ Merge failed: {str(e)}")
            except Exception as e:
                st.error(f"Error processing files: {str(e)}")
            finally:
                for tp in tmp_paths:
                    if os.path.exists(tp):
                        os.unlink(tp)

    # If data loaded, show preview
    if st.session_state["parsed_result"]:
        result = st.session_state["parsed_result"]
        # Data summary metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Samples", result.expression_df.shape[0] if result.expression_df is not None else "N/A")
        with col2:
            st.metric("Genes", result.expression_df.shape[1] if result.expression_df is not None else "N/A")
        with col3:
            st.metric("Data Type", result.data_type.value)
        
        # Data preview
        with st.expander("📋 Data Preview", expanded=True):
            st.dataframe(result.expression_df.head(10) if result.expression_df is not None else result.de_results_df.head(10))
        
        # QC plots
        if result.expression_df is not None:
            st.subheader("Quality Control")
            qc_col1, qc_col2 = st.columns(2)
            with qc_col1:
                st.plotly_chart(create_library_size_barplot(result.expression_df), use_container_width=True)
            with qc_col2:
                st.plotly_chart(create_count_distribution_boxplot(result.expression_df), use_container_width=True)
            
            qc_col3, qc_col4 = st.columns(2)
            with qc_col3:
                st.plotly_chart(create_gene_detection_plot(result.expression_df), use_container_width=True)
            with qc_col4:
                st.plotly_chart(create_sample_similarity_heatmap(result.expression_df), use_container_width=True)
        
        # Next button
        st.button("Next →", on_click=lambda: set_step(2), type="primary")

# STEP 2: Metadata
elif current_step == 2:
    st.header("Step 2: Assign Sample Metadata")
    
    result = st.session_state["parsed_result"]
    if not result or result.expression_df is None:
        st.warning("No expression data loaded. Please go back to Step 1.")
        st.button("← Back", on_click=lambda: set_step(1))
    else:
        # Show condition management
        col1, col2 = st.columns([2, 1])
        with col1:
            new_condition = st.text_input("Add new condition group", placeholder="e.g., Treatment, Control")
        with col2:
            st.markdown("")  # spacing
            if st.button("Add Condition") and new_condition:
                if new_condition not in st.session_state["conditions"]:
                    st.session_state["conditions"].append(new_condition)
        
        # Show available conditions as chips/tags
        if st.session_state["conditions"]:
            st.info(f"Available conditions: {', '.join(st.session_state['conditions'])}")
        
        # USE st.data_editor
        samples = result.expression_df.index.tolist()
        available_conditions = st.session_state["conditions"]
        
        if available_conditions:
            # Build data for editor
            editor_data = []
            for s in samples:
                current_cond = st.session_state["sample_metadata"].get(s, {}).get("condition", available_conditions[0])
                editor_data.append({"Sample": s, "Condition": current_cond})
            
            edited_df = st.data_editor(
                pd.DataFrame(editor_data),
                column_config={
                    "Sample": st.column_config.TextColumn("Sample Name", disabled=True),
                    "Condition": st.column_config.SelectboxColumn("Condition", options=available_conditions, required=True),
                },
                use_container_width=True,
                num_rows="fixed",
                hide_index=True,
            )
            
            # Update session state from editor
            for _, row in edited_df.iterrows():
                st.session_state["sample_metadata"][row["Sample"]] = {"condition": row["Condition"]}
        
        # Validation summary
        is_valid, errors = validate_metadata_assignment(st.session_state["sample_metadata"], result.expression_df)
        
        if not is_valid:
            for err in errors:
                st.error(err)
        else:
            st.success("Metadata valid! Ready to proceed.")
        
        # Navigation
        col1, col2 = st.columns(2)
        with col1:
            st.button("← Back", on_click=lambda: set_step(1))
        with col2:
            st.button("Next →", on_click=lambda: set_step(3), type="primary", disabled=not is_valid)

# STEP 3: Comparison & Analysis
elif current_step == 3:
    st.header("Step 3: Select Comparisons & Run Analysis")
    
    unique_conds = get_unique_conditions(st.session_state["sample_metadata"])
    
    st.subheader("Add Comparison")
    col1, col2, col3 = st.columns([2, 2, 1])
    with col1:
        test_cond = st.selectbox("🧪 Test Condition", unique_conds)
    with col2:
        ref_cond = st.selectbox("📋 Reference Condition", unique_conds)
    with col3:
        st.markdown("")  # spacing
        if st.button("➕ Add", use_container_width=True):
            valid, msg = validate_comparison((test_cond, ref_cond), unique_conds)
            if valid:
                add_comparison(test_cond, ref_cond)
            else:
                st.error(msg)
    
    # Show active comparisons as cards
    if st.session_state["comparisons"]:
        st.subheader("Active Comparisons")
        for i, (t, r) in enumerate(st.session_state["comparisons"]):
            col1, col2 = st.columns([5, 1])
            with col1:
                st.markdown(f"**{t}** vs **{r}**")
            with col2:
                if st.button("🗑️", key=f"del_{i}"):
                    remove_comparison((t, r))
                    st.rerun()
    else:
        st.warning("Please add at least one comparison.")
    
    st.markdown("---")
    
    # RUN ANALYSIS
    can_run = len(st.session_state["comparisons"]) > 0
    
    if st.button("🚀 Run Analysis", type="primary", use_container_width=True, disabled=not can_run):
        with st.status("Running analysis...", expanded=True) as status:
            try:
                result = st.session_state["parsed_result"]
                
                # 1. DE Analysis
                st.write("Running Differential Expression Analysis...")
                meta_df = pd.DataFrame.from_dict(st.session_state["sample_metadata"], orient="index")
                comparisons_tuple = tuple(tuple(c) for c in st.session_state["comparisons"])
                
                de_results = cached_de_analysis(
                    result.expression_df,
                    meta_df,
                    comparisons_tuple,
                    "condition",
                )
                st.session_state["de_results"] = de_results
                
                if st.session_state["comparisons"]:
                    st.session_state["active_comparison"] = st.session_state["comparisons"][0]
                
                # 2. Pathway Enrichment
                st.write("Running Pathway Enrichment...")
                pe = PathwayEnrichment()
                enrichment_results = {}
                
                for comp, de_res in st.session_state["de_results"].items():
                    if de_res.results_df is not None and not de_res.results_df.empty:
                        genes, error = pe.select_genes_for_enrichment(de_res.results_df)
                        res_obj = EnrichmentResult(
                            go_results=pd.DataFrame(),
                            kegg_results=pd.DataFrame(),
                            genes_used=genes,
                            selection_note=f"{len(genes)} genes selected",
                            error=error,
                        )
                        if not error:
                            genes_tuple = tuple(sorted(genes))
                            go_res, go_err = cached_enrichment(genes_tuple, "go_bp")
                            kegg_res, kegg_err = cached_enrichment(genes_tuple, "kegg")
                            res_obj.go_results = go_res
                            res_obj.kegg_results = kegg_res
                            if go_err or kegg_err:
                                res_obj.error = f"GO: {go_err}, KEGG: {kegg_err}"
                        enrichment_results[comp] = res_obj
                st.session_state["enrichment_results"] = enrichment_results
                
                # 2b. Additional Enrichment Databases
                st.write("Running Additional Enrichment Databases...")
                try:
                    for comp, enr_obj in enrichment_results.items():
                        if enr_obj.error and not enr_obj.genes_used:
                            continue
                        genes = enr_obj.genes_used
                        if not genes:
                            continue
                        # Additional databases
                        genes_tuple = tuple(sorted(genes))
                        for attr, db_name in [
                            ('go_mf_results', 'go_mf'),
                            ('go_cc_results', 'go_cc'),
                            ('reactome_results', 'reactome'),
                            ('wikipathway_results', 'wikipathway'),
                            ('hallmark_results', 'hallmark'),
                        ]:
                            try:
                                res_df, err = cached_enrichment(genes_tuple, db_name)
                                if not err:
                                    setattr(enr_obj, attr, res_df)
                            except Exception:
                                pass
                except Exception:
                    pass
                
                # 2c. GSEA Analysis
                st.write("Running Gene Set Enrichment Analysis (GSEA)...")
                gsea_results = {}
                try:
                    for comp, de_res in st.session_state["de_results"].items():
                        if de_res.results_df is not None and not de_res.results_df.empty:
                            try:
                                ranks = pe.create_ranking_metric(de_res.results_df)
                                if len(ranks) > 0:
                                    gsea_df, gsea_err = pe.run_gsea(ranks, database='MSigDB_Hallmark_2020')
                                    if not gsea_err:
                                        gsea_results[comp] = gsea_df
                                        if comp in enrichment_results:
                                            enrichment_results[comp].gsea_results = {'Hallmark': gsea_df}
                            except Exception:
                                pass
                except Exception:
                    pass
                st.session_state["gsea_results"] = gsea_results
                
                # 2d. Interpretation Engine
                st.write("Generating Interpretations...")
                all_interpretations = {}
                try:
                    interp_engine = InterpretationEngine()
                    for comp, de_res in st.session_state["de_results"].items():
                        if de_res.results_df is not None and not de_res.results_df.empty:
                            comp_name = f"{comp[0]} vs {comp[1]}"
                            comp_interps = interp_engine.interpret_de_results(
                                de_res.results_df, comp_name,
                                padj_threshold=padj_threshold,
                                lfc_threshold=lfc_threshold
                            )
                            # Enrichment interpretations
                            if comp in enrichment_results:
                                enr = enrichment_results[comp]
                                if not enr.go_results.empty:
                                    comp_interps.extend(interp_engine.interpret_enrichment(enr.go_results))
                            
                            # Top genes analysis
                            try:
                                comp_interps.extend(interp_engine.interpret_top_genes(de_res.results_df, comp_name))
                            except Exception:
                                pass
                            
                            # Volcano pattern analysis
                            try:
                                comp_interps.extend(interp_engine.interpret_volcano_pattern(de_res.results_df, comp_name))
                            except Exception:
                                pass
                            
                            # Misinterpretation warnings
                            try:
                                enr_df = enrichment_results[comp].go_results if comp in enrichment_results else None
                                comp_interps.extend(interp_engine.generate_misinterpretation_warnings(de_res.results_df, enr_df))
                            except Exception:
                                pass
                            
                            all_interpretations[comp] = comp_interps
                except Exception:
                    pass
                
                # Global interpretations (not per-comparison)
                global_interps = []
                
                # PCA interpretation
                try:
                    if st.session_state.get("pca_coordinates") is not None:
                        pca_interps = interp_engine.interpret_pca(
                            st.session_state["pca_coordinates"],
                            st.session_state.get("sample_metadata", {})
                        )
                        global_interps.extend(pca_interps)
                except Exception:
                    pass
                
                # GSEA interpretation
                try:
                    for comp, gsea_res in st.session_state.get("gsea_results", {}).items():
                        if gsea_res is not None and not gsea_res.empty:
                            gsea_interps = interp_engine.interpret_gsea_results(gsea_res)
                            comp_key = comp if comp in all_interpretations else list(all_interpretations.keys())[0] if all_interpretations else None
                            if comp_key:
                                all_interpretations[comp_key].extend(gsea_interps)
                except Exception:
                    pass
                
                # Deconvolution interpretation  
                try:
                    if st.session_state.get("deconvolution_results") is not None:
                        deconv_interps = interp_engine.interpret_deconvolution(
                            st.session_state["deconvolution_results"],
                            st.session_state.get("sample_metadata")
                        )
                        global_interps.extend(deconv_interps)
                except Exception:
                    pass
                
                # Add global interpretations to all comparisons
                for comp in all_interpretations:
                    all_interpretations[comp].extend(global_interps)
                
                st.session_state["interpretations"] = all_interpretations
                
                # 3. Visualizations
                st.write("Generating Visualizations...")
                figures = {}
                
                # Volcano (per comparison)
                for comp, de_res in st.session_state["de_results"].items():
                    if de_res.results_df is not None:
                        comp_name = f"{comp[0]}_vs_{comp[1]}"
                        figures[f"volcano_{comp_name}"] = create_volcano_plot(de_res.results_df)
                        figures[f"ma_{comp_name}"] = create_ma_plot(de_res.results_df)
                
                # Heatmap & PCA (Global)
                if st.session_state["de_results"]:
                    first_res = next(iter(st.session_state["de_results"].values()))
                    norm_counts = first_res.log_normalized_counts
                    
                    if norm_counts is not None:
                        sample_conds = {s: m["condition"] for s, m in st.session_state["sample_metadata"].items()}
                        figures["heatmap"] = create_clustered_heatmap(norm_counts.T, sample_conds)
                        figures["pca"] = create_pca_plot(norm_counts, sample_conds)
                        figures["correlation"] = create_correlation_heatmap(norm_counts.T, sample_conds)
                        
                        # Gene Panels
                        st.write("Analyzing Gene Panels...")
                        try:
                            analyzer = GenePanelAnalyzer(config_path="config/gene_panels.yaml")
                            for panel_name in analyzer.panels:
                                try:
                                    fig = analyzer.plot_panel(norm_counts, panel_name, sample_conds)
                                    figures[f"panel_{panel_name}"] = fig
                                except ValueError:
                                    pass
                        except FileNotFoundError:
                            pass
                        
                        # Cell Type Deconvolution
                        st.write("Running Cell Type Deconvolution...")
                        try:
                            deconv = CellTypeDeconvolution()
                            sig = deconv.get_default_skin_signature()
                            # norm_counts is samples × genes, need genes × samples
                            expr_for_deconv = norm_counts.T
                            deconv_result = deconv.run_nnls(expr_for_deconv, sig)
                            st.session_state["deconvolution_results"] = deconv_result
                        except Exception:
                            st.session_state["deconvolution_results"] = None
                        
                        # Advanced QC: Outlier Detection + PCA Coordinates
                        st.write("Running Advanced QC...")
                        try:
                            adv_qc = AdvancedQC()
                            # Compute PCA coordinates for advanced QC
                            pca_model = PCA(n_components=min(3, norm_counts.shape[0]))
                            pca_coords = pca_model.fit_transform(norm_counts.values)
                            pca_df = pd.DataFrame(
                                pca_coords[:, :min(3, pca_coords.shape[1])],
                                columns=[f'PC{i+1}' for i in range(min(3, pca_coords.shape[1]))],
                                index=norm_counts.index
                            )
                            st.session_state["pca_coordinates"] = pca_df
                            
                            outliers, distances = adv_qc.detect_outliers(pca_df)
                            st.session_state["outlier_info"] = {
                                'outliers': outliers,
                                'distances': distances,
                                'conditions': sample_conds
                            }
                        except Exception:
                            st.session_state["outlier_info"] = None
                
                st.session_state["figures"] = figures
                st.session_state["analysis_complete"] = True
                status.update(label="Analysis complete!", state="complete")
                
                st.session_state["current_step"] = 4
                st.rerun()
                
            except Exception as e:
                st.error(f"Analysis failed: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

    # Back button
    st.button("← Back", on_click=lambda: set_step(2))

# STEP 4: Results
elif current_step == 4:
    st.header("Step 4: Results")
    
    if not st.session_state["analysis_complete"]:
        st.warning("Analysis not complete. Please go back to Step 3.")
        st.button("← Back", on_click=lambda: set_step(3))
    else:
        # Comparison selector
        if len(st.session_state["comparisons"]) > 1:
            comps = st.session_state["comparisons"]
            comp_labels = [f"{t} vs {r}" for t, r in comps]
            selected_idx = st.selectbox(
                "Select Comparison to View",
                range(len(comps)),
                format_func=lambda x: comp_labels[x],
            )
            st.session_state["active_comparison"] = comps[selected_idx]
        
        active_comp = st.session_state["active_comparison"]
        res = st.session_state["de_results"].get(active_comp)
        
        if res:
            # Results tabs
            tab_summary, tab_de, tab_volcano, tab_ma, tab_heatmap, tab_pca, tab_enrichment, tab_gsea, tab_interpret, tab_deconv, tab_adv_qc, tab_panels, tab_gene_expr, tab_venn, tab_export = st.tabs([
                "📊 Summary", "📋 DE Results", "🌋 Volcano", "📈 MA Plot",
                "🗺️ Heatmap", "🔬 PCA", "🧬 Enrichment", "🔍 GSEA",
                "🧠 Interpretation", "🧫 Deconvolution", "⚠️ QC Advanced",
                "🎯 Gene Panels", "🧬 Gene Expression", "🔀 Venn Overlap", "📥 Export"
            ])
            
            with tab_summary:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["summary"])
                if res.results_df is not None:
                    summary = compute_de_summary(res.results_df)
                    col1, col2, col3, col4 = st.columns(4)
                    col1.metric("Total Genes", summary["total_genes"])
                    col2.metric("Significant", summary["significant_genes"])
                    col3.metric("Upregulated", summary["upregulated"], delta=f"↑")
                    col4.metric("Downregulated", summary["downregulated"], delta=f"↓", delta_color="inverse")
                    
                    col1, col2 = st.columns(2)
                    with col1:
                        st.subheader("Top Upregulated Genes")
                        st.dataframe(pd.DataFrame(summary["top_up_genes"], columns=["Gene", "Log2FC", "Padj"]))
                    with col2:
                        st.subheader("Top Downregulated Genes")
                        st.dataframe(pd.DataFrame(summary["top_down_genes"], columns=["Gene", "Log2FC", "Padj"]))
            
            with tab_de:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["de_results"])
                search_query = st.text_input("🔍 Search genes...", placeholder="Type gene name (e.g., COL1A1)")
                display_df = res.results_df
                if search_query:
                    display_df = search_genes(display_df, search_query)
                
                st.dataframe(display_df.sort_values("padj").head(1000))
                
                csv = display_df.to_csv().encode("utf-8")
                st.download_button(
                    "Download DE Results (CSV)",
                    csv,
                    f"de_results_{active_comp[0]}_vs_{active_comp[1]}.csv",
                    "text/csv",
                )
            
            with tab_volcano:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["volcano"])
                comp_name = f"{active_comp[0]}_vs_{active_comp[1]}"
                if f"volcano_{comp_name}" in st.session_state["figures"]:
                    plot_with_download(st.session_state["figures"][f"volcano_{comp_name}"], f"volcano_{comp_name}")
            
            with tab_ma:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["ma_plot"])
                comp_name = f"{active_comp[0]}_vs_{active_comp[1]}"
                if f"ma_{comp_name}" in st.session_state["figures"]:
                    plot_with_download(st.session_state["figures"][f"ma_{comp_name}"], f"ma_{comp_name}")
            
            with tab_heatmap:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["heatmap"])
                if "heatmap" in st.session_state["figures"]:
                    plot_with_download(st.session_state["figures"]["heatmap"], "heatmap")
                if "correlation" in st.session_state["figures"]:
                    plot_with_download(st.session_state["figures"]["correlation"], "correlation")
            
            with tab_pca:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["pca"])
                if "pca" in st.session_state["figures"]:
                    plot_with_download(st.session_state["figures"]["pca"], "pca")
            
            with tab_enrichment:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["enrichment"])
                if active_comp in st.session_state["enrichment_results"]:
                    enr = st.session_state["enrichment_results"][active_comp]
                    if enr.error:
                        st.warning(enr.error)
                    
                    # Dot plots
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("### GO Biological Process")
                        if not enr.go_results.empty:
                            plot_with_download(create_enrichment_dotplot(enr.go_results, title="GO Biological Process"), "go_bp")
                            with st.expander("View Table"):
                                st.dataframe(enr.go_results)
                        else:
                            st.info("No significant GO terms found.")
                    with col2:
                        st.write("### KEGG Pathways")
                        if not enr.kegg_results.empty:
                            plot_with_download(create_enrichment_dotplot(enr.kegg_results, title="KEGG Pathways"), "kegg")
                            with st.expander("View Table"):
                                st.dataframe(enr.kegg_results)
                        else:
                            st.info("No significant KEGG pathways found.")
                    
                    # Additional databases in expanders
                    additional_dbs = [
                        ('go_mf_results', 'GO Molecular Function'),
                        ('go_cc_results', 'GO Cellular Component'),
                        ('reactome_results', 'Reactome Pathways'),
                        ('wikipathway_results', 'WikiPathway'),
                        ('hallmark_results', 'MSigDB Hallmark'),
                    ]
                    for attr, label in additional_dbs:
                        db_df = getattr(enr, attr, pd.DataFrame())
                        if db_df is not None and not db_df.empty:
                            with st.expander(f"📂 {label} ({len(db_df)} terms)"):
                                plot_with_download(create_enrichment_dotplot(db_df, title=label), f"enrich_{attr}")
                                st.dataframe(db_df)
            
            with tab_gsea:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["gsea"])
                st.subheader("Gene Set Enrichment Analysis (GSEA)")
                if active_comp in st.session_state.get("gsea_results", {}):
                    gsea_df = st.session_state["gsea_results"][active_comp]
                    if gsea_df is not None and not gsea_df.empty:
                        st.markdown("**GSEA uses all genes ranked by expression change, not just significant ones.**")
                        
                        # Display NES bar chart
                        display_gsea = gsea_df.head(30).copy()
                        if 'NES' in display_gsea.columns:
                            term_col = 'Term' if 'Term' in display_gsea.columns else display_gsea.columns[0]
                            display_gsea = display_gsea.sort_values('NES')
                            
                            colors = ['#e74c3c' if v > 0 else '#3498db' for v in display_gsea['NES']]
                            fig_gsea = go.Figure(go.Bar(
                                x=display_gsea['NES'],
                                y=display_gsea[term_col].astype(str).apply(lambda x: x[:50] + '...' if len(x) > 50 else x),
                                orientation='h',
                                marker_color=colors,
                                hovertemplate='<b>%{y}</b><br>NES: %{x:.3f}<extra></extra>'
                            ))
                            fig_gsea.update_layout(
                                title='Normalized Enrichment Scores (NES)',
                                xaxis_title='NES',
                                yaxis_title='',
                                height=max(400, len(display_gsea) * 25 + 100),
                                margin=dict(l=300),
                            )
                            plot_with_download(fig_gsea, "gsea_nes")
                        
                        with st.expander("View Full GSEA Results Table"):
                            st.dataframe(gsea_df)
                    else:
                        st.info("GSEA analysis did not return results for this comparison.")
                else:
                    st.info("GSEA results not available. GSEA requires network access to MSigDB gene sets.")
            
            with tab_interpret:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["interpretation"])
                st.subheader("Automated Interpretations")
                comp_interps = st.session_state.get("interpretations", {}).get(active_comp, [])
                if comp_interps:
                    # Group by category
                    categories = {}
                    for interp in comp_interps:
                        cat = getattr(interp, 'category', 'general')
                        categories.setdefault(cat, []).append(interp)
                    
                    category_labels = {
                        "de": "🧬 Differential Expression Insights",
                        "enrichment": "🔬 Pathway & Enrichment Insights",
                        "qc": "📊 Quality Control Insights",
                        "deconvolution": "🧫 Cell Composition Insights",
                        "warning": "⚠️ Interpretation Cautions",
                    }
                    
                    for cat, interps in categories.items():
                        st.markdown(f"### {category_labels.get(cat, f'📝 {cat.title()} Insights')}")
                        for interp in interps:
                            icon = {"critical": "🔴", "important": "🟡", "informative": "🔵"}.get(interp.level, "⚪")
                            with st.expander(f"{icon} {interp.title}", expanded=(interp.level == "critical")):
                                st.markdown(interp.description)
                                if interp.recommendations:
                                    st.markdown("**Recommendations:**")
                                    for rec in interp.recommendations:
                                        st.markdown(f"- {rec}")
                                if interp.evidence:
                                    with st.expander("📋 Evidence Details"):
                                        st.json(interp.evidence)
                else:
                    st.info("Run analysis to generate interpretations.")
                
                # Methods text generator (keep existing)
                st.markdown("---")
                st.subheader("📝 Methods Text Generator")
                if st.button("Generate Methods Section"):
                    try:
                        interp_engine = InterpretationEngine()
                        params = {
                            'de_tool': 'PyDESeq2',
                            'design': '~ condition',
                            'padj': padj_threshold,
                            'lfc': lfc_threshold,
                            'databases': ['GO Biological Process', 'KEGG'],
                        }
                        methods_text = interp_engine.generate_methods_text(params)
                        st.markdown(methods_text)
                        st.code(methods_text, language='markdown')
                        st.download_button(
                            "Download Methods Text",
                            methods_text.encode('utf-8'),
                            "methods_section.md",
                            "text/markdown"
                        )
                    except Exception as e:
                        st.error(f"Failed to generate methods text: {str(e)}")
            
            with tab_deconv:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["deconvolution"])
                st.subheader("Cell Type Deconvolution")
                deconv_res = st.session_state.get("deconvolution_results")
                
                if deconv_res is not None:
                    st.markdown(f"**Method:** {deconv_res.method} | **Signature:** Default skin tissue (5 cell types)")
                    
                    # Stacked bar chart
                    sample_conds_for_deconv = {s: m["condition"] for s, m in st.session_state["sample_metadata"].items()}
                    plot_with_download(
                        create_deconvolution_plot(deconv_res.proportions, sample_conds_for_deconv),
                        "deconv"
                    )
                    
                    # Proportions table
                    with st.expander("View Proportions Table"):
                        st.dataframe(deconv_res.proportions.round(4))
                    
                    # Download
                    csv_deconv = deconv_res.proportions.to_csv().encode('utf-8')
                    st.download_button(
                        "Download Proportions (CSV)",
                        csv_deconv,
                        "cell_type_proportions.csv",
                        "text/csv"
                    )
                else:
                    st.info(
                        "Cell type deconvolution not available. This requires normalized expression data "
                        "and sufficient overlap with the skin tissue signature matrix."
                    )
            
            with tab_adv_qc:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["qc_advanced"])
                st.subheader("Advanced Quality Control")
                
                # Outlier Detection
                outlier_info = st.session_state.get("outlier_info")
                if outlier_info:
                    outliers = outlier_info['outliers']
                    distances = outlier_info['distances']
                    conds = outlier_info['conditions']
                    pca_df = st.session_state.get("pca_coordinates")
                    
                    if outliers:
                        st.warning(f"⚠️ Potential outlier samples detected: **{', '.join(outliers)}**")
                    else:
                        st.success("✅ No outlier samples detected (3 SD threshold)")
                    
                    if pca_df is not None:
                        adv_qc_viz = AdvancedQC()
                        fig_outlier = adv_qc_viz.create_outlier_plot(pca_df, conds, outliers, distances)
                        plot_with_download(fig_outlier, "qc_outlier")
                    
                    with st.expander("Sample Distances from Centroid"):
                        st.dataframe(distances.sort_values(ascending=False).to_frame('Distance').round(3))
                else:
                    st.info("Advanced QC requires PCA computation from normalized expression data.")
                
                # Batch Effect (if metadata has batch column)
                batch_info = st.session_state.get("batch_effects")
                if batch_info:
                    st.markdown("---")
                    st.subheader("Batch Effect Assessment")
                    adv_qc_viz = AdvancedQC()
                    fig_batch = adv_qc_viz.create_batch_effect_plot(batch_info)
                    plot_with_download(fig_batch, "qc_batch")
            
            with tab_panels:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["gene_panels"])
                panel_keys = [k for k in st.session_state["figures"].keys() if k.startswith("panel_")]
                if panel_keys:
                    for key in panel_keys:
                        plot_with_download(st.session_state["figures"][key], f"panel_{key}")
                else:
                    st.info("No gene panels available.")
            
            with tab_gene_expr:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["gene_expression"])
                st.subheader("Gene Expression Viewer")
                if res.results_df is not None:
                    all_genes = sorted(res.results_df["gene"].dropna().unique().tolist())
                    selected_gene = st.selectbox("Select Gene", all_genes, help="Type to search", key="gene_expr_select")
                    plot_type = st.radio("Plot Type", ["Box", "Violin"], horizontal=True, key="gene_expr_type")
                    norm_counts = res.log_normalized_counts
                    if norm_counts is not None and selected_gene in norm_counts.columns:
                        sample_conds = {s: m["condition"] for s, m in st.session_state["sample_metadata"].items()}
                        fig = create_gene_expression_plot(norm_counts, selected_gene, sample_conds, plot_type.lower())
                        plot_with_download(fig, f"gene_expr_{selected_gene}")
                    else:
                        st.warning(f"Gene '{selected_gene}' not found in normalized expression data.")
                else:
                    st.info("Run analysis to view gene expression.")

            with tab_venn:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["venn"])
                st.subheader("DE Gene Overlap")
                if len(st.session_state["de_results"]) >= 2:
                    gene_sets = {}
                    for comp, de_res in st.session_state["de_results"].items():
                        if de_res.results_df is not None:
                            sig = de_res.results_df[de_res.results_df["padj"] < padj_threshold]
                            sig_lfc = sig[sig["log2FoldChange"].abs() > lfc_threshold]
                            comp_label = f"{comp[0]} vs {comp[1]}"
                            gene_sets[comp_label] = set(sig_lfc["gene"].tolist())
                    if len(gene_sets) >= 2:
                        fig = create_venn_diagram(gene_sets)
                        plot_with_download(fig, "venn_overlap")
                        all_sets = list(gene_sets.values())
                        shared = set.intersection(*all_sets)
                        if shared:
                            st.markdown(f"**{len(shared)} genes shared across all comparisons:**")
                            st.dataframe(pd.DataFrame(sorted(shared), columns=["Gene"]))
                    else:
                        st.info("Need at least 2 comparisons with significant DE genes.")
                else:
                    st.info("Run multiple comparisons to see DE gene overlap.")

            with tab_export:
                with st.expander("ℹ️ How to read this tab", expanded=False):
                    st.markdown(HELP_TEXTS["export"])
                st.header("Export Results")
                col1, col2 = st.columns(2)
                
                sample_conds = {s: m["condition"] for s, m in st.session_state["sample_metadata"].items()}
                export_data = ExportData(
                    data_type=st.session_state["parsed_result"].data_type,
                    de_results=st.session_state["de_results"],
                    expression_matrix=st.session_state["expression_matrix"],
                    enrichment_results=st.session_state["enrichment_results"],
                    figures=st.session_state["figures"],
                    settings={
                        "padj_threshold": padj_threshold,
                        "lfc_threshold": lfc_threshold,
                        "comparisons": st.session_state["comparisons"],
                    },
                    sample_conditions=sample_conds,
                    active_comparison=st.session_state["active_comparison"],
                )
                engine = ExportEngine()
                
                with col1:
                    if st.button("Generate Excel Report"):
                        with st.spinner("Generating Excel..."):
                            tmp_path = None
                            try:
                                with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
                                    engine.export_excel(tmp.name, export_data)
                                    tmp_path = tmp.name
                                with open(tmp_path, "rb") as f:
                                    st.download_button("Download Excel Report", f.read(), "rnaseq_results.xlsx", "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                            finally:
                                if tmp_path and os.path.exists(tmp_path):
                                    os.unlink(tmp_path)
                
                with col2:
                    if st.button("Generate PDF Report"):
                        with st.spinner("Generating PDF..."):
                            tmp_path = None
                            try:
                                with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp:
                                    engine.export_pdf_report(tmp.name, export_data)
                                    tmp_path = tmp.name
                                with open(tmp_path, "rb") as f:
                                    st.download_button("Download PDF Report", f.read(), "rnaseq_report.pdf", "application/pdf")
                            finally:
                                if tmp_path and os.path.exists(tmp_path):
                                    os.unlink(tmp_path)
        
        st.button("← Back to Comparisons", on_click=lambda: set_step(3))
