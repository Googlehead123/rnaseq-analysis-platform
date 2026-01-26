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

# Import backend modules
from rnaseq_parser import RNASeqParser, DataType, ParseResult, ParserValidationError
from de_analysis import DEAnalysisEngine, DEResult
from pathway_enrichment import PathwayEnrichment
from visualizations import (
    create_volcano_plot,
    create_clustered_heatmap,
    create_pca_plot,
)
from gene_panels import GenePanelAnalyzer
from export_engine import ExportEngine, ExportData, EnrichmentResult

# --- Page Config ---
st.set_page_config(
    page_title="RNA-seq Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

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

# --- Constants ---
ALLOWED_GROUPS = ["control", "positive_control", "treatment"]

# --- Helper Functions ---


def get_unique_conditions(sample_metadata: Dict[str, Dict]) -> List[str]:
    """Get list of unique conditions from metadata."""
    if not sample_metadata:
        return []
    conditions = set(m["condition"] for m in sample_metadata.values())
    return sorted(list(conditions))


def count_samples_per_condition(sample_metadata: Dict[str, Dict]) -> Dict[str, int]:
    """Count how many samples are assigned to each condition."""
    counts = {}
    for meta in sample_metadata.values():
        condition = meta["condition"]
        counts[condition] = counts.get(condition, 0) + 1
    return counts


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


# --- Main App Layout ---

st.title("🧬 RNA-seq Analysis Platform")
st.markdown("---")

# Sidebar - Setup
with st.sidebar:
    st.header("1. Data Upload")
    uploaded_file = st.file_uploader(
        "Upload Counts File (CSV/TSV/Excel)", type=["csv", "tsv", "txt", "xlsx"]
    )

    if uploaded_file:
        try:
            # Save uploaded file to temp file for parser
            import tempfile

            suffix = f".{uploaded_file.name.split('.')[-1]}"
            with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
                tmp.write(uploaded_file.getvalue())
                tmp_path = tmp.name

            # Parse file
            parser = RNASeqParser()
            result = parser.parse(tmp_path)

            # Update session state if new file
            is_new_file = False
            if st.session_state["parsed_result"] is None:
                is_new_file = True
            elif hasattr(st.session_state["parsed_result"], "expression_df"):
                old_df = st.session_state["parsed_result"].expression_df
                new_df = result.expression_df
                # Handle None cases safely
                if old_df is None and new_df is None:
                    is_new_file = False  # Both None (e.g. pre-analyzed)
                    # Check if de_results_df changed for pre-analyzed
                    old_de = getattr(
                        st.session_state["parsed_result"], "de_results_df", None
                    )
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
                st.session_state["analysis_complete"] = False  # Reset analysis
                st.session_state["sample_metadata"] = {}  # Reset metadata
                st.session_state["comparisons"] = []  # Reset comparisons

                st.success(f"Loaded {result.data_type.value} data successfully!")
                if result.warnings:
                    for w in result.warnings:
                        st.warning(w)

            # Clean up temp file
            os.unlink(tmp_path)

        except ParserValidationError as e:
            st.error(f"Error parsing file: {e.message}")
            if e.details:
                st.json(e.details)
        except Exception as e:
            st.error(f"Unexpected error: {str(e)}")

    # Only show setup if data is loaded
    if st.session_state["parsed_result"]:
        result = st.session_state["parsed_result"]

        # Data Preview
        with st.expander("Data Preview"):
            if result.expression_df is not None:
                st.dataframe(result.expression_df.head())
                st.caption(
                    f"Shape: {result.expression_df.shape[0]} samples × {result.expression_df.shape[1]} genes"
                )
            elif result.de_results_df is not None:
                st.dataframe(result.de_results_df.head())

        # Metadata Assignment (Only for RAW_COUNTS)
        if result.data_type == DataType.RAW_COUNTS:
            st.header("2. Metadata")

            # Condition management
            new_condition = st.text_input("Add Condition Group")
            if st.button("Add Group") and new_condition:
                if new_condition not in st.session_state["conditions"]:
                    st.session_state["conditions"].append(new_condition)

            available_conditions = st.session_state["conditions"]

            # Sample assignment
            st.subheader("Assign Samples")
            samples = result.expression_df.index.tolist()

            # Initialize metadata if empty
            if not st.session_state["sample_metadata"]:
                for s in samples:
                    st.session_state["sample_metadata"][s] = {
                        "condition": available_conditions[0]
                    }

            # Bulk assignment helper
            col1, col2 = st.columns(2)
            with col1:
                start_idx = st.number_input(
                    "Start Index", min_value=1, max_value=len(samples), value=1
                )
            with col2:
                end_idx = st.number_input(
                    "End Index",
                    min_value=1,
                    max_value=len(samples),
                    value=min(len(samples), 5),
                )

            assign_cond = st.selectbox("Assign Range To", available_conditions)
            if st.button("Assign Range"):
                for i in range(start_idx - 1, end_idx):
                    st.session_state["sample_metadata"][samples[i]]["condition"] = (
                        assign_cond
                    )
                st.rerun()

            # Individual assignment
            with st.expander("Individual Sample Assignment"):
                for sample in samples:
                    current_cond = (
                        st.session_state["sample_metadata"]
                        .get(sample, {})
                        .get("condition", available_conditions[0])
                    )
                    new_cond = st.selectbox(
                        f"Condition for {sample}",
                        available_conditions,
                        index=available_conditions.index(current_cond)
                        if current_cond in available_conditions
                        else 0,
                        key=f"sel_{sample}",
                    )
                    st.session_state["sample_metadata"][sample] = {
                        "condition": new_cond
                    }

            # Validation
            is_valid, errors = validate_metadata_assignment(
                st.session_state["sample_metadata"], result.expression_df
            )

            if not is_valid:
                for err in errors:
                    st.error(err)
            else:
                st.success("Metadata valid!")

            # Comparisons
            st.header("3. Comparisons")
            unique_conds = get_unique_conditions(st.session_state["sample_metadata"])

            c1, c2 = st.columns(2)
            with c1:
                test_cond = st.selectbox(
                    "Test Condition", unique_conds, key="test_cond"
                )
            with c2:
                ref_cond = st.selectbox(
                    "Reference Condition", unique_conds, key="ref_cond"
                )

            if st.button("Add Comparison"):
                valid, msg = validate_comparison((test_cond, ref_cond), unique_conds)
                if valid:
                    add_comparison(test_cond, ref_cond)
                else:
                    st.error(msg)

            # List comparisons
            if st.session_state["comparisons"]:
                st.write("Active Comparisons:")
                for i, (t, r) in enumerate(st.session_state["comparisons"]):
                    col_a, col_b = st.columns([4, 1])
                    with col_a:
                        st.write(f"**{t}** vs **{r}**")
                    with col_b:
                        if st.button("🗑️", key=f"del_{i}"):
                            remove_comparison((t, r))
                            st.rerun()
            else:
                st.warning("Please add at least one comparison.")

        # Analysis Runner
        st.header("4. Analysis")

        can_run = False
        # Initialize is_valid to False by default
        is_valid = False

        if result.data_type == DataType.RAW_COUNTS:
            # Re-validate to ensure is_valid is set correctly in this scope
            is_valid, _ = validate_metadata_assignment(
                st.session_state["sample_metadata"], result.expression_df
            )
            can_run = is_valid and len(st.session_state["comparisons"]) > 0
        elif result.data_type == DataType.PRE_ANALYZED:
            can_run = True
        elif result.data_type == DataType.NORMALIZED:
            can_run = True  # Can run viz only

        if st.button("Run Analysis", disabled=not can_run):
            progress_bar = st.progress(0)
            status_text = st.empty()

            try:
                # 1. DE Analysis
                status_text.text("Running Differential Expression Analysis...")
                progress_bar.progress(20)

                de_engine = DEAnalysisEngine()

                if result.data_type == DataType.RAW_COUNTS:
                    # Convert metadata to DataFrame
                    meta_df = pd.DataFrame.from_dict(
                        st.session_state["sample_metadata"], orient="index"
                    )

                    de_results = de_engine.run_all_comparisons(
                        counts_df=result.expression_df,
                        metadata_df=meta_df,
                        comparisons=st.session_state["comparisons"],
                        design_factor="condition",
                    )
                    st.session_state["de_results"] = de_results

                    # Set active comparison to first one
                    if st.session_state["comparisons"]:
                        st.session_state["active_comparison"] = st.session_state[
                            "comparisons"
                        ][0]

                elif result.data_type == DataType.PRE_ANALYZED:
                    # Wrap pre-analyzed result in DEResult
                    # Note: Pre-analyzed doesn't have normalized counts usually, unless provided separately
                    # For now, we just store the results_df
                    dummy_comp = ("Uploaded", "File")
                    st.session_state["de_results"] = {
                        dummy_comp: DEResult(
                            results_df=result.de_results_df,
                            normalized_counts=None,
                            log_normalized_counts=None,
                            dds=None,
                            comparison=dummy_comp,
                            n_significant=len(
                                result.de_results_df[
                                    result.de_results_df["padj"] < 0.05
                                ]
                            )
                            if "padj" in result.de_results_df.columns
                            else 0,
                            warnings=[],
                        )
                    }
                    st.session_state["active_comparison"] = dummy_comp
                    st.session_state["comparisons"] = [dummy_comp]

                elif result.data_type == DataType.NORMALIZED:
                    # No DE, just viz
                    st.session_state["de_results"] = {}
                    st.session_state["active_comparison"] = None

                # 2. Pathway Enrichment
                status_text.text("Running Pathway Enrichment...")
                progress_bar.progress(50)

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
                            go_res, go_err = pe.get_go_enrichment(genes)
                            kegg_res, kegg_err = pe.get_kegg_enrichment(genes)

                            res_obj.go_results = go_res
                            res_obj.kegg_results = kegg_res
                            if go_err or kegg_err:
                                res_obj.error = f"GO: {go_err}, KEGG: {kegg_err}"

                        enrichment_results[comp] = res_obj

                st.session_state["enrichment_results"] = enrichment_results

                # 3. Visualizations
                status_text.text("Generating Visualizations...")
                progress_bar.progress(70)

                figures = {}

                # Volcano (per comparison)
                for comp, de_res in st.session_state["de_results"].items():
                    if de_res.results_df is not None:
                        comp_name = f"{comp[0]}_vs_{comp[1]}"
                        figures[f"volcano_{comp_name}"] = create_volcano_plot(
                            de_res.results_df
                        )

                # Heatmap & PCA (Global)
                # Need normalized counts.
                # For RAW_COUNTS, get from first DE result (all share same normalized counts)
                # For NORMALIZED, use expression_df directly
                norm_counts = None
                if (
                    result.data_type == DataType.RAW_COUNTS
                    and st.session_state["de_results"]
                ):
                    first_res = next(iter(st.session_state["de_results"].values()))
                    norm_counts = first_res.log_normalized_counts
                elif result.data_type == DataType.NORMALIZED:
                    import numpy as np

                    norm_counts = np.log2(result.expression_df + 1)

                if norm_counts is not None:
                    # Prepare sample conditions map
                    if result.data_type == DataType.RAW_COUNTS:
                        sample_conds = {
                            s: m["condition"]
                            for s, m in st.session_state["sample_metadata"].items()
                        }
                    else:
                        # For normalized, assume all same condition or infer from columns?
                        # Plan says "NORMALIZED input (floats) should disable DE analysis but enable viz."
                        # We'll just use dummy conditions if not provided
                        sample_conds = {s: "Sample" for s in norm_counts.index}

                    # Heatmap: genes x samples (TRANSPOSE REQUIRED)
                    figures["heatmap"] = create_clustered_heatmap(
                        norm_counts.T, sample_conds
                    )

                    # PCA: samples x genes (NO TRANSPOSE)
                    figures["pca"] = create_pca_plot(norm_counts, sample_conds)

                    # 4. Gene Panels
                    status_text.text("Analyzing Gene Panels...")
                    progress_bar.progress(90)

                    try:
                        analyzer = GenePanelAnalyzer(
                            config_path="config/gene_panels.yaml"
                        )
                        panels = analyzer.panels

                        for panel_name in panels:
                            try:
                                fig = analyzer.plot_panel(
                                    norm_counts, panel_name, sample_conds
                                )
                                figures[f"panel_{panel_name}"] = fig
                            except ValueError:
                                pass  # Skip panels with insufficient genes
                    except FileNotFoundError:
                        st.warning("Gene panels config not found.")

                st.session_state["figures"] = figures
                st.session_state["analysis_complete"] = True

                progress_bar.progress(100)
                status_text.text("Analysis Complete!")
                st.rerun()

            except Exception as e:
                st.error(f"Analysis failed: {str(e)}")
                import traceback

                st.code(traceback.format_exc())

# Main Content - Results
if st.session_state["analysis_complete"]:
    # Comparison Selector (if multiple)
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

    # Tabs
    tab1, tab2, tab3, tab4 = st.tabs(
        [
            "Differential Expression",
            "Pathway Enrichment",
            "Visualizations",
            "Gene Panels",
        ]
    )

    with tab1:
        st.subheader("Differential Expression Results")
        if active_comp and active_comp in st.session_state["de_results"]:
            res = st.session_state["de_results"][active_comp]
            if res.results_df is not None:
                st.dataframe(res.results_df.sort_values("padj").head(1000))

                # Download CSV
                csv = res.results_df.to_csv().encode("utf-8")
                st.download_button(
                    "Download DE Results (CSV)",
                    csv,
                    f"de_results_{active_comp[0]}_vs_{active_comp[1]}.csv",
                    "text/csv",
                )
            else:
                st.warning("No results for this comparison.")
        else:
            st.info("No active comparison selected.")

    with tab2:
        st.subheader("Pathway Enrichment")
        if active_comp and active_comp in st.session_state["enrichment_results"]:
            enr = st.session_state["enrichment_results"][active_comp]

            if enr.error:
                st.warning(enr.error)

            col1, col2 = st.columns(2)
            with col1:
                st.write("### GO Biological Process")
                if not enr.go_results.empty:
                    st.dataframe(enr.go_results)
                else:
                    st.info("No significant GO terms found.")

            with col2:
                st.write("### KEGG Pathways")
                if not enr.kegg_results.empty:
                    st.dataframe(enr.kegg_results)
                else:
                    st.info("No significant KEGG pathways found.")
        else:
            st.info("No enrichment results available.")

    with tab3:
        st.subheader("Visualizations")

        # Volcano
        if active_comp:
            comp_name = f"{active_comp[0]}_vs_{active_comp[1]}"
            volcano_key = f"volcano_{comp_name}"
            if volcano_key in st.session_state["figures"]:
                st.plotly_chart(
                    st.session_state["figures"][volcano_key], use_container_width=True
                )

        col1, col2 = st.columns(2)
        with col1:
            if "heatmap" in st.session_state["figures"]:
                st.plotly_chart(
                    st.session_state["figures"]["heatmap"], use_container_width=True
                )
        with col2:
            if "pca" in st.session_state["figures"]:
                st.plotly_chart(
                    st.session_state["figures"]["pca"], use_container_width=True
                )

    with tab4:
        st.subheader("Gene Panels")
        # Display all panel figures
        panel_keys = [
            k for k in st.session_state["figures"].keys() if k.startswith("panel_")
        ]
        if panel_keys:
            for key in panel_keys:
                st.plotly_chart(
                    st.session_state["figures"][key], use_container_width=True
                )
        else:
            st.info(
                "No gene panels available (check if config exists and genes are present in data)."
            )

    # Export Section
    st.markdown("---")
    st.header("5. Export Results")

    col1, col2 = st.columns(2)

    # Prepare ExportData
    # Need to reconstruct sample_conditions for export
    sample_conds = {}
    if st.session_state["sample_metadata"]:
        sample_conds = {
            s: m["condition"] for s, m in st.session_state["sample_metadata"].items()
        }

    export_data = ExportData(
        data_type=st.session_state["parsed_result"].data_type,
        de_results=st.session_state["de_results"],
        expression_matrix=st.session_state[
            "expression_matrix"
        ],  # This is raw counts for RAW_COUNTS type
        enrichment_results=st.session_state["enrichment_results"],
        figures=st.session_state["figures"],
        settings={
            "padj_threshold": 0.05,
            "lfc_threshold": 1.0,
            "comparisons": st.session_state["comparisons"],
        },
        sample_conditions=sample_conds,
        active_comparison=st.session_state["active_comparison"],
    )

    engine = ExportEngine()

    with col1:
        if st.button("Generate Excel Report"):
            with st.spinner("Generating Excel..."):
                buffer = io.BytesIO()
                # We need to patch export_excel to accept a buffer or file-like object
                # The current implementation takes a filepath.
                # Let's use a temp file.
                import tempfile

                with tempfile.NamedTemporaryFile(delete=False, suffix=".xlsx") as tmp:
                    engine.export_excel(tmp.name, export_data)
                    tmp_path = tmp.name

                with open(tmp_path, "rb") as f:
                    st.download_button(
                        "Download Excel Report",
                        f.read(),
                        "rnaseq_results.xlsx",
                        "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )
                os.unlink(tmp_path)

    with col2:
        if st.button("Generate PDF Report"):
            with st.spinner("Generating PDF..."):
                import tempfile

                with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp:
                    engine.export_pdf_report(tmp.name, export_data)
                    tmp_path = tmp.name

                with open(tmp_path, "rb") as f:
                    st.download_button(
                        "Download PDF Report",
                        f.read(),
                        "rnaseq_report.pdf",
                        "application/pdf",
                    )
                os.unlink(tmp_path)

else:
    if not st.session_state["parsed_result"]:
        st.info("👋 Welcome! Please upload a counts file in the sidebar to begin.")
    else:
        st.info("Please configure metadata and comparisons, then click 'Run Analysis'.")
