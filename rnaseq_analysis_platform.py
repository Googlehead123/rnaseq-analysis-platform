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
if "current_step" not in st.session_state:
    st.session_state["current_step"] = 1

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

# --- Sidebar ---
with st.sidebar:
    st.title("🧬 RNA-seq Analysis")
    st.markdown("---")
    
    # Step indicator
    current_step = st.session_state["current_step"]
    steps = ["📤 Upload", "🏷️ Metadata", "🔬 Compare", "📊 Results"]
    for i, step_name in enumerate(steps, 1):
        if i < current_step:
            st.markdown(f"✅ {step_name}")
        elif i == current_step:
            st.markdown(f"🔵 **{step_name}**")
        else:
            st.markdown(f"⚪ {step_name}")
    
    st.markdown("---")
    # Global settings only
    st.subheader("Settings")
    padj_threshold = st.slider("P-adj threshold", 0.001, 0.1, 0.05, key="padj_thresh")
    lfc_threshold = st.slider("Log2FC threshold", 0.5, 3.0, 1.0, key="lfc_thresh")

# --- Main Area ---

# STEP 1: Upload & Preview
if current_step == 1:
    st.header("Step 1: Upload Data")
    st.markdown("Upload your RNA-seq count matrix or load demo data to explore the platform.")
    
    col1, col2 = st.columns([3, 1])
    with col1:
        uploaded_file = st.file_uploader("Upload Counts File (CSV/TSV/Excel)", type=["csv", "tsv", "txt", "xlsx"])
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
    if uploaded_file:
        tmp_path = None
        try:
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
                st.success(f"Loaded {result.data_type.value} data successfully!")
                
                # Auto-detect metadata logic (copied from original)
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
                de_engine = DEAnalysisEngine()
                
                meta_df = pd.DataFrame.from_dict(st.session_state["sample_metadata"], orient="index")
                
                de_results = de_engine.run_all_comparisons(
                    counts_df=result.expression_df,
                    metadata_df=meta_df,
                    comparisons=st.session_state["comparisons"],
                    design_factor="condition",
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
                            go_res, go_err = pe.get_go_enrichment(genes)
                            kegg_res, kegg_err = pe.get_kegg_enrichment(genes)
                            res_obj.go_results = go_res
                            res_obj.kegg_results = kegg_res
                            if go_err or kegg_err:
                                res_obj.error = f"GO: {go_err}, KEGG: {kegg_err}"
                        enrichment_results[comp] = res_obj
                st.session_state["enrichment_results"] = enrichment_results
                
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
            tab_summary, tab_de, tab_volcano, tab_ma, tab_heatmap, tab_pca, tab_enrichment, tab_panels, tab_export = st.tabs([
                "📊 Summary", "📋 DE Results", "🌋 Volcano", "📈 MA Plot", 
                "🗺️ Heatmap", "🔬 PCA", "🧬 Enrichment", "🎯 Gene Panels", "📥 Export"
            ])
            
            with tab_summary:
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
                comp_name = f"{active_comp[0]}_vs_{active_comp[1]}"
                if f"volcano_{comp_name}" in st.session_state["figures"]:
                    st.plotly_chart(st.session_state["figures"][f"volcano_{comp_name}"], use_container_width=True)
            
            with tab_ma:
                comp_name = f"{active_comp[0]}_vs_{active_comp[1]}"
                if f"ma_{comp_name}" in st.session_state["figures"]:
                    st.plotly_chart(st.session_state["figures"][f"ma_{comp_name}"], use_container_width=True)
            
            with tab_heatmap:
                if "heatmap" in st.session_state["figures"]:
                    st.plotly_chart(st.session_state["figures"]["heatmap"], use_container_width=True)
                if "correlation" in st.session_state["figures"]:
                    st.plotly_chart(st.session_state["figures"]["correlation"], use_container_width=True)
            
            with tab_pca:
                if "pca" in st.session_state["figures"]:
                    st.plotly_chart(st.session_state["figures"]["pca"], use_container_width=True)
            
            with tab_enrichment:
                if active_comp in st.session_state["enrichment_results"]:
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
            
            with tab_panels:
                panel_keys = [k for k in st.session_state["figures"].keys() if k.startswith("panel_")]
                if panel_keys:
                    for key in panel_keys:
                        st.plotly_chart(st.session_state["figures"][key], use_container_width=True)
                else:
                    st.info("No gene panels available.")
            
            with tab_export:
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
