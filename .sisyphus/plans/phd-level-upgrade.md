# PhD-Level Platform Upgrade

## TL;DR

> **Quick Summary**: Upgrade the RNA-seq Streamlit platform to PhD/expert level by adding professional CSS theming, per-plot downloads, gene expression viewer, caching, inline help, Venn diagrams, and session save/load — closing gaps vs iDEP/Degust/NASQAR.
> 
> **Deliverables**:
> - Professional CSS theme (genomics blues/teals)
> - Per-plot PNG/SVG download buttons on every chart
> - Gene Expression box/violin plot tab
> - `@st.cache_data` on expensive operations
> - Inline help tooltips on all major widgets/tabs
> - Venn/UpSet diagram for multi-comparison DE overlap
> - Session save/load (JSON export/import)
> - (Bonus) Normalization comparison view, styled step indicator
> 
> **Estimated Effort**: Large
> **Parallel Execution**: YES - 3 waves
> **Critical Path**: Task 1 (CSS) and Task 3 (visualizations.py new functions) → Task 6 (main app integration of new tabs) → Task 7 (caching)

---

## Context

### Original Request
Upgrade existing RNA-seq platform to PhD/expert level based on competitive analysis of 8 top platforms. Focus on 7 HIGH PRIORITY gaps and 3 MEDIUM PRIORITY items.

### Architecture Constraints
- Single main file: `rnaseq_analysis_platform.py` (962 lines)
- ALL Plotly, NO matplotlib
- NO custom CSS currently — default Streamlit gray
- Must NOT break 147 passing tests
- Backend modules are stable — don't modify
- New file allowed: `session_manager.py`

---

## Work Objectives

### Core Objective
Close the top 7 feature gaps between our platform and professional tools (iDEP, Degust, NASQAR, BioJupies).

### Definition of Done
- [ ] Platform has professional CSS theme with branded colors
- [ ] Every Plotly chart has PNG/SVG download button
- [ ] Gene expression box/violin plot tab works
- [ ] Expensive ops are cached with `@st.cache_data`
- [ ] Key widgets have `help=` tooltips
- [ ] Venn diagram shows DE overlap for multi-comparison
- [ ] Session state can be saved/loaded as JSON
- [ ] All 147 existing tests still pass

### Must NOT Have (Guardrails)
- NO matplotlib — all Plotly
- NO modifications to backend modules (de_analysis.py, pathway_enrichment.py, etc.)
- NO modifications to test files
- NO breaking changes to existing tab content
- NO external CSS files — all injected via `st.markdown(unsafe_allow_html=True)`
- NO new pip dependencies beyond what's installed

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest, 150 tests)
- **User wants tests**: Tests-after (don't break existing; new features verified manually via Streamlit)
- **Framework**: pytest

### Verification Approach
- `pytest tests/ -v` must show 147+ passes (same as baseline)
- New `visualizations.py` functions get unit tests in existing test files
- UI features verified via `streamlit run` + browser inspection

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately — all independent):
├── Task 1: CSS Theme (rnaseq_analysis_platform.py only)
├── Task 2: visualizations.py — add 3 new functions
├── Task 3: session_manager.py — new file
└── Task 4: Inline help text constants

Wave 2 (After Wave 1):
├── Task 5: Per-plot download buttons (needs existing figures)
├── Task 6: New tabs — Gene Expression + Venn (needs Task 2 functions)
└── Task 7: Session save/load UI (needs Task 3)

Wave 3 (After Wave 2):
├── Task 8: Caching (@st.cache_data — touches analysis pipeline)
├── Task 9: Step indicator styling (needs Task 1 CSS)
└── Task 10: Verification & test pass confirmation
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 5, 9 | 2, 3, 4 |
| 2 | None | 6 | 1, 3, 4 |
| 3 | None | 7 | 1, 2, 4 |
| 4 | None | 5, 6, 7 | 1, 2, 3 |
| 5 | 1 | 10 | 6, 7 |
| 6 | 2 | 10 | 5, 7 |
| 7 | 3 | 10 | 5, 6 |
| 8 | 6 | 10 | 9 |
| 9 | 1 | 10 | 8 |
| 10 | 5-9 | None | None (final) |

---

## TODOs

- [ ] 1. Professional CSS Theme Injection

  **What to do**:
  - Add a CSS string constant at top of `rnaseq_analysis_platform.py` (after imports, before session state init)
  - Inject via `st.markdown(f"<style>{CSS}</style>", unsafe_allow_html=True)` right after `st.set_page_config()`
  - CSS must cover:
    - Sidebar: dark genomics-blue background (#1a1a2e or #0f3460), white text, styled nav
    - Metric cards: bordered cards with subtle shadow, colored delta indicators
    - Step indicator: pill-shaped steps with active/completed states
    - Tab headers: teal accent color (#16a085 or #1abc9c) for active tab
    - Buttons: primary button teal, hover state
    - Headers: professional font sizing, slight letter-spacing
    - Expander headers: subtle background tint
    - DataFrames: alternating row colors
  - Color palette: primary=#0f3460 (dark blue), accent=#16a085 (teal), bg=#f8f9fa (light gray), text=#2c3e50

  **Must NOT do**:
  - Don't use external CSS files
  - Don't override Streamlit's core layout (keep `layout="wide"`)
  - Don't add custom fonts via @import (slow loading)

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]
    - `frontend-ui-ux`: CSS theming is pure UI/UX work

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2, 3, 4)
  - **Blocks**: Tasks 5, 9
  - **Blocked By**: None

  **References**:
  - `rnaseq_analysis_platform.py:44-49` — `st.set_page_config()` block, inject CSS right after this
  - `rnaseq_analysis_platform.py:167-182` — Sidebar section with step indicator (will be styled by CSS)
  - `rnaseq_analysis_platform.py:661-668` — Metric cards in Summary tab (style these)

  **Acceptance Criteria**:
  ```bash
  # Verify CSS is injected (string present in file)
  grep -c "unsafe_allow_html=True" rnaseq_analysis_platform.py
  # Should be >= 1 (was 0 before)
  
  grep -c "<style>" rnaseq_analysis_platform.py
  # Should be >= 1
  
  # Verify no matplotlib imports crept in
  grep -c "matplotlib" rnaseq_analysis_platform.py
  # Should be 0
  
  # Tests still pass
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES
  - Message: `feat(ui): add professional CSS theme with genomics color palette`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 2. New Visualization Functions in visualizations.py

  **What to do**:
  - Add 3 new functions to `visualizations.py`:
  
  **2a. `create_gene_expression_plot()`**:
  ```python
  def create_gene_expression_plot(
      expression_df: pd.DataFrame,  # samples × genes
      gene: str,
      sample_conditions: Dict[str, str],
      plot_type: str = "box",  # "box" or "violin"
  ) -> go.Figure:
  ```
  - Extract the gene's expression across all samples
  - Group by condition from `sample_conditions`
  - Box plot: `go.Box()` with individual points overlaid (`boxpoints='all'`)
  - Violin plot: `go.Violin()` with points
  - Color by condition (use Plotly qualitative palette)
  - Title: f"Expression of {gene}"
  - Hover shows sample name and expression value
  
  **2b. `create_venn_diagram()`**:
  ```python
  def create_venn_diagram(
      gene_sets: Dict[str, set],  # comparison_name → set of DE gene names
      title: str = "DE Gene Overlap",
  ) -> go.Figure:
  ```
  - For 2 sets: classic 2-circle Venn using `go.Scatter` with circle shapes
  - For 3 sets: 3-circle Venn
  - For 4+ sets: UpSet-style horizontal bar chart (intersection sizes)
  - Show counts in each region
  - Use Plotly shapes for circles, annotations for counts
  
  **2c. `create_normalization_comparison_plot()`** (Medium priority bonus):
  ```python
  def create_normalization_comparison_plot(
      raw_counts: pd.DataFrame,  # genes × samples
      normalized_counts: pd.DataFrame,  # genes × samples
      sample_conditions: Dict[str, str],
  ) -> go.Figure:
  ```
  - Side-by-side box plots: raw log2(counts+1) vs normalized
  - Use `plotly.subplots.make_subplots(1, 2)` 
  - Color boxes by condition

  **Must NOT do**:
  - Don't modify existing functions
  - Don't import matplotlib
  - Don't change function signatures of existing functions

  **Recommended Agent Profile**:
  - **Category**: `ultrabrain`
  - **Skills**: []
    - Pure Python/Plotly work, no special skills needed

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 3, 4)
  - **Blocks**: Task 6
  - **Blocked By**: None

  **References**:
  - `visualizations.py:251-346` — `create_pca_plot()` — pattern for using `sample_conditions` dict and `px.scatter`
  - `visualizations.py:421-471` — `create_correlation_heatmap()` — pattern for sample grouping
  - `visualizations.py:623-650` — `create_deconvolution_plot()` — pattern for stacked bar with conditions
  - `visualizations.py:1-16` — imports section (add `from plotly.subplots import make_subplots` if needed)

  **Acceptance Criteria**:
  ```bash
  # Verify new functions exist
  python3 -c "from visualizations import create_gene_expression_plot, create_venn_diagram, create_normalization_comparison_plot; print('OK')"
  # Should print OK
  
  # Quick smoke test for gene expression plot
  python3 -c "
  import pandas as pd
  from visualizations import create_gene_expression_plot
  df = pd.DataFrame({'GENE1': [10, 12, 8, 15]}, index=['S1','S2','S3','S4'])
  conds = {'S1':'Ctrl','S2':'Ctrl','S3':'Treat','S4':'Treat'}
  fig = create_gene_expression_plot(df, 'GENE1', conds, 'box')
  assert fig is not None
  print('PASS')
  "
  
  # Quick smoke test for venn
  python3 -c "
  from visualizations import create_venn_diagram
  sets = {'A': {'g1','g2','g3'}, 'B': {'g2','g3','g4'}}
  fig = create_venn_diagram(sets)
  assert fig is not None
  print('PASS')
  "
  
  # No matplotlib
  grep -c "matplotlib" visualizations.py
  # Should be 0
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES
  - Message: `feat(viz): add gene expression plot, Venn diagram, and normalization comparison`
  - Files: `visualizations.py`

---

- [ ] 3. Session Manager Module

  **What to do**:
  - Create new file `session_manager.py`
  - Implement save/load of analysis session state as JSON:
  
  ```python
  import json
  import pandas as pd
  from typing import Dict, Any
  from datetime import datetime
  
  class SessionManager:
      """Save and load RNA-seq analysis sessions."""
      
      @staticmethod
      def save_session(session_state: dict) -> str:
          """Serialize session state to JSON string.
          Saves: settings, comparisons, metadata, DE results summary (not full DataFrames).
          """
          
      @staticmethod
      def load_session(json_str: str) -> dict:
          """Deserialize JSON back to session state dict."""
          
      @staticmethod
      def get_session_summary(session_state: dict) -> dict:
          """Get human-readable summary of session."""
  ```
  
  - **What to serialize**: comparisons list, sample_metadata dict, conditions, padj/lfc thresholds, active_comparison, analysis_complete flag, timestamp
  - **What NOT to serialize**: full DataFrames (too large), figures (not serializable), parsed_result (re-upload needed)
  - Include metadata: platform version, save timestamp, data shape info
  - Handle edge cases: empty state, partial state

  **Must NOT do**:
  - Don't serialize Plotly figures or large DataFrames
  - Don't use pickle (security risk)
  - Don't add new pip dependencies

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2, 4)
  - **Blocks**: Task 7
  - **Blocked By**: None

  **References**:
  - `rnaseq_analysis_platform.py:52-85` — Session state keys (this is the full list of what's in session_state)
  - `rnaseq_analysis_platform.py:918-931` — ExportData construction shows what data is available

  **Acceptance Criteria**:
  ```bash
  python3 -c "
  from session_manager import SessionManager
  state = {
      'comparisons': [('Treatment', 'Control')],
      'sample_metadata': {'S1': {'condition': 'Control'}},
      'conditions': ['Control', 'Treatment'],
      'analysis_complete': True,
  }
  json_str = SessionManager.save_session(state)
  loaded = SessionManager.load_session(json_str)
  assert loaded['comparisons'] == [['Treatment', 'Control']]  # tuples become lists in JSON
  assert loaded['analysis_complete'] == True
  print('PASS')
  "
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES
  - Message: `feat: add session save/load manager`
  - Files: `session_manager.py`

---

- [ ] 4. Inline Help Text Constants

  **What to do**:
  - Add a `HELP_TEXTS` dict constant in `rnaseq_analysis_platform.py` (after imports, near top)
  - Contains help strings for every major widget and tab:
  
  ```python
  HELP_TEXTS = {
      "padj_threshold": "Adjusted p-value cutoff. Genes with padj below this are considered significant. Default 0.05 controls false discovery rate at 5%.",
      "lfc_threshold": "Log2 fold change cutoff. Genes must exceed this magnitude to be called up/down regulated. 1.0 = 2-fold change.",
      "volcano": "Volcano plots show statistical significance (-log10 p-value) vs biological significance (fold change). Points in upper corners are most interesting.",
      "ma_plot": "MA plots show fold change vs average expression. Helps identify expression-dependent bias.",
      "pca": "Principal Component Analysis projects high-dimensional expression data into 2D. Samples that cluster together have similar expression profiles.",
      "heatmap": "Clustered heatmap of top DE genes. Z-score normalized per gene. Rows clustered by similarity.",
      "enrichment": "Pathway enrichment tests if DE genes are overrepresented in known biological pathways (GO, KEGG, Reactome).",
      "gsea": "Gene Set Enrichment Analysis uses ALL genes ranked by expression change, not just significant ones. More sensitive than ORA.",
      "deconvolution": "Estimates cell type proportions from bulk RNA-seq using reference signatures. NNLS method.",
      "gene_expression": "View normalized expression of individual genes across conditions. Useful for validating DE results.",
      "venn": "Overlap between DE gene sets from different comparisons. Shared genes respond to multiple conditions.",
      "upload": "Upload a gene count matrix (CSV/TSV/Excel). Rows = genes, columns = samples. Raw counts preferred for DE analysis.",
      "metadata": "Assign each sample to an experimental condition (e.g., Control, Treatment). Need ≥2 samples per condition.",
  }
  ```

  **Must NOT do**:
  - Don't add help text to backend modules
  - Don't modify widget keys or IDs

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2, 3)
  - **Blocks**: Tasks 5, 6, 7
  - **Blocked By**: None

  **References**:
  - `rnaseq_analysis_platform.py:185-186` — Sidebar sliders (add `help=HELP_TEXTS["padj_threshold"]`)
  - `rnaseq_analysis_platform.py:654-659` — Tab definitions (add about sections inside each tab)

  **Acceptance Criteria**:
  ```bash
  grep -c "HELP_TEXTS" rnaseq_analysis_platform.py
  # Should be >= 10 (dict definition + usages)
  
  grep -c 'help=' rnaseq_analysis_platform.py
  # Should be >= 2 (at minimum the sliders)
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES (group with Task 5 or 6)
  - Message: `feat(ui): add inline help tooltips and contextual guidance`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 5. Per-Plot Download Buttons

  **What to do**:
  - Create a helper function in `rnaseq_analysis_platform.py`:
  
  ```python
  def plot_with_download(fig: go.Figure, key: str, title: str = ""):
      """Display Plotly chart with PNG/SVG download buttons."""
      st.plotly_chart(fig, use_container_width=True, key=f"chart_{key}")
      col1, col2, _ = st.columns([1, 1, 6])
      with col1:
          try:
              png_bytes = fig.to_image(format="png", scale=2, width=1200, height=800)
              st.download_button("📷 PNG", png_bytes, f"{key}.png", "image/png", key=f"dl_png_{key}")
          except Exception:
              pass  # kaleido not available
      with col2:
          try:
              svg_bytes = fig.to_image(format="svg", width=1200, height=800)
              st.download_button("📐 SVG", svg_bytes, f"{key}.svg", "image/svg+xml", key=f"dl_svg_{key}")
          except Exception:
              pass
  ```
  
  - Replace ALL `st.plotly_chart(fig, use_container_width=True)` calls in Step 4 with `plot_with_download(fig, unique_key)`
  - Affected locations (~12 calls):
    - Line 697: volcano
    - Line 702: MA plot
    - Line 706: heatmap
    - Line 708: correlation
    - Line 712: PCA
    - Line 725: GO enrichment dotplot
    - Line 733: KEGG enrichment dotplot
    - Line 751: additional DB dotplots
    - Line 782: GSEA NES bar
    - Line 848-850: deconvolution
    - Line 889: outlier plot
    - Line 909: gene panels

  **Must NOT do**:
  - Don't add downloads to non-results charts
  - Don't break existing chart display
  - Don't remove `use_container_width=True`

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 6, 7)
  - **Blocks**: Task 10
  - **Blocked By**: Task 1 (CSS should be in place for button styling)

  **References**:
  - `rnaseq_analysis_platform.py:694-712` — Volcano/MA/Heatmap/PCA chart displays
  - `rnaseq_analysis_platform.py:721-752` — Enrichment chart displays
  - `rnaseq_analysis_platform.py:768-782` — GSEA chart display
  - `rnaseq_analysis_platform.py:846-850` — Deconvolution chart display
  - `export_engine.py` — Has `export_figure()` but we don't need it; Plotly's `to_image()` is simpler per-plot

  **Acceptance Criteria**:
  ```bash
  # Verify helper function exists
  grep -c "def plot_with_download" rnaseq_analysis_platform.py
  # Should be 1
  
  # Verify it's used (replaced st.plotly_chart calls)
  grep -c "plot_with_download" rnaseq_analysis_platform.py
  # Should be >= 10
  
  # Old direct plotly_chart calls in Step 4 should be minimal
  # (some may remain for non-result charts)
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES
  - Message: `feat(ui): add per-plot PNG/SVG download buttons`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 6. New Tabs: Gene Expression + Venn Diagram

  **What to do**:
  - Add 2 new tabs to the results tabs list (line 654):
    - "🧬 Gene Expression" — after Gene Panels
    - "🔀 Venn Overlap" — after Gene Expression
  - Update the `st.tabs()` call to include new tab variables
  
  **Gene Expression Tab**:
  ```python
  with tab_gene_expr:
      st.subheader("Gene Expression Viewer")
      with st.expander("ℹ️ About this analysis", expanded=False):
          st.markdown(HELP_TEXTS["gene_expression"])
      
      # Gene selector — searchable selectbox from DE results
      if res.results_df is not None:
          all_genes = sorted(res.results_df["gene"].dropna().unique().tolist())
          selected_gene = st.selectbox("Select Gene", all_genes, 
                                        help="Type to search")
          plot_type = st.radio("Plot Type", ["Box", "Violin"], horizontal=True)
          
          # Get normalized counts
          norm_counts = res.log_normalized_counts
          if norm_counts is not None and selected_gene in norm_counts.columns:
              sample_conds = {s: m["condition"] for s, m in st.session_state["sample_metadata"].items()}
              fig = create_gene_expression_plot(
                  norm_counts, selected_gene, sample_conds, plot_type.lower()
              )
              plot_with_download(fig, f"gene_expr_{selected_gene}")
          else:
              st.warning(f"Gene '{selected_gene}' not found in normalized expression data.")
  ```
  
  **Venn Tab**:
  ```python
  with tab_venn:
      st.subheader("DE Gene Overlap")
      with st.expander("ℹ️ About this analysis", expanded=False):
          st.markdown(HELP_TEXTS["venn"])
      
      if len(st.session_state["de_results"]) >= 2:
          # Build gene sets from DE results
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
              
              # Show shared genes table
              all_sets = list(gene_sets.values())
              shared = set.intersection(*all_sets)
              if shared:
                  st.markdown(f"**{len(shared)} genes shared across all comparisons:**")
                  st.dataframe(pd.DataFrame(sorted(shared), columns=["Gene"]))
          else:
              st.info("Need at least 2 comparisons with significant DE genes.")
      else:
          st.info("Run multiple comparisons to see overlap. Currently only 1 comparison.")
  ```
  
  - Also add import for new visualization functions at top of file

  **Must NOT do**:
  - Don't modify existing tab content
  - Don't change tab order of existing tabs
  - Don't modify DE results data structure

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 5, 7)
  - **Blocks**: Task 8, 10
  - **Blocked By**: Task 2 (needs new viz functions), Task 4 (needs HELP_TEXTS)

  **References**:
  - `rnaseq_analysis_platform.py:654-659` — Current `st.tabs()` call — append new tabs here
  - `rnaseq_analysis_platform.py:18-39` — Import block — add new function imports
  - `rnaseq_analysis_platform.py:559-564` — `norm_counts` and `sample_conds` pattern — reuse in gene expr tab
  - `de_analysis.py` — `DEResult` dataclass has `log_normalized_counts` attribute (samples × genes)
  - `visualizations.py` — new `create_gene_expression_plot()` and `create_venn_diagram()` from Task 2

  **Acceptance Criteria**:
  ```bash
  # Verify new tabs exist
  grep -c "Gene Expression" rnaseq_analysis_platform.py
  # Should be >= 1
  
  grep -c "Venn" rnaseq_analysis_platform.py
  # Should be >= 1
  
  # Verify imports
  grep "create_gene_expression_plot" rnaseq_analysis_platform.py
  # Should appear in imports and usage
  
  grep "create_venn_diagram" rnaseq_analysis_platform.py
  # Should appear in imports and usage
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES
  - Message: `feat: add Gene Expression viewer and Venn diagram tabs`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 7. Session Save/Load UI

  **What to do**:
  - Add session save/load to sidebar, below the Settings section:
  
  ```python
  st.markdown("---")
  st.subheader("💾 Session")
  
  # Save
  if st.session_state.get("analysis_complete"):
      from session_manager import SessionManager
      json_str = SessionManager.save_session(dict(st.session_state))
      st.download_button(
          "Save Session",
          json_str.encode("utf-8"),
          f"rnaseq_session_{datetime.now().strftime('%Y%m%d_%H%M')}.json",
          "application/json",
          key="save_session"
      )
  
  # Load
  session_file = st.file_uploader("Load Session", type=["json"], key="load_session")
  if session_file:
      from session_manager import SessionManager
      loaded = SessionManager.load_session(session_file.read().decode("utf-8"))
      for key, val in loaded.items():
          if key in st.session_state:
              st.session_state[key] = val
      st.success("Session loaded! Note: re-upload data file to restore full analysis.")
      st.rerun()
  ```

  **Must NOT do**:
  - Don't auto-restore figures or DataFrames (too large for JSON)
  - Don't use pickle

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Tasks 5, 6)
  - **Blocks**: Task 10
  - **Blocked By**: Task 3 (needs SessionManager)

  **References**:
  - `rnaseq_analysis_platform.py:167-186` — Sidebar section where save/load goes
  - `session_manager.py` — SessionManager class from Task 3
  - `rnaseq_analysis_platform.py:52-85` — All session state keys

  **Acceptance Criteria**:
  ```bash
  grep -c "session_manager" rnaseq_analysis_platform.py
  # Should be >= 1
  
  grep -c "Save Session" rnaseq_analysis_platform.py
  # Should be >= 1
  
  grep -c "Load Session" rnaseq_analysis_platform.py
  # Should be >= 1
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES (group with Task 3)
  - Message: `feat: add session save/load to sidebar`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 8. Caching with @st.cache_data

  **What to do**:
  - Add `@st.cache_data` decorators to expensive operations. Since the analysis pipeline runs inside a button callback (lines 429-625), we need to extract cacheable functions:
  
  **8a. Cache file parsing** — wrap parser call:
  ```python
  @st.cache_data
  def cached_parse(file_bytes: bytes, filename: str) -> ParseResult:
      parser = RNASeqParser()
      return parser.parse_bytes(file_bytes, filename)
  ```
  Note: Check if `RNASeqParser` has a `parse_bytes` or similar. If not, use `io.BytesIO`.
  
  **8b. Cache DE analysis**:
  ```python
  @st.cache_data
  def cached_de_analysis(_counts_df, _meta_df, comparisons, design_factor):
      engine = DEAnalysisEngine()
      return engine.run_all_comparisons(_counts_df, _meta_df, comparisons, design_factor)
  ```
  Note: Use `_` prefix for unhashable args that should be hashed by id.
  
  **8c. Cache PCA computation**:
  ```python
  @st.cache_data  
  def cached_pca(expression_values, n_components):
      pca = PCA(n_components=n_components)
      return pca.fit_transform(expression_values), pca.explained_variance_ratio_
  ```
  
  **8d. Cache enrichment** — each API call:
  ```python
  @st.cache_data(ttl=3600)  # 1hr cache for API calls
  def cached_enrichment(genes_tuple, db_name):
      pe = PathwayEnrichment()
      method = getattr(pe, f'get_{db_name}_enrichment')
      return method(list(genes_tuple))
  ```
  
  - Replace inline calls in the analysis pipeline (lines 429-625) with cached versions
  - Use `hash_funcs` or `_` prefix for DataFrames

  **Must NOT do**:
  - Don't cache visualization generation (cheap and depends on settings)
  - Don't cache session state
  - Don't modify backend module files

  **Recommended Agent Profile**:
  - **Category**: `ultrabrain`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Task 9)
  - **Blocks**: Task 10
  - **Blocked By**: Task 6 (main app should be stable before adding caching)

  **References**:
  - `rnaseq_analysis_platform.py:429-625` — Full analysis pipeline block (wrap functions here)
  - `rnaseq_analysis_platform.py:200-230` — File upload section (cache parsing here)
  - `rnaseq_parser.py` — `RNASeqParser.parse()` method signature (check if accepts bytes/BytesIO)
  - `de_analysis.py` — `DEAnalysisEngine.run_all_comparisons()` signature
  - Streamlit docs: `@st.cache_data` supports `hash_funcs`, `ttl`, `_` prefix for unhashable params

  **Acceptance Criteria**:
  ```bash
  grep -c "st.cache_data" rnaseq_analysis_platform.py
  # Should be >= 3
  
  grep -c "def cached_" rnaseq_analysis_platform.py
  # Should be >= 2
  
  # Ensure no circular imports or runtime errors
  python3 -c "import rnaseq_analysis_platform" 2>&1 | head -5
  # Should not show import errors (may show Streamlit context warning, that's OK)
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES
  - Message: `perf: add st.cache_data to parsing, DE analysis, and enrichment`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 9. Styled Step Indicator (Medium Priority Bonus)

  **What to do**:
  - Replace the plain markdown step indicator (lines 172-180) with HTML-styled version:
  ```python
  # Replace plain markdown with styled HTML
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
  ```
  - CSS classes for `.step-indicator`, `.step.completed`, `.step.active`, `.step.pending` should already exist from Task 1

  **Must NOT do**:
  - Don't break step navigation logic
  - Don't change step numbering

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 3 (with Task 8)
  - **Blocks**: Task 10
  - **Blocked By**: Task 1 (needs CSS classes)

  **References**:
  - `rnaseq_analysis_platform.py:172-180` — Current step indicator code to replace
  - Task 1 CSS — should define `.step-indicator` classes

  **Acceptance Criteria**:
  ```bash
  grep -c "step-indicator" rnaseq_analysis_platform.py
  # Should be >= 1
  
  pytest tests/ -x -q 2>&1 | tail -5
  ```

  **Commit**: YES (group with Task 1)
  - Message: `feat(ui): styled step indicator with progress states`
  - Files: `rnaseq_analysis_platform.py`

---

- [ ] 10. Final Verification

  **What to do**:
  - Run full test suite: `pytest tests/ -v`
  - Verify 147+ tests pass (same baseline)
  - Verify no matplotlib imports: `grep -r "matplotlib" *.py` (should be 0 in main app and visualizations)
  - Verify file count: only `session_manager.py` is new, rest are edits
  - Quick manual check: `streamlit run rnaseq_analysis_platform.py` → load demo → verify CSS, tabs, downloads

  **Must NOT do**:
  - Don't modify test files
  - Don't skip failing tests

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`playwright`]
    - `playwright`: For browser verification of the running Streamlit app

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (final, sequential)
  - **Blocks**: None (final task)
  - **Blocked By**: Tasks 5, 6, 7, 8, 9

  **References**:
  - `tests/` directory — all test files
  - `rnaseq_analysis_platform.py` — final state of main app

  **Acceptance Criteria**:
  ```bash
  # Full test suite
  pytest tests/ -v 2>&1 | tail -20
  # Should show 147+ passed, 3 or fewer failures (pre-existing PyDESeq2 issues)
  
  # No matplotlib in app or viz
  grep -c "matplotlib" rnaseq_analysis_platform.py visualizations.py
  # Both should be 0
  
  # New file exists
  test -f session_manager.py && echo "EXISTS" || echo "MISSING"
  
  # CSS is injected
  grep -c "<style>" rnaseq_analysis_platform.py
  # >= 1
  
  # Cache is active
  grep -c "cache_data" rnaseq_analysis_platform.py
  # >= 3
  
  # Download buttons exist
  grep -c "plot_with_download" rnaseq_analysis_platform.py
  # >= 10
  
  # New tabs exist
  grep -c "Gene Expression" rnaseq_analysis_platform.py
  # >= 1
  ```

  **Commit**: NO (verification only)

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `feat(ui): add professional CSS theme with genomics color palette` | rnaseq_analysis_platform.py | pytest tests/ -x |
| 2 | `feat(viz): add gene expression plot, Venn diagram, and normalization comparison` | visualizations.py | pytest tests/ -x |
| 3+7 | `feat: add session save/load` | session_manager.py, rnaseq_analysis_platform.py | pytest tests/ -x |
| 4+5 | `feat(ui): add help tooltips and per-plot download buttons` | rnaseq_analysis_platform.py | pytest tests/ -x |
| 6 | `feat: add Gene Expression viewer and Venn diagram tabs` | rnaseq_analysis_platform.py | pytest tests/ -x |
| 8 | `perf: add st.cache_data caching` | rnaseq_analysis_platform.py | pytest tests/ -x |
| 9 | `feat(ui): styled step indicator` | rnaseq_analysis_platform.py | pytest tests/ -x |

---

## Success Criteria

### Verification Commands
```bash
pytest tests/ -v                    # 147+ pass
grep -c "matplotlib" rnaseq_analysis_platform.py visualizations.py  # 0, 0
grep -c "cache_data" rnaseq_analysis_platform.py    # >= 3
grep -c "plot_with_download" rnaseq_analysis_platform.py  # >= 10
grep -c "<style>" rnaseq_analysis_platform.py       # >= 1
grep -c "HELP_TEXTS" rnaseq_analysis_platform.py    # >= 10
python3 -c "from session_manager import SessionManager; print('OK')"  # OK
python3 -c "from visualizations import create_gene_expression_plot, create_venn_diagram; print('OK')"  # OK
```

### Final Checklist
- [ ] Professional CSS theme applied
- [ ] Every chart has PNG/SVG download
- [ ] Gene expression viewer works
- [ ] Venn diagram shows for multi-comparison
- [ ] Session save/load functional
- [ ] Caching active on expensive ops
- [ ] Help tooltips on key widgets
- [ ] Step indicator styled
- [ ] 147+ tests pass
- [ ] Zero matplotlib in UI code
