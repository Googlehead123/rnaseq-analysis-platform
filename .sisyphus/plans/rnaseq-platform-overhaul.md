# RNA-seq Analysis Platform Overhaul

## TL;DR

> **Quick Summary**: Fix critical KeyError bug blocking analysis, overhaul UI from sidebar-centric to wizard-style main-area workflow inspired by iDEP/DEBrowser/Shiny-Seq, and enhance scientific visualization/QC patterns.
> 
> **Deliverables**:
> - Bug fix for KeyError: 'gene' in DE analysis pipeline
> - Wizard-style UI with step progression in main area
> - Spreadsheet-style metadata editor (st.data_editor)
> - Pre-analysis QC visualizations
> - Demo dataset, gene search, inline help
> - Enhanced scientific visualizations and summary dashboard
> - Updated test suite
> 
> **Estimated Effort**: Large
> **Parallel Execution**: YES - 3 waves
> **Critical Path**: Task 1 (bug fix) → Task 3 (UI overhaul) → Task 6 (integration testing)

---

## Context

### Original Request
User wants 4 things: (1) apply scientific best practices for analysis/viz, (2) UI overhaul making condition selection prominent, (3) fix KeyError: 'gene' bug, (4) use iDEP/DEBrowser/Shiny-Seq/GENAVi/START as reference.

### Research Findings
- Bug is in `de_analysis.py:128` — `reset_index()` preserves named index, breaking downstream column lookups for "gene"
- UI is 100% sidebar-centric; main area is blank until analysis completes
- Reference platforms all use wizard/stepper patterns with QC before analysis
- Streamlit has `st.data_editor`, `st.tabs`, `st.status` that map perfectly to needed patterns

---

## Work Objectives

### Core Objective
Transform the platform from a sidebar-heavy, bug-ridden prototype into a polished wizard-style RNA-seq analysis tool.

### Definition of Done
- [ ] `pytest tests/ -v` — all tests pass
- [ ] Upload test data → assign metadata → select comparison → run analysis → view results: no errors
- [ ] KeyError: 'gene' never occurs regardless of index naming
- [ ] Condition selection is visible in main area without scrolling

### Must NOT Have (Guardrails)
- No new analysis engines or DE methods
- No new frameworks (stay Streamlit/Plotly/PyDESeq2/GSEApy)
- No database or auth features
- No unnecessary abstractions — keep single-app simplicity
- No breaking changes to module APIs (de_analysis, pathway_enrichment, visualizations, gene_panels)

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest)
- **User wants tests**: YES (tests-after)
- **Framework**: pytest

### Automated Verification

All acceptance criteria use automated commands:
```bash
# Regression suite
pytest tests/ -v

# Smoke test — app launches
timeout 10 streamlit run rnaseq_analysis_platform.py --server.headless true 2>&1 | grep "You can now view"
```

For UI verification: Playwright browser automation via playwright skill.

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately):
├── Task 1: Fix KeyError bug (de_analysis.py)
├── Task 2: Add demo dataset + gene search utility
└── Task 4: Enhanced visualizations module

Wave 2 (After Wave 1):
├── Task 3: UI overhaul (depends: 1, 2)
└── Task 5: Scientific QC module (depends: 4)

Wave 3 (After Wave 2):
└── Task 6: Test updates + integration verification (depends: all)
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 3, 6 | 2, 4 |
| 2 | None | 3 | 1, 4 |
| 3 | 1, 2 | 6 | 5 |
| 4 | None | 5, 6 | 1, 2 |
| 5 | 4 | 6 | 3 |
| 6 | All | None | None |

---

## TODOs

- [ ] 1. Fix KeyError: 'gene' bug in DE analysis pipeline

  **What to do**:
  - In `de_analysis.py` line ~128, replace the fragile rename:
    ```python
    # BEFORE (broken):
    results_df = stat_res.results_df.reset_index().rename(columns={"index": "gene"})
    
    # AFTER (robust):
    results_df = stat_res.results_df.copy()
    results_df.index.name = None
    results_df = results_df.reset_index()
    results_df.columns = ["gene"] + list(results_df.columns[1:])
    ```
  - Audit ALL downstream consumers for hardcoded "gene" column access and add defensive lowercasing:
    - `pathway_enrichment.py:73` — `sig["gene"]`
    - `visualizations.py:39` — volcano plot gene column
    - `gene_panels.py` — any gene column references
  - Add a helper function `normalize_gene_column(df)` that ensures the gene column is always named "gene" regardless of input

  **Must NOT do**:
  - Do not change the DEResult dataclass structure
  - Do not change function signatures

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]
    - `git-master`: Atomic commit for critical bug fix

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2, 4)
  - **Blocks**: Tasks 3, 6
  - **Blocked By**: None

  **References**:
  - `de_analysis.py:128` — Bug location: `reset_index().rename(columns={"index": "gene"})`
  - `pathway_enrichment.py:73` — Crash point: `sig["gene"].tolist()`
  - `visualizations.py:39` — Volcano plot gene column usage
  - `gene_panels.py` — Gene column references in panel analysis
  - `tests/test_integration.py` — Existing integration tests to verify fix

  **Acceptance Criteria**:
  ```bash
  # Existing tests still pass
  pytest tests/ -v
  
  # Specific bug verification — create a test that uses named index
  pytest tests/ -v -k "gene" 
  ```

  **Commit**: YES
  - Message: `fix(de-analysis): handle named DataFrame index in DE results to prevent KeyError`
  - Files: `de_analysis.py`, `pathway_enrichment.py`, `visualizations.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 2. Add demo dataset and gene search utility

  **What to do**:
  - Create `demo_data.py` module:
    - Function `load_demo_dataset()` → returns a count matrix DataFrame + sample metadata DataFrame
    - Use `generate_test_data.py` as reference for realistic RNA-seq data generation
    - Include 6 samples (3 treatment, 3 control), ~500 genes with realistic distributions
    - Include known DE genes so demo results are meaningful
  - Create `gene_search.py` utility:
    - Function `search_genes(results_df, query)` → filtered DataFrame matching gene name pattern
    - Support partial matching, case-insensitive
    - Return highlighted matches for display

  **Must NOT do**:
  - Do not embed large data files — generate programmatically
  - Do not add external data dependencies

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 4)
  - **Blocks**: Task 3
  - **Blocked By**: None

  **References**:
  - `generate_test_data.py` — Existing test data generation patterns
  - `de_analysis.py:DEResult` — Data contract for downstream compatibility
  - `rnaseq_parser.py` — Parser output format that demo data must match

  **Acceptance Criteria**:
  ```bash
  # Demo data loads without error
  python -c "from demo_data import load_demo_dataset; df, meta = load_demo_dataset(); print(df.shape, meta.shape)"
  # Assert: prints shapes like (500, 6) (6, 2)
  
  # Gene search works
  python -c "from gene_search import search_genes; import pandas as pd; df = pd.DataFrame({'gene':['TP53','BRCA1','MYC']}); print(search_genes(df, 'br'))"
  # Assert: returns DataFrame with BRCA1
  ```

  **Commit**: YES
  - Message: `feat: add demo dataset and gene search utility`
  - Files: `demo_data.py`, `gene_search.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 3. UI overhaul — wizard-style main area workflow

  **What to do**:
  
  Restructure `rnaseq_analysis_platform.py` from sidebar-centric to wizard-style main area:

  **Step 1 — Upload & Preview** (main area):
  - File uploader in main area (prominent, centered)
  - "Load Demo Data" button alongside uploader
  - On upload: show data preview table (first 10 rows, basic stats)
  - Show sample count, gene count, data format detected

  **Step 2 — Sample Metadata** (main area):
  - Replace N selectboxes with `st.data_editor` showing editable table:
    - Columns: Sample Name | Condition (dropdown)
    - Bulk tools: "Set all selected to..." dropdown above table
  - Condition summary: show count per condition
  - Validation: highlight if <2 samples per condition

  **Step 3 — Comparison & Parameters** (main area, PROMINENT):
  - Large, clear comparison selector: Test Condition vs Reference Condition
  - Visual indication of which samples are in each group
  - Analysis parameters (FC threshold, p-value) with inline help tooltips
  - "Run Analysis" button — large, primary-styled

  **Step 4 — Results** (main area, tabbed):
  - Use `st.tabs`: Summary | DE Results | Volcano | Heatmap | PCA | Enrichment | Gene Panels | Export
  - Summary tab: key metrics dashboard (total DE genes, up/down counts, top genes)
  - Gene search bar at top of DE Results tab
  - All existing result displays reorganized into tabs

  **Sidebar** (minimal):
  - App title/logo
  - Step indicator (Steps 1-4 with current step highlighted)
  - Global settings only (theme, export preferences)

  **Step navigation**:
  - Use `st.session_state` to track current step
  - "Next" / "Back" buttons between steps
  - Step indicator showing progress (custom CSS or emoji-based: ✅ Step 1 → 🔵 Step 2 → ⚪ Step 3 → ⚪ Step 4)
  - Validate before allowing next step

  **Additional UI improvements**:
  - `st.status` container for analysis progress feedback
  - `st.cache_data` on DE analysis and enrichment calls
  - Error messages in main area (not buried in sidebar)
  - Inline help with `st.popover` or tooltips for technical parameters

  **Must NOT do**:
  - Do not change module APIs (de_analysis, pathway_enrichment, etc.)
  - Do not remove any existing functionality
  - Do not change the export engine interface
  - Do not add multi-page Streamlit — keep single page with wizard steps

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: [`frontend-ui-ux`, `playwright`]
    - `frontend-ui-ux`: UI/UX design decisions, layout, visual hierarchy
    - `playwright`: Browser-based verification of the wizard flow

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 2
  - **Blocks**: Task 6
  - **Blocked By**: Tasks 1, 2

  **References**:
  - `rnaseq_analysis_platform.py:146-263` — Current upload logic to preserve
  - `rnaseq_analysis_platform.py:316-389` — Current metadata assignment to replace with st.data_editor
  - `rnaseq_analysis_platform.py:391-424` — Current comparison selection to move to main area
  - `rnaseq_analysis_platform.py:426-642` — Current analysis runner to integrate into step 3
  - `rnaseq_analysis_platform.py:644-843` — Current results display to reorganize into tabs
  - `export_engine.py` — Export interface (do not change, just wire into Export tab)
  - `demo_data.py` (Task 2) — Demo dataset loading for "Load Demo Data" button
  - `gene_search.py` (Task 2) — Gene search for DE Results tab

  **Acceptance Criteria**:
  ```
  # Agent executes via playwright browser automation:
  1. Run: streamlit run rnaseq_analysis_platform.py --server.headless true
  2. Navigate to: http://localhost:8501
  3. Assert: Main area shows upload step (not blank welcome page)
  4. Assert: "Load Demo Data" button visible
  5. Click: "Load Demo Data"
  6. Assert: Data preview table appears with gene count and sample count
  7. Click: "Next" button
  8. Assert: Metadata editor table visible (st.data_editor)
  9. Assert: Condition dropdowns work in table
  10. Click: "Next" button  
  11. Assert: Comparison selector visible in MAIN AREA (not sidebar)
  12. Assert: Test vs Reference dropdowns visible and prominent
  13. Click: "Run Analysis"
  14. Assert: Results tabs appear (Summary, DE Results, Volcano, etc.)
  15. Assert: No KeyError or any error during entire flow
  16. Screenshot: .sisyphus/evidence/task-3-wizard-flow.png
  ```

  **Commit**: YES
  - Message: `feat(ui): overhaul to wizard-style main area workflow with step navigation`
  - Files: `rnaseq_analysis_platform.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 4. Enhanced visualizations and scientific best practices

  **What to do**:
  - Enhance `visualizations.py`:
    - Volcano plot: add gene labels for top N significant genes, color by significance category
    - Heatmap: add dendrogram clustering, color scale legend, sample grouping bars
    - PCA: add variance explained labels on axes, confidence ellipses per group, loadings option
    - NEW: MA plot (log2FC vs mean expression) — standard RNA-seq QC plot
    - NEW: Sample correlation heatmap
  - Add summary statistics functions:
    - `compute_de_summary(results_df)` → dict with total DE, up/down counts, top genes by FC and significance
    - `format_summary_dashboard(summary)` → Streamlit-renderable metrics

  **Must NOT do**:
  - Do not change existing function signatures (add new optional parameters only)
  - Do not remove existing plot functions

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2)
  - **Blocks**: Tasks 5, 6
  - **Blocked By**: None

  **References**:
  - `visualizations.py:1-285` — All current visualization functions
  - `de_analysis.py:DEResult` — Data contract (results_df columns, normalized_counts format)
  - Plotly docs for dendrogram, annotations, confidence ellipses

  **Acceptance Criteria**:
  ```bash
  # Import and call enhanced functions
  python -c "
  from visualizations import create_volcano_plot, create_heatmap, create_pca_plot
  import pandas as pd, numpy as np
  # Create test data
  df = pd.DataFrame({'gene': [f'G{i}' for i in range(100)], 'log2FoldChange': np.random.randn(100), 'padj': np.random.uniform(0,1,100)})
  fig = create_volcano_plot(df)
  print('Volcano OK:', fig is not None)
  "
  
  # New functions exist
  python -c "from visualizations import create_ma_plot, create_correlation_heatmap, compute_de_summary; print('All imports OK')"
  ```

  **Commit**: YES
  - Message: `feat(viz): enhance plots with labels, MA plot, correlation heatmap, summary stats`
  - Files: `visualizations.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 5. Pre-analysis QC visualizations

  **What to do**:
  - Add QC functions to `visualizations.py` (or new `qc_plots.py` if cleaner):
    - `create_library_size_barplot(count_df)` — bar plot of total counts per sample
    - `create_count_distribution_boxplot(count_df)` — box plots of log counts per sample
    - `create_gene_detection_plot(count_df)` — genes detected (count > 0) per sample
    - `create_sample_similarity_heatmap(count_df)` — Spearman correlation between samples
  - These are shown in Step 1 of the wizard (after upload, before metadata)

  **Must NOT do**:
  - Do not duplicate logic already in visualizations.py

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Task 3)
  - **Blocks**: Task 6
  - **Blocked By**: Task 4

  **References**:
  - `visualizations.py` — Existing Plotly patterns to follow
  - `de_analysis.py` — Normalized counts format for QC input
  - iDEP/DEBrowser QC patterns: library size, count distribution, gene detection are standard

  **Acceptance Criteria**:
  ```bash
  python -c "
  from qc_plots import create_library_size_barplot, create_count_distribution_boxplot, create_gene_detection_plot, create_sample_similarity_heatmap
  import pandas as pd, numpy as np
  df = pd.DataFrame(np.random.poisson(100, (500,6)), columns=[f'S{i}' for i in range(6)], index=[f'G{i}' for i in range(500)])
  for fn in [create_library_size_barplot, create_count_distribution_boxplot, create_gene_detection_plot, create_sample_similarity_heatmap]:
      fig = fn(df)
      print(f'{fn.__name__} OK:', fig is not None)
  "
  ```

  **Commit**: YES
  - Message: `feat: add pre-analysis QC visualizations (library size, distributions, detection, similarity)`
  - Files: `qc_plots.py` or `visualizations.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 6. Update tests and final integration verification

  **What to do**:
  - Add test for the KeyError bug fix:
    - Test DE analysis with named index DataFrame (e.g., `df.index.name = "Gene"`)
    - Assert results_df has "gene" column regardless of input index name
  - Add tests for new modules:
    - `test_demo_data.py` — demo dataset loads, has correct shape, correct columns
    - `test_gene_search.py` — search returns correct matches, case-insensitive
    - `test_qc_plots.py` — QC plot functions return Plotly figures
  - Add test for enhanced visualizations:
    - MA plot, correlation heatmap, summary stats produce valid output
  - Run full integration test:
    - Load demo data → assign metadata → run DE → check all result types present
  - Update `test_comparison_ui.py` if UI flow changed

  **Must NOT do**:
  - Do not modify test infrastructure (conftest, pytest config)
  - Do not remove existing tests

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]
    - `git-master`: Final commit with all tests passing

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (sequential, after all other tasks)
  - **Blocks**: None
  - **Blocked By**: Tasks 1, 2, 3, 4, 5

  **References**:
  - `tests/test_integration.py` — Existing integration test patterns
  - `tests/test_comparison_ui.py` — UI test patterns
  - `tests/conftest.py` — Shared fixtures
  - `tests/data/` — Test data files

  **Acceptance Criteria**:
  ```bash
  # All tests pass including new ones
  pytest tests/ -v
  # Assert: 0 failures, 0 errors
  
  # Specifically the bug fix test
  pytest tests/ -v -k "named_index or gene_column"
  # Assert: PASSED
  ```

  **Commit**: YES
  - Message: `test: add tests for bug fix, demo data, gene search, QC plots, and enhanced viz`
  - Files: `tests/test_de_bugfix.py`, `tests/test_demo_data.py`, `tests/test_gene_search.py`, `tests/test_qc_plots.py`
  - Pre-commit: `pytest tests/ -v`

---

## Commit Strategy

| After Task | Message | Key Files | Verification |
|------------|---------|-----------|--------------|
| 1 | `fix(de-analysis): handle named DataFrame index` | de_analysis.py, pathway_enrichment.py, visualizations.py | pytest tests/ -v |
| 2 | `feat: add demo dataset and gene search utility` | demo_data.py, gene_search.py | pytest tests/ -v |
| 4 | `feat(viz): enhance plots, add MA plot, summary stats` | visualizations.py | pytest tests/ -v |
| 3 | `feat(ui): wizard-style main area workflow` | rnaseq_analysis_platform.py | pytest tests/ -v + playwright |
| 5 | `feat: pre-analysis QC visualizations` | qc_plots.py | pytest tests/ -v |
| 6 | `test: comprehensive test updates` | tests/*.py | pytest tests/ -v |

---

## Success Criteria

### Verification Commands
```bash
# All tests pass
pytest tests/ -v

# App launches
timeout 10 streamlit run rnaseq_analysis_platform.py --server.headless true 2>&1 | grep "You can now view"

# No import errors for new modules
python -c "import demo_data, gene_search, qc_plots; print('All modules OK')"
```

### Final Checklist
- [ ] KeyError: 'gene' fixed — works with any index naming
- [ ] Wizard UI with 4 steps in main area
- [ ] Condition selection prominent in Step 3 (not buried in sidebar)
- [ ] Demo data loadable with one click
- [ ] QC plots shown before analysis
- [ ] Results in tabbed layout
- [ ] Gene search functional in results
- [ ] All existing functionality preserved
- [ ] All tests pass
