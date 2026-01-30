# RNA-seq Platform Enhancements: Multi-file Upload, Tab Explanations, Enhanced Interpretations

## TL;DR

> **Quick Summary**: Add multi-file upload with merge, educational explanations to all 15 tabs, and deeply enhanced InterpretationEngine with context-aware biological insights.
> 
> **Deliverables**:
> - Multi-file upload UI + merge utility (`multi_file_merger.py`)
> - Educational HELP_TEXTS wired to all 15 tabs
> - Enhanced `interpretation_engine.py` with 6+ new interpretation methods
> 
> **Estimated Effort**: Medium
> **Parallel Execution**: YES - 3 waves
> **Critical Path**: Task 1 (merger utility) → Task 2 (upload UI wiring) → Task 6 (integration test)

---

## Context

### Original Request
Three improvements to the RNA-seq Streamlit platform:
1. Multi-file upload with merge into single count matrix
2. Rich educational tab explanations for first-time RNA-seq users
3. Enhanced InterpretationEngine with deeper, data-specific insights

### Research Findings
- `st.file_uploader(accept_multiple_files=True)` returns `list[UploadedFile]`, each BytesIO-like
- `RNASeqParser.parse()` requires a file path (tempfile needed per file)
- `ParseResult.expression_df` is samples × genes (canonical shape)
- HELP_TEXTS has 13 keys; only 2 wired to tabs (gene_expression, venn via `st.caption`); 7 keys exist unused; 6 tabs have no keys at all
- InterpretationEngine has only `interpret_de_results()` (4 insights) and `interpret_enrichment()` (4 keyword categories)
- Available session_state for richer interpretations: `de_results`, `pca_coordinates`, `outlier_info`, `deconvolution_results`, `gsea_results`

### Gap Analysis (Metis-equivalent review)
**Addressed in plan:**
- Gene ID mismatch across files → validation + user warning
- Duplicate sample names → suffix with filename
- Mixed data types (Raw vs Normalized) → enforce same DataType
- Orientation conflicts → parse each file individually first (parser handles this)
- Session state clearing on re-upload → explicit reset of downstream state
- `conftest.py` mocks single UploadedFile → new tests needed

---

## Work Objectives

### Core Objective
Enhance the RNA-seq platform to support multi-file input workflows, provide educational context for every analysis tab, and generate deeper automated biological insights from actual data.

### Concrete Deliverables
- `multi_file_merger.py`: New utility for merging multiple ParseResults
- Modified `rnaseq_analysis_platform.py`: Multi-file upload UI, expanded HELP_TEXTS, tab explanation blocks, wiring enhanced interpretations
- Modified `interpretation_engine.py`: 6+ new interpretation methods (top genes, QC, PCA, deconvolution, GSEA, misinterpretation warnings)

### Definition of Done
- [ ] Multiple CSV/Excel files can be uploaded and merged into single expression matrix
- [ ] All 15 tabs display educational explanation blocks
- [ ] InterpretationEngine generates 10+ insight types per comparison
- [ ] `pytest tests/ -v` passes 142+ tests (no regressions)

### Must Have
- Inner join merge by gene ID with overlap report
- Duplicate sample name handling (suffix with filename)
- DataType consistency validation across files
- Educational content: purpose, "how to read", "what to look for", pitfalls per tab
- Top DE genes by name in interpretation text
- Misinterpretation warnings (9 common pitfalls)

### Must NOT Have (Guardrails)
- Do NOT modify `de_analysis.py`, `pathway_enrichment.py`, `visualizations.py`, `gene_panels.py`, `export_engine.py`, `rnaseq_parser.py`
- Do NOT use matplotlib — all Plotly
- Do NOT add new Python package dependencies
- Do NOT change the existing 15-tab structure or tab order
- Do NOT break the single-file upload path (must still work)
- Do NOT over-abstract — keep merge utility simple and focused

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest, 12 test files)
- **User wants tests**: Tests-after (must not break 142 existing)
- **Framework**: pytest

### Verification Approach
Each task includes:
1. `pytest tests/ -v` — regression check (142+ pass)
2. Manual Streamlit verification via `streamlit run rnaseq_analysis_platform.py`
3. Specific functional checks described per task

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately — all independent):
├── Task 1: multi_file_merger.py (new file, no deps)
├── Task 3: Expand HELP_TEXTS dict + tab explanation content
└── Task 4: Enhanced interpretation_engine.py methods

Wave 2 (After Wave 1):
├── Task 2: Wire multi-file upload UI (depends: Task 1)
└── Task 5: Wire tab explanations + interpretations into UI (depends: Tasks 3, 4)

Wave 3 (After Wave 2):
└── Task 6: Integration testing + regression check (depends: all)
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 2 | 3, 4 |
| 3 | None | 5 | 1, 4 |
| 4 | None | 5 | 1, 3 |
| 2 | 1 | 6 | 5 |
| 5 | 3, 4 | 6 | 2 |
| 6 | 2, 5 | None | None (final) |

---

## TODOs

- [ ] 1. Create multi-file merge utility (`multi_file_merger.py`)

  **What to do**:
  - Create new file `multi_file_merger.py` in project root
  - Implement `merge_parse_results(results: List[ParseResult], filenames: List[str]) -> ParseResult`
  - Logic:
    1. Validate all ParseResults have same `data_type` — raise ValueError if mixed
    2. Filter out results with `expression_df is None`
    3. For each result, get `expression_df` (samples × genes)
    4. Handle duplicate sample names: if sample "X" appears in multiple files, rename to "X_filename1", "X_filename2"
    5. Concatenate all expression_dfs along axis=0 (stack samples)
    6. For gene overlap: inner join on columns (genes). Warn if overlap < 80%
    7. Return merged ParseResult with combined expression_df, merged warnings, data_type from first file
  - Implement `validate_merge_compatibility(results: List[ParseResult]) -> Dict` returning validation report
  - Handle edge cases: empty files (skip with warning), single file (return as-is), zero gene overlap (raise ValueError)

  **Must NOT do**:
  - Do NOT modify `rnaseq_parser.py`
  - Do NOT add external dependencies
  - Do NOT handle gene ID mapping/conversion (out of scope)

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: Pure Python utility, ~100 lines, no UI or complex tooling needed

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 3, 4)
  - **Blocks**: Task 2
  - **Blocked By**: None

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py` (lines 12-30) — `ParseResult` dataclass definition with all fields needed for return value
  - `rnaseq_parser.py` (lines 35-50) — `DataType` enum definition (RAW_COUNTS, NORMALIZED, PRE_ANALYZED)

  **API/Type References**:
  - `ParseResult.expression_df` — DataFrame with samples as rows (index), genes as columns
  - `ParseResult.data_type` — `DataType` enum value
  - `ParseResult.warnings` — `List[str]`
  - `ParseResult.can_run_de` — bool (True only for RAW_COUNTS)

  **WHY Each Reference Matters**:
  - `ParseResult` dataclass: You must construct a valid return object matching all expected fields
  - `DataType` enum: Needed for consistency validation — all files must share same type
  - `expression_df` shape: Must understand samples × genes orientation to concatenate correctly (axis=0 = stack samples)

  **Acceptance Criteria**:

  ```bash
  # Verify file exists and imports
  python3 -c "from multi_file_merger import merge_parse_results, validate_merge_compatibility; print('Import OK')"
  # Assert: prints "Import OK"

  # Verify merge logic with inline test
  python3 -c "
  import pandas as pd
  from rnaseq_parser import ParseResult, DataType
  from multi_file_merger import merge_parse_results
  
  df1 = pd.DataFrame({'GeneA': [100, 200], 'GeneB': [300, 400]}, index=['S1', 'S2'])
  df2 = pd.DataFrame({'GeneA': [500, 600], 'GeneC': [700, 800]}, index=['S3', 'S4'])
  r1 = ParseResult(expression_df=df1, normalized_df=None, de_results_df=None, data_type=DataType.RAW_COUNTS, can_run_de=True, warnings=[], dropped_columns=[], gene_column_source='test', needs_user_input=False, gene_column_candidates=[], data_types_detected=[DataType.RAW_COUNTS])
  r2 = ParseResult(expression_df=df2, normalized_df=None, de_results_df=None, data_type=DataType.RAW_COUNTS, can_run_de=True, warnings=[], dropped_columns=[], gene_column_source='test', needs_user_input=False, gene_column_candidates=[], data_types_detected=[DataType.RAW_COUNTS])
  merged = merge_parse_results([r1, r2], ['file1.csv', 'file2.csv'])
  assert merged.expression_df.shape == (4, 1), f'Expected (4,1) got {merged.expression_df.shape}'  # inner join: only GeneA
  assert list(merged.expression_df.index) == ['S1', 'S2', 'S3', 'S4']
  print('Merge test PASSED')
  "
  # Assert: prints "Merge test PASSED"
  ```

  **Commit**: YES
  - Message: `feat(merger): add multi-file merge utility for combining count matrices`
  - Files: `multi_file_merger.py`
  - Pre-commit: `python3 -c "from multi_file_merger import merge_parse_results; print('OK')"`

---

- [ ] 2. Wire multi-file upload UI in main app

  **What to do**:
  - In `rnaseq_analysis_platform.py`, modify Step 1 upload section (lines 441-520):
    1. Change `st.file_uploader` to `accept_multiple_files=True`
    2. Add logic to handle `list[UploadedFile]`:
       - If 1 file: existing single-file path (unchanged behavior)
       - If 2+ files: parse each via tempfile → `RNASeqParser().parse()`, then call `merge_parse_results()`
    3. Show merge report: number of files, samples per file, gene overlap percentage, any warnings
    4. Display merged preview with `st.dataframe()`
    5. Clear downstream state on re-upload: `de_results`, `enrichment_results`, `gsea_results`, `interpretations`, `deconvolution_results`, `sample_metadata`, `comparisons`
  - Add import for `merge_parse_results` from `multi_file_merger`
  - Add progress indicator during multi-file parsing (`st.progress()`)

  **Must NOT do**:
  - Do NOT change the demo data loading path
  - Do NOT modify RNASeqParser itself
  - Do NOT break single-file upload (if user uploads 1 file, behavior must be identical to current)

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: []
  - **Reason**: Streamlit UI modifications with careful state management

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Task 5)
  - **Parallel Group**: Wave 2
  - **Blocks**: Task 6
  - **Blocked By**: Task 1

  **References**:

  **Pattern References**:
  - `rnaseq_analysis_platform.py:441-520` — Current single-file upload flow (the code to modify)
  - `rnaseq_analysis_platform.py:449-475` — Demo data loading pattern (DO NOT TOUCH, but follow similar session_state reset pattern)

  **API/Type References**:
  - `multi_file_merger.merge_parse_results(results, filenames)` — The new utility from Task 1
  - `RNASeqParser().parse(file_path)` — Returns `ParseResult`
  - `st.file_uploader(label, type=[], accept_multiple_files=True)` — Returns `list[UploadedFile]`

  **WHY Each Reference Matters**:
  - Lines 441-520: This is the exact code being modified; understand the tempfile pattern and session_state reset
  - Demo loading (449-475): Shows the pattern for resetting all session_state keys — multi-file must do the same
  - `merge_parse_results`: The new API to call after parsing all files individually

  **Acceptance Criteria**:

  ```bash
  # Verify import works
  python3 -c "
  import ast, sys
  with open('rnaseq_analysis_platform.py') as f:
      src = f.read()
  assert 'accept_multiple_files' in src, 'accept_multiple_files not found'
  assert 'merge_parse_results' in src, 'merge_parse_results import not found'
  print('Upload UI wiring check PASSED')
  "
  # Assert: prints "Upload UI wiring check PASSED"
  
  # Regression test
  pytest tests/ -v --tb=short 2>&1 | tail -5
  # Assert: 142+ passed
  ```

  **Commit**: YES
  - Message: `feat(upload): support multi-file upload with automatic merge`
  - Files: `rnaseq_analysis_platform.py`
  - Pre-commit: `pytest tests/ -v --tb=short`

---

- [ ] 3. Expand HELP_TEXTS with educational tab explanations

  **What to do**:
  - In `rnaseq_analysis_platform.py`, expand the `HELP_TEXTS` dict (line 251) to include all 15 tabs
  - Add these missing keys with rich educational content:
    - `summary` — What the summary metrics mean, what's normal
    - `de_results` — How to read a DE results table, what padj/LFC mean, sorting strategies
    - `interpretation` — What automated interpretations are, how they're generated, limitations
    - `qc_advanced` — Why QC matters, what outliers mean, batch effects
    - `gene_panels` — What curated gene panels are, how to use dermatology panels
    - `export` — What formats are available, when to use each
  - Enhance existing 7 unused keys (volcano, ma_plot, pca, heatmap, enrichment, gsea, deconvolution) with richer content:
    - Each should include: **Purpose**, **How to Read This**, **What to Look For**, **Common Pitfalls**
    - Include a log2FC reference table in volcano/MA help: `|LFC| = 1 → 2-fold, = 2 → 4-fold, = 3 → 8-fold`
  - Format: Multi-line strings with markdown formatting (headers, bullets, bold) since they'll be rendered with `st.markdown()`
  - Target audience: molecular biologist doing RNA-seq for the first time

  **Must NOT do**:
  - Do NOT wire the help texts into tabs yet (Task 5 does that)
  - Do NOT modify any tab rendering code
  - Do NOT add non-English content

  **Recommended Agent Profile**:
  - **Category**: `writing`
  - **Skills**: []
  - **Reason**: Content writing task — requires biological domain knowledge and clear explanatory writing

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 4)
  - **Blocks**: Task 5
  - **Blocked By**: None

  **References**:

  **Pattern References**:
  - `rnaseq_analysis_platform.py:251-265` — Current HELP_TEXTS dict with 13 keys (style reference)
  - `rnaseq_analysis_platform.py:1165` — `st.caption(HELP_TEXTS["gene_expression"])` — How help is currently displayed (simple caption)

  **Documentation References**:
  - README.md — Feature descriptions for each analysis type
  - `interpretation_engine.py:42-50` — DE results column descriptions (gene, log2FoldChange, padj, pvalue, baseMean)

  **WHY Each Reference Matters**:
  - Current HELP_TEXTS: Shows existing style/length conventions (currently brief tooltips — we're expanding to educational blocks)
  - README: Contains accurate feature descriptions to base educational content on
  - DE column descriptions: Need accurate technical details for the DE Results tab explanation

  **Acceptance Criteria**:

  ```bash
  # Verify all 19 keys exist (13 existing + 6 new)
  python3 -c "
  import ast
  with open('rnaseq_analysis_platform.py') as f:
      src = f.read()
  tree = ast.parse(src)
  for node in ast.walk(tree):
      if isinstance(node, ast.Assign):
          for target in node.targets:
              if hasattr(target, 'id') and target.id == 'HELP_TEXTS':
                  if isinstance(node.value, ast.Dict):
                      keys = [k.value for k in node.value.keys if isinstance(k, ast.Constant)]
                      required = ['summary', 'de_results', 'interpretation', 'qc_advanced', 'gene_panels', 'export',
                                  'volcano', 'ma_plot', 'pca', 'heatmap', 'enrichment', 'gsea', 'deconvolution',
                                  'gene_expression', 'venn', 'padj_threshold', 'lfc_threshold', 'upload', 'metadata']
                      missing = [k for k in required if k not in keys]
                      assert not missing, f'Missing keys: {missing}'
                      print(f'All {len(required)} HELP_TEXTS keys present. PASSED')
  "
  # Assert: prints "All 19 HELP_TEXTS keys present. PASSED"
  
  # Verify educational content is substantial (not just tooltips)
  python3 -c "
  with open('rnaseq_analysis_platform.py') as f:
      src = f.read()
  # Check that new keys have multi-line content (>100 chars)
  for key in ['summary', 'de_results', 'interpretation', 'qc_advanced', 'gene_panels', 'export']:
      idx = src.find(f'\"{key}\"')
      assert idx > 0, f'Key {key} not found'
  print('Content length check PASSED')
  "
  ```

  **Commit**: YES
  - Message: `feat(help): add educational explanations for all 15 analysis tabs`
  - Files: `rnaseq_analysis_platform.py`
  - Pre-commit: `python3 -c "import rnaseq_analysis_platform; print('OK')"`

---

- [ ] 4. Enhance InterpretationEngine with deep context-aware methods

  **What to do**:
  - Modify `interpretation_engine.py` to add these new methods:

  1. **`interpret_top_genes(de_results_df, contrast_name, top_n=10)`**:
     - Identify top 10 up/down genes by |LFC| among significant
     - Name them explicitly: "Top upregulated: COL1A1 (LFC=3.2), ELN (LFC=2.8)..."
     - Flag known dermatology-relevant genes (from gene_panels.yaml categories)

  2. **`interpret_volcano_pattern(de_results_df, contrast_name)`**:
     - Classify response pattern: symmetric, skewed-up, skewed-down, bimodal
     - Note outlier genes with extreme LFC (>5) or extreme significance (-log10p > 50)

  3. **`interpret_pca(pca_coordinates, sample_metadata)`**:
     - Assess clustering quality: do conditions separate on PC1/PC2?
     - Flag potential outlier samples (distance from group centroid > 2 SD)
     - Report variance explained by PC1/PC2 if available

  4. **`interpret_gsea_results(gsea_df)`**:
     - Analyze NES (Normalized Enrichment Score) for directionality
     - Identify top activated vs suppressed pathways
     - Compare with ORA results if available (overlap of findings)

  5. **`interpret_deconvolution(deconv_results)`**:
     - Report dominant cell types
     - Flag unexpected cell type proportions
     - Suggest implications for interpretation

  6. **`generate_misinterpretation_warnings(de_results_df, enrichment_df=None)`**:
     - Generate warnings for 9 common pitfalls:
       1. "More significant ≠ more biologically important" (if many genes at padj~0.049)
       2. "Absence of evidence ≠ evidence of absence" (if many genes near threshold)
       3. "Fold change alone is insufficient" (if high LFC but low expression genes)
       4. "Correlation ≠ causation" (always include)
       5. "Enrichment bias" (if DE list is very large, enrichment loses specificity)
       6. "Multiple testing burden" (if >10k genes tested)
       7. "Replicate concern" (if <3 samples per condition detected)
       8. "Batch effect warning" (if PCA PC1 separates by non-biological factor)
       9. "Platform/normalization caveat" (if normalized data used instead of raw counts)
     - Only trigger warnings that are relevant to the actual data

  - Enhance existing `interpret_enrichment()`:
    - Add top 3 enriched terms by name in the description
    - Add overlap count: "X of your DE genes appear in this pathway"
    - Add skin biology keywords: `keratin`, `melanin`, `sebaceous`, `wound`, `filaggrin`, `barrier`

  **Must NOT do**:
  - Do NOT change the `Interpretation` dataclass structure (keep same fields)
  - Do NOT modify `generate_methods_text()` 
  - Do NOT import from `de_analysis.py` or `pathway_enrichment.py`
  - Do NOT add external API calls

  **Recommended Agent Profile**:
  - **Category**: `ultrabrain`
  - **Skills**: []
  - **Reason**: Complex biological domain logic, multiple interacting methods, needs careful threshold tuning

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 3)
  - **Blocks**: Task 5
  - **Blocked By**: None

  **References**:

  **Pattern References**:
  - `interpretation_engine.py:35-177` — `interpret_de_results()` full implementation (pattern for new methods: filter data → compute metrics → create Interpretation objects with evidence dict)
  - `interpretation_engine.py:179-308` — `interpret_enrichment()` full implementation (keyword matching pattern to extend)
  - `interpretation_engine.py:13-21` — `Interpretation` dataclass (category, level, title, description, evidence, recommendations)

  **API/Type References**:
  - `de_results_df` columns: `gene`, `log2FoldChange`, `padj`, `pvalue`, `baseMean`
  - `pca_coordinates`: DataFrame with PC1, PC2 columns, sample names as index
  - `gsea_df`: DataFrame from gseapy with `Term`, `NES`, `NOM p-val`, `FDR q-val` columns
  - `deconv_results`: DataFrame with cell types as columns, samples as rows, proportions as values

  **Documentation References**:
  - `config/gene_panels.yaml` — Dermatology gene panel definitions (gene lists for skin biology detection)

  **WHY Each Reference Matters**:
  - `interpret_de_results()`: This is your template — every new method should follow the same pattern of creating `Interpretation` objects with structured evidence
  - `Interpretation` dataclass: You MUST return objects matching this exact structure
  - DE results columns: Exact column names needed for DataFrame operations
  - `gene_panels.yaml`: Gene lists for detecting skin-biology-relevant genes in top DE results

  **Acceptance Criteria**:

  ```bash
  # Verify all new methods exist
  python3 -c "
  from interpretation_engine import InterpretationEngine
  engine = InterpretationEngine()
  methods = ['interpret_top_genes', 'interpret_volcano_pattern', 'interpret_pca', 
             'interpret_gsea_results', 'interpret_deconvolution', 'generate_misinterpretation_warnings']
  for m in methods:
      assert hasattr(engine, m), f'Missing method: {m}'
  print('All 6 new methods present. PASSED')
  "

  # Verify interpret_top_genes produces named genes
  python3 -c "
  import pandas as pd
  from interpretation_engine import InterpretationEngine
  df = pd.DataFrame({
      'gene': ['COL1A1', 'IL6', 'TYR', 'FLG', 'GAPDH'],
      'log2FoldChange': [3.2, -2.1, 1.5, -0.3, 0.1],
      'padj': [0.001, 0.01, 0.03, 0.5, 0.9],
      'pvalue': [0.0001, 0.001, 0.005, 0.1, 0.5],
      'baseMean': [1500, 450, 120, 800, 5000]
  })
  engine = InterpretationEngine()
  results = engine.interpret_top_genes(df, 'Test vs Control')
  assert len(results) > 0, 'No interpretations generated'
  assert 'COL1A1' in results[0].description, 'Top gene not named in description'
  print('Top genes interpretation PASSED')
  "

  # Verify misinterpretation warnings
  python3 -c "
  import pandas as pd
  from interpretation_engine import InterpretationEngine
  df = pd.DataFrame({
      'gene': [f'Gene{i}' for i in range(100)],
      'log2FoldChange': [1.5]*50 + [-1.5]*50,
      'padj': [0.01]*100,
      'pvalue': [0.001]*100,
      'baseMean': [500]*100
  })
  engine = InterpretationEngine()
  warnings = engine.generate_misinterpretation_warnings(df)
  assert len(warnings) > 0, 'No warnings generated'
  assert all(w.category == 'de' or w.category == 'warning' for w in warnings)
  print('Misinterpretation warnings PASSED')
  "

  # Regression
  pytest tests/ -v --tb=short 2>&1 | tail -5
  # Assert: 142+ passed
  ```

  **Commit**: YES
  - Message: `feat(interpret): add 6 new context-aware interpretation methods with misinterpretation warnings`
  - Files: `interpretation_engine.py`
  - Pre-commit: `pytest tests/ -v --tb=short`

---

- [ ] 5. Wire tab explanations and enhanced interpretations into UI

  **What to do**:
  - In `rnaseq_analysis_platform.py`, add educational explanation blocks to ALL 15 tabs:
    1. At the top of each `with tab_X:` block, add an expandable explanation:
       ```python
       with st.expander("ℹ️ How to read this tab", expanded=False):
           st.markdown(HELP_TEXTS["key_name"])
       ```
    2. Wire all 15 tabs: summary, de_results, volcano, ma_plot, heatmap, pca, enrichment, gsea, interpretation, deconvolution, qc_advanced, gene_panels, gene_expression, venn, export
    3. Replace the 2 existing `st.caption()` calls (lines 1165, 1182) with the new expander pattern for consistency
  
  - Wire enhanced InterpretationEngine into the interpretation tab:
    1. After existing `interpret_de_results()` call, add calls to new methods:
       - `interpret_top_genes()`
       - `interpret_volcano_pattern()`
       - `generate_misinterpretation_warnings()`
    2. In the interpretation tab rendering, also call context-specific methods when data is available:
       - If `pca_coordinates` in session_state → `interpret_pca()`
       - If `gsea_results` available → `interpret_gsea_results()`
       - If `deconvolution_results` available → `interpret_deconvolution()`
    3. Group interpretations by category with section headers in the UI

  - Find where `InterpretationEngine` is called to generate interpretations (likely in the analysis run flow) and add the new method calls there

  **Must NOT do**:
  - Do NOT change tab order or tab names
  - Do NOT modify the analysis pipeline itself
  - Do NOT add new tabs

  **Recommended Agent Profile**:
  - **Category**: `visual-engineering`
  - **Skills**: []
  - **Reason**: Streamlit UI integration, needs understanding of session_state flow and tab structure

  **Parallelization**:
  - **Can Run In Parallel**: YES (with Task 2)
  - **Parallel Group**: Wave 2
  - **Blocks**: Task 6
  - **Blocked By**: Tasks 3, 4

  **References**:

  **Pattern References**:
  - `rnaseq_analysis_platform.py:904-909` — Tab creation with `st.tabs()` (all 15 tab variable names)
  - `rnaseq_analysis_platform.py:911-915` — `tab_summary` block start (pattern for where to insert expander)
  - `rnaseq_analysis_platform.py:1041-1087` — Current `tab_interpret` rendering (interpretation display pattern)
  - `rnaseq_analysis_platform.py:1165` — Existing `st.caption(HELP_TEXTS["gene_expression"])` (to replace)
  - `rnaseq_analysis_platform.py:1182` — Existing `st.caption(HELP_TEXTS["venn"])` (to replace)

  **API/Type References**:
  - `InterpretationEngine.interpret_top_genes(df, name)` — New from Task 4
  - `InterpretationEngine.interpret_volcano_pattern(df, name)` — New from Task 4
  - `InterpretationEngine.interpret_pca(coords, metadata)` — New from Task 4
  - `InterpretationEngine.interpret_gsea_results(gsea_df)` — New from Task 4
  - `InterpretationEngine.interpret_deconvolution(deconv)` — New from Task 4
  - `InterpretationEngine.generate_misinterpretation_warnings(df, enrichment_df)` — New from Task 4

  **WHY Each Reference Matters**:
  - Tab creation line 904: Maps tab variable names to their string labels — need exact variable names
  - Tab_interpret rendering: This is where interpretations are displayed; need to add new interpretation calls and grouping
  - Existing caption calls: Must be replaced with consistent expander pattern
  - New method signatures: Need exact API to call correctly

  **Acceptance Criteria**:

  ```bash
  # Verify all 15 tabs have help expanders
  python3 -c "
  with open('rnaseq_analysis_platform.py') as f:
      src = f.read()
  # Count 'How to read this tab' occurrences (should be 15)
  count = src.count('How to read this tab')
  assert count >= 15, f'Only {count}/15 tabs have explanations'
  print(f'Found {count} tab explanations. PASSED')
  "

  # Verify old st.caption calls replaced
  python3 -c "
  with open('rnaseq_analysis_platform.py') as f:
      src = f.read()
  # Should NOT have st.caption(HELP_TEXTS anymore
  import re
  old_pattern = re.findall(r'st\.caption\(HELP_TEXTS', src)
  assert len(old_pattern) == 0, f'Found {len(old_pattern)} old st.caption(HELP_TEXTS calls'
  print('Old caption calls replaced. PASSED')
  "

  # Verify enhanced interpretation wiring
  python3 -c "
  with open('rnaseq_analysis_platform.py') as f:
      src = f.read()
  for method in ['interpret_top_genes', 'interpret_volcano_pattern', 'generate_misinterpretation_warnings']:
      assert method in src, f'{method} not wired in UI'
  print('Enhanced interpretations wired. PASSED')
  "

  # Regression
  pytest tests/ -v --tb=short 2>&1 | tail -5
  # Assert: 142+ passed
  ```

  **Commit**: YES
  - Message: `feat(ui): wire tab explanations and enhanced interpretations across all 15 tabs`
  - Files: `rnaseq_analysis_platform.py`
  - Pre-commit: `pytest tests/ -v --tb=short`

---

- [ ] 6. Integration testing and regression verification

  **What to do**:
  - Run full test suite: `pytest tests/ -v`
  - Verify no regressions (142+ tests pass, same 3 pre-existing PyDESeq2 failures)
  - Write a quick smoke test for the new merger:
    - Create `tests/test_multi_file_merger.py` with:
      1. Test single file passthrough
      2. Test two-file merge with overlapping genes
      3. Test duplicate sample name handling
      4. Test DataType mismatch rejection
      5. Test empty file handling
  - Verify Streamlit app launches without errors:
    ```bash
    timeout 10 streamlit run rnaseq_analysis_platform.py --server.headless true 2>&1 | head -20
    ```
  - Verify all imports resolve correctly

  **Must NOT do**:
  - Do NOT modify existing test files
  - Do NOT change test configuration

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: []
  - **Reason**: Test writing and verification, straightforward pytest patterns

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (sequential, final)
  - **Blocks**: None (final task)
  - **Blocked By**: Tasks 2, 5

  **References**:

  **Pattern References**:
  - `tests/test_integration.py` — Existing integration test pattern (fixtures, assertions)
  - `conftest.py` — Test fixtures and mocks (especially `mock_streamlit`)
  - `tests/test_multiformat_parser.py` — Parser test patterns

  **API/Type References**:
  - `multi_file_merger.merge_parse_results(results, filenames)` — From Task 1
  - `multi_file_merger.validate_merge_compatibility(results)` — From Task 1

  **WHY Each Reference Matters**:
  - Existing test patterns: Follow same fixture/assertion style for consistency
  - conftest.py: May need to understand mock patterns but should NOT modify them

  **Acceptance Criteria**:

  ```bash
  # Full test suite
  pytest tests/ -v 2>&1 | tail -10
  # Assert: 142+ passed (plus new merger tests), ≤3 xfail/skip for PyDESeq2

  # New merger tests specifically
  pytest tests/test_multi_file_merger.py -v
  # Assert: 5+ tests pass

  # App launch check
  timeout 10 streamlit run rnaseq_analysis_platform.py --server.headless true 2>&1 | grep -E "You can now view|error"
  # Assert: "You can now view" appears, no "error"
  ```

  **Commit**: YES
  - Message: `test(merger): add integration tests for multi-file merge functionality`
  - Files: `tests/test_multi_file_merger.py`
  - Pre-commit: `pytest tests/ -v --tb=short`

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `feat(merger): add multi-file merge utility for combining count matrices` | `multi_file_merger.py` | import check |
| 2 | `feat(upload): support multi-file upload with automatic merge` | `rnaseq_analysis_platform.py` | `pytest tests/ -v` |
| 3 | `feat(help): add educational explanations for all 15 analysis tabs` | `rnaseq_analysis_platform.py` | import check |
| 4 | `feat(interpret): add 6 new context-aware interpretation methods` | `interpretation_engine.py` | `pytest tests/ -v` |
| 5 | `feat(ui): wire tab explanations and enhanced interpretations` | `rnaseq_analysis_platform.py` | `pytest tests/ -v` |
| 6 | `test(merger): add integration tests for multi-file merge` | `tests/test_multi_file_merger.py` | `pytest tests/ -v` |

---

## Success Criteria

### Verification Commands
```bash
# Full regression
pytest tests/ -v                    # Expected: 147+ passed (142 existing + 5+ new)

# Import check
python3 -c "from multi_file_merger import merge_parse_results; from interpretation_engine import InterpretationEngine; print('All imports OK')"

# App launch
timeout 10 streamlit run rnaseq_analysis_platform.py --server.headless true 2>&1 | head -5
# Expected: "You can now view your Streamlit app"
```

### Final Checklist
- [ ] Multi-file upload works with 2+ CSV/Excel files
- [ ] Single-file upload still works identically to before
- [ ] All 15 tabs have educational explanation expanders
- [ ] InterpretationEngine generates 10+ insight types
- [ ] Misinterpretation warnings appear when relevant
- [ ] Top DE genes named explicitly in interpretations
- [ ] 142+ existing tests pass (no regressions)
- [ ] New merger tests pass
- [ ] App launches without errors
