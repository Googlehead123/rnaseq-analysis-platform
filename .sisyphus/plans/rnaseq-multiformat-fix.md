# RNA-seq Multi-Format Parsing Fix

## Context

### Original Request
User gets "Analysis failed: unsupported operand type(s) for +: 'NoneType' and 'int'" when uploading RNA-seq reference data files (data3_*.xlsx). The root cause is that these files contain multiple data types (DE results + raw counts + normalized values) but the parser only extracts one, leaving `expression_df=None` for PRE_ANALYZED detection while visualization code assumes it's populated.

### Interview Summary
**Key Discussions**:
- **Data Extraction Strategy**: User chose Smart Multi-Extract - extract ALL three data types into separate fields
- **Column Mapping**: Auto-detect with Confirmation - regex patterns detect non-standard columns, show confirmation message
- **Sample Detection**: Auto-detect Conditions from column names (e.g., `231222_none_3d_Read_Count` → condition: `none`)
- **Test Strategy**: TDD - write failing tests first with reference files as fixtures

**Research Findings**:
- Reference file structure: 42 columns with metadata (10), DE results (6), raw counts (6), FPKM (6), TPM (6)
- Column naming pattern: `{comparison}.{stat}` (e.g., `Bt10U/none.fc`, `Bt10U/none.bh.pval`)
- Sample column pattern: `{date}_{condition}_{timepoint}_{metric}` (e.g., `231222_none_3d_Read_Count`)
- 8 critical None-check gaps identified across codebase
- 2 silent failure points (lines 509, 288)
- Existing test infrastructure: pytest with fixtures in `tests/`

### Guardrails Applied
- **Backward Compatibility**: Existing simple CSV/Excel files MUST continue to work exactly as before
- **No PyDESeq2 Changes**: DE analysis engine remains untouched
- **Single Comparison Per File**: Initially support ONE comparison per multi-format file (future enhancement for multiple)
- **Graceful Degradation**: If pattern detection fails, fall back to current behavior with clear error message

---

## Work Objectives

### Core Objective
Fix the NoneType crash and enable the RNA-seq platform to intelligently parse multi-format Excel files that contain DE results, raw counts, and normalized values simultaneously, extracting each data type for appropriate use in visualizations.

### Concrete Deliverables
- Fixed crash: No more `NoneType + int` error on reference file upload
- Extended `ParseResult` dataclass with `normalized_df` field
- Pattern-based column detection for non-standard names
- Sample/condition auto-detection from column name patterns
- Comprehensive None guards across all visualization and export paths
- Test suite covering multi-format parsing scenarios

### Definition of Done
- [ ] Upload `data3_Bt10U_vs_none_fc2_&_raw.p.xlsx` → no crash
- [ ] Volcano plot renders with DE results from `.fc` and `.bh.pval` columns
- [ ] Heatmap renders with normalized values from `*_FPKM` or `*_TPM` columns
- [ ] PCA plot renders with samples auto-detected from column names
- [ ] All existing tests pass: `pytest tests/ -v`
- [ ] New tests cover multi-format scenarios

### Must Have
- P0: Immediate crash fix (None guard at line 474)
- P1: Multi-data-type extraction from single file
- P1: Pattern-based column name detection
- P1: Sample/condition auto-detection
- P1: Comprehensive None checks

### Must NOT Have (Guardrails)
- NO changes to `de_analysis.py` or PyDESeq2 integration
- NO breaking changes to existing CSV/TSV parsing
- NO support for multiple comparisons per file (scope for future)
- NO new dependencies (use existing pandas, numpy, re)
- NO over-engineering: patterns must be specific to observed data, not hypothetical formats

---

## Verification Strategy (MANDATORY)

### Test Decision
- **Infrastructure exists**: YES (pytest in `tests/`)
- **User wants tests**: YES (TDD)
- **Framework**: pytest

### TDD Workflow
Each TODO follows RED-GREEN-REFACTOR:
1. **RED**: Write failing test first
2. **GREEN**: Implement minimum code to pass
3. **REFACTOR**: Clean up while keeping tests green

### Test Files
- `tests/test_multiformat_parser.py` - New file for multi-format parsing tests
- `tests/test_none_guards.py` - New file for None-check validation tests
- Use `Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx` as test fixture

---

## Task Flow

```
[P0: Immediate Fix]
    Task 1 (None guard)
           ↓
[P1: Parser Enhancement - Wave 1]
    Task 2 (ParseResult extension)
           ↓
    Task 3 (Column patterns) ──────┐
    Task 4 (Sample detection) ─────┤ (parallel)
           ↓                       ↓
    Task 5 (Multi-format parser) ←─┘
           ↓
[P1: Integration - Wave 2]
    Task 6 (Platform integration)
           ↓
[P1: Comprehensive None Guards - Wave 3]
    Task 7 (Visualization guards) ─┐
    Task 8 (Export guards) ────────┤ (parallel)
    Task 9 (Gene panel guards) ────┘
           ↓
[P2: UX Enhancement - Wave 4]
    Task 10 (Confirmation messages)
    Task 11 (Error message improvements)
           ↓
[P3: Cleanup - Wave 5]
    Task 12 (Silent failure fixes)
    Task 13 (Final integration test)
```

## Parallelization

| Group | Tasks | Reason |
|-------|-------|--------|
| A | 3, 4 | Independent pattern detection modules |
| B | 7, 8, 9 | Independent None guard additions per module |

| Task | Depends On | Reason |
|------|------------|--------|
| 5 | 2, 3, 4 | Multi-format parser needs extended ParseResult and detection utilities |
| 6 | 5 | Platform integration needs complete parser |
| 7, 8, 9 | 6 | Guards need integration to verify behavior |
| 10, 11 | 6 | UX needs working integration |
| 13 | 1-12 | Final verification needs all fixes |

---

## TODOs

### P0: IMMEDIATE FIX (Unblocks User)

- [ ] 1. Add None guard at line 474 crash location

  **What to do**:
  - Write test that loads PRE_ANALYZED data and attempts visualization path
  - Test should fail with current code (NoneType + int error)
  - Add None check before `np.log2(result.expression_df + 1)`
  - For PRE_ANALYZED, skip heatmap/PCA generation (they require expression data)
  - Add user-facing info message explaining why heatmap/PCA unavailable

  **Must NOT do**:
  - Do NOT change the overall flow structure
  - Do NOT add complex fallback logic yet (that's Task 5)

  **Parallelizable**: NO (first task, unblocks everything)

  **References**:

  **Pattern References**:
  - `rnaseq_analysis_platform.py:216-223` - Existing None check pattern for `expression_df` in data preview
  - `rnaseq_analysis_platform.py:174-186` - Pattern for handling None DataFrames in comparison

  **Error Location**:
  - `rnaseq_analysis_platform.py:471-474` - THE CRASH LOCATION:
    ```python
    elif result.data_type == DataType.NORMALIZED:
        import numpy as np
        norm_counts = np.log2(result.expression_df + 1)  # CRASHES HERE
    ```

  **Context**:
  - `rnaseq_analysis_platform.py:464-476` - Full visualization block that needs guarding
  - `rnaseq_parser.py:510-522` - Why expression_df is None for PRE_ANALYZED

  **Test References**:
  - `tests/test_integration.py:67-85` - Existing pattern for testing parsed results

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test file created: `tests/test_none_guards.py`
  - [ ] Test: `test_preanalyzed_data_no_crash_on_visualization`
  - [ ] `pytest tests/test_none_guards.py -v` → PASS (1 test)

  **Manual Verification**:
  - [ ] Using Streamlit app:
    - Upload: `Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx`
    - Click: "Run Analysis"
    - Verify: No crash, info message displayed if heatmap/PCA unavailable
    - Screenshot saved to `.sisyphus/evidence/task-1-no-crash.png`

  **Commit**: YES
  - Message: `fix(platform): guard None expression_df in visualization path`
  - Files: `rnaseq_analysis_platform.py`, `tests/test_none_guards.py`
  - Pre-commit: `pytest tests/test_none_guards.py -v`

---

### P1: PARSER ENHANCEMENT

- [ ] 2. Extend ParseResult dataclass to support multiple data types

  **What to do**:
  - Write test that checks ParseResult can hold expression_df, normalized_df, and de_results_df simultaneously
  - Add `normalized_df: Optional[pd.DataFrame]` field to ParseResult
  - Add `data_types_detected: List[DataType]` field to track what was found
  - Update docstring to document new field population rules
  - Ensure backward compatibility: existing code paths still work

  **Must NOT do**:
  - Do NOT change any parsing logic yet (just the dataclass)
  - Do NOT break existing ParseResult usage patterns

  **Parallelizable**: NO (foundational change, blocks 3, 4, 5)

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py:85-128` - Current ParseResult definition with field documentation

  **Type References**:
  - `rnaseq_parser.py:55-60` - DataType enum definition

  **Test References**:
  - `tests/test_integration.py:45-65` - How ParseResult is currently tested

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_parseresult_supports_multiple_datatypes`
  - [ ] Test: `test_parseresult_backward_compatible`
  - [ ] `pytest tests/test_multiformat_parser.py::test_parseresult* -v` → PASS

  **Commit**: YES
  - Message: `feat(parser): extend ParseResult to support multi-format files`
  - Files: `rnaseq_parser.py`, `tests/test_multiformat_parser.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 3. Implement pattern-based column name detection for DE results

  **What to do**:
  - Write tests for detecting `{comparison}.fc` → `log2FoldChange` mapping
  - Write tests for detecting `{comparison}.bh.pval` → `padj` mapping
  - Create `DE_COLUMN_PATTERNS` dict with regex patterns:
    ```python
    DE_COLUMN_PATTERNS = {
        "log2FoldChange": [r"^.+\.fc$", r"^.+\.logFC$", r"^.+\.log2FC$"],
        "padj": [r"^.+\.bh\.pval$", r"^.+\.adj\.pval$", r"^.+\.FDR$", r"^.+\.qvalue$"],
        "pvalue": [r"^.+\.raw\.pval$", r"^.+\.pval$", r"^.+\.PValue$"],
        "baseMean": [r"^.+\.baseMean$", r"^.+\.AveExpr$"],
    }
    ```
  - Implement `detect_de_columns(df: DataFrame) -> Dict[str, str]` function
  - Return mapping of canonical name → actual column name
  - Extract comparison name from pattern (e.g., "Bt10U/none" from "Bt10U/none.fc")

  **Must NOT do**:
  - Do NOT modify existing COLUMN_ALIASES (keep for backward compat)
  - Do NOT make detection mandatory (graceful fallback to current behavior)

  **Parallelizable**: YES (with Task 4)

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py:46-52` - Existing COLUMN_ALIASES dict
  - `rnaseq_parser.py:472-489` - Existing normalize_de_columns function

  **Data References**:
  - Reference file columns: `Bt10U/none.fc`, `Bt10U/none.baseMean`, `Bt10U/none.lfcSE`, `Bt10U/none.stat`, `Bt10U/none.raw.pval`, `Bt10U/none.bh.pval`

  **External References**:
  - Python re module: https://docs.python.org/3/library/re.html

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_detect_de_columns_pattern_fc`
  - [ ] Test: `test_detect_de_columns_pattern_pval`
  - [ ] Test: `test_detect_de_columns_extracts_comparison_name`
  - [ ] Test: `test_detect_de_columns_fallback_to_standard`
  - [ ] `pytest tests/test_multiformat_parser.py::test_detect_de* -v` → PASS

  **Commit**: YES
  - Message: `feat(parser): add pattern-based DE column detection`
  - Files: `rnaseq_parser.py`, `tests/test_multiformat_parser.py`
  - Pre-commit: `pytest tests/test_multiformat_parser.py -v`

---

- [ ] 4. Implement sample/condition auto-detection from column names

  **What to do**:
  - Write tests for detecting count columns: `{id}_{condition}_{timepoint}_Read_Count`
  - Write tests for detecting normalized columns: `{id}_{condition}_{timepoint}_FPKM`, `*_TPM`
  - Create pattern constants:
    ```python
    COUNT_COLUMN_PATTERN = r"^(.+)_(.+)_(.+)_Read_Count$"  # (date, condition, timepoint)
    FPKM_COLUMN_PATTERN = r"^(.+)_(.+)_(.+)_FPKM$"
    TPM_COLUMN_PATTERN = r"^(.+)_(.+)_(.+)_TPM$"
    ```
  - Implement `detect_sample_columns(df: DataFrame) -> SampleColumnInfo`:
    ```python
    @dataclass
    class SampleColumnInfo:
        count_columns: List[str]
        fpkm_columns: List[str]
        tpm_columns: List[str]
        conditions_detected: Set[str]
        sample_to_condition: Dict[str, str]  # column_name → condition
    ```
  - Extract conditions from column names (e.g., "none", "Bt10U")

  **Must NOT do**:
  - Do NOT make auto-detection mandatory (user can still manually assign)
  - Do NOT assume specific date formats (be flexible)

  **Parallelizable**: YES (with Task 3)

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py:131-142` - Existing `looks_like_sample_names` heuristic

  **Data References**:
  - Reference file columns: `231222_none_3d_Read_Count`, `240203_none_3d_Read_Count`, `240726_none_3d_Read_Count`, `231222_Bt10U_3d_Read_Count`, etc.

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_detect_sample_columns_read_count`
  - [ ] Test: `test_detect_sample_columns_fpkm`
  - [ ] Test: `test_detect_sample_columns_tpm`
  - [ ] Test: `test_detect_conditions_from_columns`
  - [ ] Test: `test_sample_to_condition_mapping`
  - [ ] `pytest tests/test_multiformat_parser.py::test_detect_sample* -v` → PASS

  **Commit**: YES
  - Message: `feat(parser): add sample/condition auto-detection from column names`
  - Files: `rnaseq_parser.py`, `tests/test_multiformat_parser.py`
  - Pre-commit: `pytest tests/test_multiformat_parser.py -v`

---

- [ ] 5. Implement multi-format file parsing in RNASeqParser

  **What to do**:
  - Write integration test that parses `data3_Bt10U_vs_none_fc2_&_raw.p.xlsx`
  - Test should verify: de_results_df populated, expression_df populated (from counts), normalized_df populated (from FPKM/TPM)
  - Add `parse_multiformat` method to RNASeqParser:
    ```python
    def parse_multiformat(self, file_path: str) -> ParseResult:
        """Parse multi-format file with DE results + counts + normalized."""
    ```
  - Detection logic:
    1. Load Excel file
    2. Call `detect_de_columns()` - if found, extract DE results
    3. Call `detect_sample_columns()` - if count/FPKM/TPM columns found, extract expression matrices
    4. Populate ParseResult with ALL available data types
    5. Set `data_types_detected` list
  - Modify main `parse()` method to try `parse_multiformat` when:
    - File is Excel AND
    - Has both DE-like columns AND expression-like columns
  - Build expression_df from count columns (genes × samples → samples × genes)
  - Build normalized_df from FPKM or TPM columns (prefer TPM if both exist)
  - Build de_results_df with canonical column names

  **Must NOT do**:
  - Do NOT break existing simple file parsing
  - Do NOT modify parse_csv (Excel-only feature)
  - Do NOT require all three data types (any combination is valid)

  **Parallelizable**: NO (depends on 2, 3, 4)

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py:495-572` - Current `parse()` method structure
  - `rnaseq_parser.py:428-469` - Current `parse_excel()` function
  - `rnaseq_parser.py:296-359` - `convert_to_canonical_shape()` for orientation handling

  **Data References**:
  - Reference file: `Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx`
  - Expected output: de_results_df with ~500 genes, expression_df with 6 samples × ~500 genes

  **Test References**:
  - `tests/test_integration.py:150-200` - Integration test patterns

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_parse_multiformat_extracts_de_results`
  - [ ] Test: `test_parse_multiformat_extracts_counts`
  - [ ] Test: `test_parse_multiformat_extracts_normalized`
  - [ ] Test: `test_parse_multiformat_all_datatypes_populated`
  - [ ] Test: `test_parse_multiformat_backward_compatible`
  - [ ] `pytest tests/test_multiformat_parser.py -v` → PASS (all tests)

  **Manual Verification**:
  - [ ] Python REPL:
    ```python
    from rnaseq_parser import RNASeqParser
    parser = RNASeqParser()
    result = parser.parse("Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx")
    print(f"DE results: {result.de_results_df.shape if result.de_results_df is not None else None}")
    print(f"Expression: {result.expression_df.shape if result.expression_df is not None else None}")
    print(f"Normalized: {result.normalized_df.shape if result.normalized_df is not None else None}")
    print(f"Data types: {result.data_types_detected}")
    ```
  - Expected output:
    ```
    DE results: (~500, 7)  # genes × DE columns
    Expression: (6, ~500)  # samples × genes
    Normalized: (6, ~500)  # samples × genes
    Data types: [DataType.PRE_ANALYZED, DataType.RAW_COUNTS, DataType.NORMALIZED]
    ```

  **Commit**: YES
  - Message: `feat(parser): implement multi-format Excel parsing`
  - Files: `rnaseq_parser.py`, `tests/test_multiformat_parser.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 6. Integrate multi-format parser with platform UI

  **What to do**:
  - Modify `rnaseq_analysis_platform.py` to use new multi-format parsing
  - When multi-format file detected:
    1. Use `de_results_df` for volcano plots
    2. Use `normalized_df` (or `expression_df` if normalized unavailable) for heatmap/PCA
    3. Auto-populate sample metadata from `sample_to_condition` mapping
    4. Skip metadata assignment UI for auto-detected conditions
    5. Show confirmation message: "Detected X conditions: {list}. DE results: {comparison}"
  - Update visualization path (lines 464-518):
    - Check for `normalized_df` first, then `expression_df`
    - Use appropriate source for each visualization type
  - Update session state to track multiple data sources

  **Must NOT do**:
  - Do NOT remove manual metadata assignment (keep as override option)
  - Do NOT change DE analysis path for RAW_COUNTS files
  - Do NOT modify export logic yet (Task 8)

  **Parallelizable**: NO (depends on 5)

  **References**:

  **Pattern References**:
  - `rnaseq_analysis_platform.py:152-208` - File upload and parsing flow
  - `rnaseq_analysis_platform.py:225-333` - Metadata assignment section
  - `rnaseq_analysis_platform.py:384-407` - PRE_ANALYZED handling

  **UI References**:
  - `rnaseq_analysis_platform.py:195-198` - Success message and warnings pattern
  - `rnaseq_analysis_platform.py:536-542` - Comparison selector UI pattern

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_platform_multiformat_file_upload`
  - [ ] Test: `test_platform_auto_metadata_population`
  - [ ] Test: `test_platform_visualization_uses_correct_source`
  - [ ] `pytest tests/test_integration.py -v` → PASS

  **Manual Verification (Playwright browser)**:
  - [ ] Start app: `streamlit run rnaseq_analysis_platform.py`
  - [ ] Navigate to: `http://localhost:8501`
  - [ ] Upload: `Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx`
  - [ ] Verify: Success message shows "Detected conditions: none, Bt10U"
  - [ ] Verify: Metadata section shows auto-populated conditions
  - [ ] Click: "Run Analysis"
  - [ ] Verify: Volcano plot appears with DE results
  - [ ] Verify: Heatmap appears with normalized data
  - [ ] Verify: PCA plot appears
  - [ ] Screenshot: `.sisyphus/evidence/task-6-full-workflow.png`

  **Commit**: YES
  - Message: `feat(platform): integrate multi-format parser with UI`
  - Files: `rnaseq_analysis_platform.py`, `tests/test_integration.py`
  - Pre-commit: `pytest tests/ -v`

---

### P1: COMPREHENSIVE NONE GUARDS

- [ ] 7. Add None guards to visualization functions

  **What to do**:
  - Write tests that pass None or empty DataFrames to each visualization function
  - Add input validation to `create_volcano_plot()`:
    - Check results_df is not None
    - Check required columns exist: gene, log2FoldChange, padj
    - Raise ValueError with helpful message if validation fails
  - Add input validation to `create_clustered_heatmap()`:
    - Check expression_df is not None and not empty
    - Check sample_conditions is not empty
    - Validate samples in expression_df match sample_conditions keys
  - Add input validation to `create_pca_plot()`:
    - Check expression_df is not None
    - Check sufficient samples for PCA (≥2)
    - Check sufficient features for requested components

  **Must NOT do**:
  - Do NOT change return types (still return go.Figure)
  - Do NOT add new optional parameters

  **Parallelizable**: YES (with Tasks 8, 9)

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py:362-379` - `validate_for_de()` as validation pattern
  - `gene_panels.py:100-122` - Validation pattern in `score_panel()`

  **Function Locations**:
  - `visualizations.py:17-71` - `create_volcano_plot()`
  - `visualizations.py:74-160` - `create_clustered_heatmap()`
  - `visualizations.py:163-204` - `create_pca_plot()`

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_volcano_plot_none_input_raises`
  - [ ] Test: `test_volcano_plot_missing_columns_raises`
  - [ ] Test: `test_heatmap_none_input_raises`
  - [ ] Test: `test_heatmap_empty_conditions_raises`
  - [ ] Test: `test_pca_none_input_raises`
  - [ ] Test: `test_pca_insufficient_samples_raises`
  - [ ] `pytest tests/test_none_guards.py::test_*_plot* -v` → PASS

  **Commit**: YES
  - Message: `fix(viz): add input validation to visualization functions`
  - Files: `visualizations.py`, `tests/test_none_guards.py`
  - Pre-commit: `pytest tests/test_none_guards.py -v`

---

- [ ] 8. Add None guards to export functions

  **What to do**:
  - Write tests that pass incomplete ExportData to export functions
  - Add validation to `export_excel()`:
    - Check ExportData has at least one data source
    - Validate filepath is writable (parent dir exists)
    - Handle None de_results gracefully (skip that sheet)
    - Handle None expression_matrix gracefully (skip that sheet)
  - Add validation to `export_pdf_report()`:
    - Same data validation as Excel
    - Check figures dict is not empty or handle gracefully
  - Update `ExportData` dataclass with validation method

  **Must NOT do**:
  - Do NOT change export file formats
  - Do NOT add new export types

  **Parallelizable**: YES (with Tasks 7, 9)

  **References**:

  **Pattern References**:
  - `export_engine.py:1-50` - ExportData and ExportEngine class structure
  - `rnaseq_analysis_platform.py:655-670` - How ExportData is constructed

  **Function Locations**:
  - `export_engine.py` - `export_excel()`, `export_pdf_report()` methods

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_export_excel_none_de_results_skips_sheet`
  - [ ] Test: `test_export_excel_none_expression_skips_sheet`
  - [ ] Test: `test_export_pdf_empty_figures_handled`
  - [ ] Test: `test_export_invalid_path_raises`
  - [ ] `pytest tests/test_none_guards.py::test_export* -v` → PASS

  **Commit**: YES
  - Message: `fix(export): add None guards to export functions`
  - Files: `export_engine.py`, `tests/test_none_guards.py`
  - Pre-commit: `pytest tests/test_none_guards.py -v`

---

- [ ] 9. Add None guards to gene panel functions

  **What to do**:
  - Write tests for edge cases in gene panel functions
  - Add validation to `score_panel()` (already has some, enhance):
    - Validate expression_df is not None before column access
    - Better error message when panel not found
  - Add validation to `plot_panel()`:
    - Check expression_df is not None
    - Check sample_conditions is not empty
    - Handle case where NO genes from panel are in expression data
  - Update platform code that calls gene panels (lines 501-516):
    - Catch ValueError exceptions properly
    - Track which panels failed and why
    - Show user-friendly message listing unavailable panels

  **Must NOT do**:
  - Do NOT change gene panel config format
  - Do NOT add new panels

  **Parallelizable**: YES (with Tasks 7, 8)

  **References**:

  **Pattern References**:
  - `gene_panels.py:100-139` - Existing validation in `score_panel()`
  - `gene_panels.py:162-175` - Existing validation in `plot_panel()`

  **Silent Failure Location**:
  - `rnaseq_analysis_platform.py:508-514` - Silent `pass` on ValueError

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_score_panel_none_expression_raises`
  - [ ] Test: `test_plot_panel_none_expression_raises`
  - [ ] Test: `test_plot_panel_empty_conditions_raises`
  - [ ] Test: `test_platform_tracks_failed_panels`
  - [ ] `pytest tests/test_none_guards.py::test_*panel* -v` → PASS

  **Commit**: YES
  - Message: `fix(panels): add None guards and improve error tracking`
  - Files: `gene_panels.py`, `rnaseq_analysis_platform.py`, `tests/test_none_guards.py`
  - Pre-commit: `pytest tests/test_none_guards.py -v`

---

### P2: UX ENHANCEMENT

- [ ] 10. Add confirmation messages for auto-detected columns/conditions

  **What to do**:
  - Add UI feedback when pattern-based detection succeeds:
    - "Detected fold change column: Bt10U/none.fc"
    - "Detected adjusted p-value column: Bt10U/none.bh.pval"
    - "Detected comparison: Bt10U vs none"
    - "Detected conditions from sample columns: none (3 samples), Bt10U (3 samples)"
  - Add expander showing detailed column mapping:
    ```
    Column Mapping:
    - Bt10U/none.fc → log2FoldChange
    - Bt10U/none.bh.pval → padj
    - Bt10U/none.raw.pval → pvalue
    ```
  - Add override option: "Use different column mapping" button that shows manual dropdowns

  **Must NOT do**:
  - Do NOT make messages dismissible (they're important context)
  - Do NOT add complex modal dialogs

  **Parallelizable**: NO (depends on 6)

  **References**:

  **Pattern References**:
  - `rnaseq_analysis_platform.py:195-198` - Success/warning message pattern
  - `rnaseq_analysis_platform.py:215-223` - Expander pattern for data preview

  **Acceptance Criteria**:

  **Manual Verification (Playwright browser)**:
  - [ ] Upload multi-format file
  - [ ] Verify: Confirmation messages appear in sidebar
  - [ ] Verify: Expander shows column mapping details
  - [ ] Screenshot: `.sisyphus/evidence/task-10-confirmation-ui.png`

  **Commit**: YES
  - Message: `feat(ui): add confirmation messages for auto-detection`
  - Files: `rnaseq_analysis_platform.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 11. Improve error messages with actionable suggestions

  **What to do**:
  - Enhance ParserValidationError messages with suggestions:
    - "Column 'log2FoldChange' not found. Did you mean 'Bt10U/none.fc'? The parser detected this as a fold change column."
    - "No gene column detected. Your file should have a column named 'Gene', 'Gene_Symbol', or similar."
  - Add context to visualization errors:
    - "Cannot create heatmap: expression data not available. This file contains only pre-analyzed DE results."
  - Add suggestions for common issues:
    - Empty DataFrame: "No genes passed the significance threshold. Try adjusting padj < {threshold} or |log2FC| > {threshold}."

  **Must NOT do**:
  - Do NOT create overly long error messages
  - Do NOT suggest impossible actions

  **Parallelizable**: YES (with Task 10)

  **References**:

  **Pattern References**:
  - `rnaseq_parser.py:63-69` - ParserValidationError with details
  - `rnaseq_parser.py:481-487` - Error message with available columns

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_error_message_suggests_similar_column`
  - [ ] Test: `test_error_message_includes_context`
  - [ ] `pytest tests/test_none_guards.py::test_error_message* -v` → PASS

  **Commit**: YES
  - Message: `feat(ux): improve error messages with actionable suggestions`
  - Files: `rnaseq_parser.py`, `visualizations.py`, `tests/test_none_guards.py`
  - Pre-commit: `pytest tests/ -v`

---

### P3: CLEANUP & FINAL VERIFICATION

- [ ] 12. Fix silent failure points

  **What to do**:
  - Fix line 509 (gene panels silent failure):
    - Track failed panels in a list
    - After loop, show info message: "Gene panels unavailable: {list}. Reason: insufficient genes in your data."
  - Fix line 288 (conftest ImportError):
    - Add logging or warning instead of silent pass
  - Audit for any other silent failures introduced

  **Must NOT do**:
  - Do NOT make failures block the entire analysis
  - Do NOT create noisy warnings for expected conditions

  **Parallelizable**: NO (cleanup task)

  **References**:

  **Silent Failure Locations**:
  - `rnaseq_analysis_platform.py:508-514`:
    ```python
    except ValueError:
        pass  # Skip panels with insufficient genes
    ```
  - `conftest.py:288` - ImportError silent pass

  **Acceptance Criteria**:

  **Manual Verification**:
  - [ ] Upload file with genes NOT in any panel
  - [ ] Verify: Info message appears listing unavailable panels
  - [ ] Verify: No silent failures in logs

  **Commit**: YES
  - Message: `fix(platform): replace silent failures with user feedback`
  - Files: `rnaseq_analysis_platform.py`, `conftest.py`
  - Pre-commit: `pytest tests/ -v`

---

- [ ] 13. Final integration test with all reference files

  **What to do**:
  - Create comprehensive integration test that:
    1. Loads each of the 4 reference files
    2. Verifies parsing succeeds for all
    3. Verifies visualizations render for all
    4. Verifies export works for all
  - Run full test suite and verify 100% pass
  - Manual smoke test with Streamlit app

  **Must NOT do**:
  - Do NOT add tests for features not implemented
  - Do NOT create flaky tests dependent on external services

  **Parallelizable**: NO (final verification, depends on all)

  **References**:

  **Reference Files**:
  - `Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx`
  - `Reference sequencing data/data3_Bu10_vs_none_fc2_&_raw.p.xlsx`
  - `Reference sequencing data/data3_my45_vs_none_fc2_&_raw.p.xlsx`
  - `Reference sequencing data/data3_SwX15_vs_none_fc2_&_raw.p.xlsx`

  **Test References**:
  - `tests/test_integration.py` - Existing integration test patterns

  **Acceptance Criteria**:

  **TDD**:
  - [ ] Test: `test_all_reference_files_parse_successfully`
  - [ ] Test: `test_all_reference_files_visualize_successfully`
  - [ ] `pytest tests/ -v` → ALL PASS (0 failures)

  **Manual Verification (Playwright browser)**:
  - [ ] Upload each of the 4 reference files
  - [ ] Verify: All parse without error
  - [ ] Verify: All show volcano plot, heatmap, PCA
  - [ ] Verify: Export to Excel works
  - [ ] Screenshots: `.sisyphus/evidence/task-13-{filename}.png`

  **Commit**: YES
  - Message: `test: add comprehensive integration tests for multi-format files`
  - Files: `tests/test_integration.py`
  - Pre-commit: `pytest tests/ -v`

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `fix(platform): guard None expression_df in visualization path` | rnaseq_analysis_platform.py, tests/test_none_guards.py | pytest tests/test_none_guards.py |
| 2 | `feat(parser): extend ParseResult to support multi-format files` | rnaseq_parser.py, tests/test_multiformat_parser.py | pytest tests/ |
| 3 | `feat(parser): add pattern-based DE column detection` | rnaseq_parser.py, tests/test_multiformat_parser.py | pytest tests/test_multiformat_parser.py |
| 4 | `feat(parser): add sample/condition auto-detection from column names` | rnaseq_parser.py, tests/test_multiformat_parser.py | pytest tests/test_multiformat_parser.py |
| 5 | `feat(parser): implement multi-format Excel parsing` | rnaseq_parser.py, tests/test_multiformat_parser.py | pytest tests/ |
| 6 | `feat(platform): integrate multi-format parser with UI` | rnaseq_analysis_platform.py, tests/test_integration.py | pytest tests/ |
| 7 | `fix(viz): add input validation to visualization functions` | visualizations.py, tests/test_none_guards.py | pytest tests/test_none_guards.py |
| 8 | `fix(export): add None guards to export functions` | export_engine.py, tests/test_none_guards.py | pytest tests/test_none_guards.py |
| 9 | `fix(panels): add None guards and improve error tracking` | gene_panels.py, rnaseq_analysis_platform.py, tests/test_none_guards.py | pytest tests/test_none_guards.py |
| 10 | `feat(ui): add confirmation messages for auto-detection` | rnaseq_analysis_platform.py | pytest tests/ |
| 11 | `feat(ux): improve error messages with actionable suggestions` | rnaseq_parser.py, visualizations.py, tests/test_none_guards.py | pytest tests/ |
| 12 | `fix(platform): replace silent failures with user feedback` | rnaseq_analysis_platform.py, conftest.py | pytest tests/ |
| 13 | `test: add comprehensive integration tests for multi-format files` | tests/test_integration.py | pytest tests/ |

---

## Success Criteria

### Verification Commands
```bash
# Run all tests
pytest tests/ -v

# Run specific test categories
pytest tests/test_multiformat_parser.py -v  # Parser tests
pytest tests/test_none_guards.py -v         # None guard tests
pytest tests/test_integration.py -v         # Integration tests

# Verify reference file parsing
python -c "
from rnaseq_parser import RNASeqParser
parser = RNASeqParser()
result = parser.parse('Reference sequencing data/data3_Bt10U_vs_none_fc2_&_raw.p.xlsx')
assert result.de_results_df is not None, 'DE results missing'
assert result.expression_df is not None, 'Expression data missing'
print('SUCCESS: Multi-format parsing works!')
"

# Launch app for manual testing
streamlit run rnaseq_analysis_platform.py
```

### Final Checklist
- [ ] All "Must Have" present:
  - [ ] P0 crash fix implemented
  - [ ] Multi-format parsing works
  - [ ] Pattern-based column detection works
  - [ ] Sample auto-detection works
  - [ ] None guards comprehensive
- [ ] All "Must NOT Have" absent:
  - [ ] No changes to de_analysis.py
  - [ ] No breaking changes to existing parsing
  - [ ] No new dependencies
- [ ] All tests pass: `pytest tests/ -v` → 0 failures
- [ ] All 4 reference files work end-to-end
