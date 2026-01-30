# RNA-seq Platform PhD-Level Feature Upgrade

## TL;DR

> **Quick Summary**: Add 8 advanced analysis features (GSEA, interpretation engine, cell deconvolution, expanded enrichment DBs, outlier detection, batch effect assessment, methods text, enrichment dot plots) to the existing Streamlit RNA-seq platform by creating 3 new modules and extending 3 existing ones.
> 
> **Deliverables**:
> - New module: `interpretation_engine.py` (automated insights + methods text generation)
> - New module: `cell_deconvolution.py` (NNLS-based cell type estimation)
> - New module: `advanced_qc.py` (outlier detection + batch effect assessment)
> - Extended: `pathway_enrichment.py` (GSEA + 5 additional databases)
> - Extended: `visualizations.py` (enrichment dot plot + deconvolution stacked bar)
> - Extended: `rnaseq_analysis_platform.py` (new tabs + pipeline integration)
> - Tests for all new modules
> 
> **Estimated Effort**: Large
> **Parallel Execution**: YES - 3 waves
> **Critical Path**: Task 1 (enrichment expansion) → Task 5 (GSEA) → Task 8 (UI integration)

---

## Context

### Original Request
User wants to upgrade an existing, working RNA-seq Streamlit platform to PhD/expert level by integrating 8 new features from a reference implementation (`rnaseq_platform_implementation.py`).

### Research Findings
- **Main app**: 4-step wizard (Upload→Metadata→Compare→Results) with 9 tabs in Results (lines 532-534)
- **DEResult dataclass** (de_analysis.py:69-82): `results_df`, `normalized_counts`, `log_normalized_counts`, `dds`, `comparison`, `n_significant`, `warnings`
- **EnrichmentResult** (export_engine.py:25-32): `go_results`, `kegg_results`, `genes_used`, `selection_note`, `error` — used by both main app and export engine
- **PathwayEnrichment** (pathway_enrichment.py:20): Only GO_BP_2023 + KEGG_2021 via Enrichr
- **All visualizations**: Plotly exclusively (`go.Figure` return pattern)
- **Reference enrichment dotplot uses matplotlib** — must be converted to Plotly
- **Reference cell proportion plot uses matplotlib** — must be converted to Plotly
- **gseapy, scipy, sklearn**: Already installed as dependencies
- **Analysis pipeline** (lines 411-498): Runs DE → Enrichment → Visualizations sequentially in Step 3

### Key Architecture Constraints
- EnrichmentResult is defined in `export_engine.py` and imported by main app — adding new enrichment types (GSEA, additional DBs) must extend this dataclass
- The analysis pipeline in Step 3 (lines 411-498) needs extension for new analyses
- Session state keys must be added for new results
- New Results tabs must be added to the tab list at line 532

---

## Work Objectives

### Core Objective
Add 8 PhD-level analysis features to the existing platform without breaking current functionality.

### Concrete Deliverables
- 3 new Python modules
- 3 modified Python modules  
- 1 new test file
- New Results tabs: "🔍 GSEA", "🧠 Interpretation", "🧫 Deconvolution", "⚠️ QC Advanced"

### Definition of Done
- [ ] `pytest tests/ -v` — all existing + new tests pass
- [ ] Existing workflow (Upload → Metadata → Compare → Results) unchanged
- [ ] All 9 original tabs still work identically
- [ ] New tabs appear and render without errors
- [ ] All new visualizations use Plotly (zero matplotlib in UI)
- [ ] GSEA runs successfully with gseapy.prerank
- [ ] Enrichment tab shows all available databases (GO_BP, GO_MF, GO_CC, KEGG, Reactome, WikiPathway, Hallmark)

### Must Have
- GSEA with ranked gene lists
- Enrichment dot plot (Plotly)
- Interpretation engine with automated insights
- Cell type deconvolution with built-in skin signature matrix
- Outlier detection in PCA space
- Batch effect variance decomposition
- Methods text generator
- 5 additional enrichment databases

### Must NOT Have (Guardrails)
- NO matplotlib/seaborn in any UI-facing code — Plotly only
- NO changes to the DEResult dataclass
- NO changes to the existing 4-step wizard flow
- NO changes to existing visualization functions (volcano, heatmap, PCA, MA, correlation)
- NO breaking changes to ExportEngine or existing export functionality
- NO new pip dependencies beyond what's already installed
- NO removal of existing Enrichr-based ORA — GSEA is additive
- NO over-engineered abstractions — keep module simplicity consistent with existing codebase
- NO interactive/blocking operations during analysis (all background computation)

---

## Verification Strategy

### Test Decision
- **Infrastructure exists**: YES (pytest, tests/ directory with 8 test files)
- **User wants tests**: YES (tests after implementation)
- **Framework**: pytest

### Automated Verification

Each TODO includes verification via:
- `pytest tests/ -v` for unit tests
- `python -c "from module import Class; ..."` for import verification
- `streamlit run rnaseq_analysis_platform.py` + manual spot-check for UI integration

---

## Execution Strategy

### Parallel Execution Waves

```
Wave 1 (Start Immediately — independent backend modules):
├── Task 1: Extend pathway_enrichment.py (add databases + GSEA)
├── Task 2: Create interpretation_engine.py
├── Task 3: Create cell_deconvolution.py
└── Task 4: Create advanced_qc.py

Wave 2 (After Wave 1 — depends on new modules):
├── Task 5: Add enrichment dot plot + deconvolution bar chart to visualizations.py
└── Task 6: Extend EnrichmentResult dataclass in export_engine.py

Wave 3 (After Wave 2 — full integration):
├── Task 7: Integrate all into rnaseq_analysis_platform.py
└── Task 8: Write tests for all new modules
```

### Dependency Matrix

| Task | Depends On | Blocks | Can Parallelize With |
|------|------------|--------|---------------------|
| 1 | None | 5, 6, 7 | 2, 3, 4 |
| 2 | None | 7 | 1, 3, 4 |
| 3 | None | 5, 7 | 1, 2, 4 |
| 4 | None | 7 | 1, 2, 3 |
| 5 | 1, 3 | 7 | 6 |
| 6 | 1 | 7 | 5 |
| 7 | 1-6 | 8 | None |
| 8 | 1-7 | None | None |

---

## TODOs

- [ ] 1. Extend `pathway_enrichment.py` — Add GSEA + Additional Databases

  **What to do**:
  - Add `run_gsea(self, gene_ranks: pd.Series, database: str, min_size: int = 15, max_size: int = 500, permutations: int = 1000) -> Tuple[pd.DataFrame, Optional[str]]` method to `PathwayEnrichment` class. Follow the existing error-handling pattern (return tuple with optional error string). Use `gseapy.prerank`.
  - Add `create_ranking_metric(de_results: pd.DataFrame) -> pd.Series` as a `@staticmethod`. Ranking = `-log10(pvalue) * sign(log2FoldChange)`. Handle NaN/inf by dropping.
  - Add convenience methods following the existing `get_go_enrichment`/`get_kegg_enrichment` pattern for each new database:
    - `get_go_mf_enrichment(genes)` → `GO_Molecular_Function_2023`
    - `get_go_cc_enrichment(genes)` → `GO_Cellular_Component_2023`
    - `get_reactome_enrichment(genes)` → `Reactome_2022`
    - `get_wikipathway_enrichment(genes)` → `WikiPathway_2023_Human`
    - `get_hallmark_enrichment(genes)` → `MSigDB_Hallmark_2020`
  - Add a `AVAILABLE_DATABASES` class constant dict mapping display names to Enrichr library names.
  - GSEA results DataFrame should have columns: `Term`, `NES`, `P-value`, `FDR q-value`, `Lead_genes`.

  **Must NOT do**:
  - Do NOT modify existing `get_go_enrichment` or `get_kegg_enrichment` methods
  - Do NOT change the `format_results` method behavior
  - Do NOT change the `select_genes_for_enrichment` logic

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 2, 3, 4)
  - **Blocks**: Tasks 5, 6, 7
  - **Blocked By**: None

  **References**:
  - `pathway_enrichment.py:20-181` — Existing class structure, error handling pattern (return `Tuple[pd.DataFrame, Optional[str]]`)
  - `pathway_enrichment.py:83-121` — `run_enrichment` method pattern to follow for new database methods
  - `pathway_enrichment.py:123-154` — `get_go_enrichment` and `get_kegg_enrichment` as templates for new convenience methods
  - `rnaseq_platform_implementation.py:585-638` — Reference GSEA implementation (adapt `run_gsea` and `create_ranking_metric`)
  - `de_analysis.py:69-82` — DEResult dataclass showing `results_df` columns (`gene`, `baseMean`, `log2FoldChange`, `pvalue`, `padj`)

  **Acceptance Criteria**:
  ```bash
  # Verify import works
  python -c "from pathway_enrichment import PathwayEnrichment; pe = PathwayEnrichment(); print('AVAILABLE_DATABASES:', list(pe.AVAILABLE_DATABASES.keys())); print('OK')"
  # Assert: Output includes GO_BP, GO_MF, GO_CC, KEGG, Reactome, WikiPathway, Hallmark
  
  # Verify GSEA method exists and is callable
  python -c "
  from pathway_enrichment import PathwayEnrichment
  import pandas as pd, numpy as np
  pe = PathwayEnrichment()
  # Test create_ranking_metric with mock data
  mock_de = pd.DataFrame({'pvalue': [0.001, 0.05, 0.5], 'log2FoldChange': [2.0, -1.5, 0.1]})
  ranks = pe.create_ranking_metric(mock_de)
  assert len(ranks) == 3, 'Should have 3 ranked genes'
  assert ranks.iloc[0] > 0, 'First rank should be positive (upregulated)'
  print('Ranking metric OK')
  "
  
  # Verify new database methods exist
  python -c "
  from pathway_enrichment import PathwayEnrichment
  pe = PathwayEnrichment()
  assert hasattr(pe, 'get_go_mf_enrichment')
  assert hasattr(pe, 'get_go_cc_enrichment')
  assert hasattr(pe, 'get_reactome_enrichment')
  assert hasattr(pe, 'get_wikipathway_enrichment')
  assert hasattr(pe, 'get_hallmark_enrichment')
  assert hasattr(pe, 'run_gsea')
  print('All new methods exist')
  "
  ```

  **Commit**: YES
  - Message: `feat(enrichment): add GSEA and 5 additional enrichment databases`
  - Files: `pathway_enrichment.py`
  - Pre-commit: `python -c "from pathway_enrichment import PathwayEnrichment; print('OK')"`

---

- [ ] 2. Create `interpretation_engine.py` — Automated Insights + Methods Text

  **What to do**:
  - Create new file `interpretation_engine.py` with:
  - `Interpretation` dataclass: `category` (str), `level` (str: 'critical'/'important'/'informative'), `title` (str), `description` (str), `evidence` (Dict), `recommendations` (List[str])
  - `InterpretationEngine` class with:
    - `interpret_de_results(de_results_df: pd.DataFrame, contrast_name: str, padj_threshold: float = 0.05, lfc_threshold: float = 1.0) -> List[Interpretation]` — Generates insights about DE summary (total sig genes, up/down counts), effect size assessment (modest if median |LFC| < 0.5), and balance assessment
    - `interpret_enrichment(enrichment_df: pd.DataFrame) -> List[Interpretation]` — Detects immune-related, metabolic, signaling pathway enrichment via keyword matching
    - `generate_methods_text(params: Dict) -> str` — Returns publication-ready methods section as markdown string
  - Accept `pd.DataFrame` inputs (NOT the DEResult dataclass) so the engine is decoupled from DE module
  - Follow existing codebase pattern: no class-level state beyond `self.interpretations` list

  **Must NOT do**:
  - Do NOT import DEResult or depend on de_analysis.py
  - Do NOT use any visualization code
  - Do NOT add overly specific dermatology keywords — keep generic for any biology

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 3, 4)
  - **Blocks**: Task 7
  - **Blocked By**: None

  **References**:
  - `rnaseq_platform_implementation.py:111-119` — `Interpretation` dataclass definition
  - `rnaseq_platform_implementation.py:1055-1189` — Full reference implementation of InterpretationEngine (adapt logic, keep interface)
  - `de_analysis.py:69-82` — DEResult dataclass fields (note: `results_df` has columns `gene`, `baseMean`, `log2FoldChange`, `lfcSE`, `stat`, `pvalue`, `padj`)
  - `pathway_enrichment.py:157-181` — Enrichment result DataFrame format: columns `Term`, `Overlap`, `P-value`, `Adjusted P-value`, `Genes`

  **Acceptance Criteria**:
  ```bash
  python -c "
  from interpretation_engine import InterpretationEngine, Interpretation
  import pandas as pd
  
  engine = InterpretationEngine()
  
  # Test DE interpretation
  mock_de = pd.DataFrame({
      'gene': ['A', 'B', 'C', 'D'],
      'log2FoldChange': [3.0, -2.5, 0.1, 1.5],
      'padj': [0.001, 0.01, 0.8, 0.03],
      'pvalue': [0.0001, 0.001, 0.5, 0.01],
      'baseMean': [100, 200, 50, 150]
  })
  results = engine.interpret_de_results(mock_de, 'Treatment_vs_Control')
  assert len(results) >= 1, 'Should produce at least 1 interpretation'
  assert results[0].category == 'de'
  assert 'upregulated' in results[0].description.lower() or 'differentially expressed' in results[0].description.lower()
  print('DE interpretation OK')
  
  # Test methods text
  text = engine.generate_methods_text({'de_tool': 'PyDESeq2', 'padj': 0.05, 'lfc': 1.0})
  assert 'PyDESeq2' in text
  print('Methods text OK')
  "
  ```

  **Commit**: YES
  - Message: `feat: add interpretation engine with automated insights and methods text generation`
  - Files: `interpretation_engine.py`
  - Pre-commit: `python -c "from interpretation_engine import InterpretationEngine; print('OK')"`

---

- [ ] 3. Create `cell_deconvolution.py` — NNLS Cell Type Estimation

  **What to do**:
  - Create new file `cell_deconvolution.py` with:
  - `DeconvolutionResult` dataclass: `proportions` (pd.DataFrame, samples × cell types), `method` (str), `goodness_of_fit` (Optional[pd.Series])
  - `CellTypeDeconvolution` class with:
    - `run_nnls(expression: pd.DataFrame, signature: pd.DataFrame) -> DeconvolutionResult` — NNLS deconvolution using `scipy.optimize.nnls`. Finds common genes between expression and signature, requires ≥50 overlap. Normalizes coefficients to sum to 1 per sample.
    - `get_default_skin_signature() -> pd.DataFrame` — Returns a built-in signature matrix for skin tissue with cell types: Keratinocyte, Fibroblast, Melanocyte, Immune (T-cell), Endothelial. Use ~20-30 well-known marker genes per type with realistic relative expression values.
  - Expression input: genes × samples (gene names as index, sample names as columns)
  - Signature input: genes × cell types (gene names as index, cell type names as columns)

  **Must NOT do**:
  - Do NOT require external signature matrix files — built-in default is mandatory
  - Do NOT import matplotlib/seaborn
  - Do NOT make this Streamlit-dependent

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2, 4)
  - **Blocks**: Tasks 5, 7
  - **Blocked By**: None

  **References**:
  - `rnaseq_platform_implementation.py:104-108` — `DeconvolutionResult` dataclass
  - `rnaseq_platform_implementation.py:779-849` — Reference NNLS implementation
  - `de_analysis.py:69-82` — DEResult has `normalized_counts` (samples × genes) and `log_normalized_counts` — note: deconvolution needs genes × samples, so transpose needed at call site
  - `rnaseq_analysis_platform.py:472-478` — How normalized counts are accessed: `first_res.log_normalized_counts` is samples × genes

  **Acceptance Criteria**:
  ```bash
  python -c "
  from cell_deconvolution import CellTypeDeconvolution, DeconvolutionResult
  import pandas as pd, numpy as np
  
  deconv = CellTypeDeconvolution()
  
  # Test default signature
  sig = deconv.get_default_skin_signature()
  assert sig.shape[0] >= 20, f'Signature should have >=20 genes, got {sig.shape[0]}'
  assert sig.shape[1] >= 4, f'Signature should have >=4 cell types, got {sig.shape[1]}'
  print(f'Signature: {sig.shape[0]} genes x {sig.shape[1]} cell types')
  
  # Test NNLS with mock data using signature genes
  genes = sig.index.tolist()
  mock_expr = pd.DataFrame(
      np.random.rand(len(genes), 3) * 100,
      index=genes,
      columns=['S1', 'S2', 'S3']
  )
  result = deconv.run_nnls(mock_expr, sig)
  assert isinstance(result, DeconvolutionResult)
  assert result.proportions.shape == (3, sig.shape[1])
  assert np.allclose(result.proportions.sum(axis=1), 1.0, atol=0.01), 'Proportions should sum to ~1'
  print('NNLS deconvolution OK')
  "
  ```

  **Commit**: YES
  - Message: `feat: add NNLS cell type deconvolution with skin tissue signature`
  - Files: `cell_deconvolution.py`
  - Pre-commit: `python -c "from cell_deconvolution import CellTypeDeconvolution; print('OK')"`

---

- [ ] 4. Create `advanced_qc.py` — Outlier Detection + Batch Effect Assessment

  **What to do**:
  - Create new file `advanced_qc.py` with:
  - `AdvancedQC` class with:
    - `detect_outliers(pca_coords: pd.DataFrame, threshold_sd: float = 3.0) -> Tuple[List[str], pd.Series]` — Input is DataFrame with PC1, PC2 columns and sample names as index. Calculates Euclidean distance from centroid, flags samples > threshold_sd standard deviations away. Returns (list of outlier sample names, Series of all distances).
    - `assess_batch_effects(pca_coords: pd.DataFrame, batch_labels: pd.Series) -> Dict[str, float]` — Calculates variance explained by batch variable for each PC using sum-of-squares decomposition. Returns dict mapping PC name to fraction of variance explained by batch (0-1).
    - `create_outlier_plot(pca_coords: pd.DataFrame, sample_conditions: Dict[str, str], outliers: List[str], distances: pd.Series) -> go.Figure` — Plotly PCA scatter plot highlighting outlier samples with red markers and distance annotations.
    - `create_batch_effect_plot(batch_results: Dict[str, float]) -> go.Figure` — Plotly bar chart showing variance explained by batch per PC.
  - All plots must return `plotly.graph_objects.Figure`

  **Must NOT do**:
  - Do NOT import matplotlib/seaborn
  - Do NOT duplicate existing PCA computation — take pre-computed PCA coordinates as input
  - Do NOT modify existing `qc_plots.py`

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 1 (with Tasks 1, 2, 3)
  - **Blocks**: Task 7
  - **Blocked By**: None

  **References**:
  - `rnaseq_platform_implementation.py:736-772` — Reference outlier detection and batch effect assessment logic
  - `visualizations.py:251-347` — Existing `create_pca_plot` uses `sklearn.decomposition.PCA`, computes PC1/PC2 and returns `go.Figure` — note how PCA coordinates are computed (lines 270-290) so the advanced QC can reuse same coordinates
  - `qc_plots.py:9-153` — Existing QC plots pattern (all return `go.Figure`)
  - `visualizations.py:18-132` — Plotly scatter plot patterns (volcano plot uses `go.Scatter` with custom colors, markers)

  **Acceptance Criteria**:
  ```bash
  python -c "
  from advanced_qc import AdvancedQC
  import pandas as pd, numpy as np
  import plotly.graph_objects as go
  
  qc = AdvancedQC()
  
  # Mock PCA coordinates (one outlier at [10, 10])
  pca_df = pd.DataFrame({
      'PC1': [1.0, 1.2, 0.8, 1.1, 10.0],
      'PC2': [0.5, 0.6, 0.4, 0.5, 10.0]
  }, index=['S1', 'S2', 'S3', 'S4', 'Outlier'])
  
  outliers, distances = qc.detect_outliers(pca_df, threshold_sd=2.0)
  assert 'Outlier' in outliers, f'Should detect Outlier, got: {outliers}'
  assert len(outliers) == 1, f'Should find exactly 1 outlier, got {len(outliers)}'
  print(f'Outlier detection OK: {outliers}')
  
  # Batch effect assessment
  batch = pd.Series(['A', 'A', 'B', 'B', 'A'], index=pca_df.index)
  batch_results = qc.assess_batch_effects(pca_df, batch)
  assert 'PC1' in batch_results
  assert 0 <= batch_results['PC1'] <= 1
  print(f'Batch effect OK: {batch_results}')
  
  # Verify plots return Plotly figures
  conds = {s: 'ctrl' for s in pca_df.index}
  fig1 = qc.create_outlier_plot(pca_df, conds, outliers, distances)
  assert isinstance(fig1, go.Figure), 'Outlier plot must be Plotly Figure'
  
  fig2 = qc.create_batch_effect_plot(batch_results)
  assert isinstance(fig2, go.Figure), 'Batch plot must be Plotly Figure'
  print('All plots are Plotly Figures')
  "
  ```

  **Commit**: YES
  - Message: `feat: add advanced QC with outlier detection and batch effect assessment`
  - Files: `advanced_qc.py`
  - Pre-commit: `python -c "from advanced_qc import AdvancedQC; print('OK')"`

---

- [ ] 5. Add Plotly enrichment dot plot + deconvolution bar chart to `visualizations.py`

  **What to do**:
  - Add `create_enrichment_dotplot(enrichment_df: pd.DataFrame, top_n: int = 20, title: str = "Enrichment Results") -> go.Figure` to `visualizations.py`:
    - X-axis: `-log10(Adjusted P-value)` (using column name from existing format: `Adjusted P-value`)
    - Y-axis: Term names (sorted by significance)
    - Dot size: Proportional to gene count (parse from `Overlap` column: "3/150" → 3)
    - Dot color: `-log10(Adjusted P-value)` with viridis-like colorscale
    - Show top_n terms, most significant at top
    - Handle empty DataFrame gracefully (return empty figure with message)
  - Add `create_deconvolution_plot(proportions_df: pd.DataFrame, sample_conditions: Optional[Dict[str, str]] = None) -> go.Figure`:
    - Stacked bar chart with samples on x-axis, proportions on y-axis
    - Each cell type as a different colored bar segment
    - Use Plotly's built-in color sequences
    - If sample_conditions provided, group/sort samples by condition

  **Must NOT do**:
  - Do NOT import matplotlib or seaborn
  - Do NOT modify any existing visualization functions
  - Do NOT change existing function signatures

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Task 6)
  - **Blocks**: Task 7
  - **Blocked By**: Tasks 1, 3

  **References**:
  - `visualizations.py:18-132` — Existing Plotly patterns (`go.Scatter`, `fig.update_layout`)
  - `visualizations.py:134-249` — Heatmap pattern (complex Plotly figure construction)
  - `rnaseq_platform_implementation.py:996-1048` — Reference matplotlib dotplot and stacked bar (adapt logic to Plotly)
  - `pathway_enrichment.py:157-181` — Enrichment result DataFrame columns: `Term`, `Overlap`, `P-value`, `Adjusted P-value`, `Genes` — the dot plot must work with these exact column names
  - `cell_deconvolution.py` (Task 3 output) — `DeconvolutionResult.proportions` is samples × cell types DataFrame

  **Acceptance Criteria**:
  ```bash
  python -c "
  from visualizations import create_enrichment_dotplot, create_deconvolution_plot
  import pandas as pd, numpy as np
  import plotly.graph_objects as go
  
  # Test enrichment dotplot
  mock_enr = pd.DataFrame({
      'Term': ['Pathway A', 'Pathway B', 'Pathway C'],
      'Overlap': ['5/100', '3/80', '10/200'],
      'P-value': [0.001, 0.01, 0.05],
      'Adjusted P-value': [0.01, 0.05, 0.1],
      'Genes': ['A;B;C;D;E', 'F;G;H', 'I;J;K;L;M;N;O;P;Q;R']
  })
  fig = create_enrichment_dotplot(mock_enr)
  assert isinstance(fig, go.Figure), 'Dotplot must be Plotly Figure'
  print('Enrichment dotplot OK')
  
  # Test empty enrichment
  empty_fig = create_enrichment_dotplot(pd.DataFrame())
  assert isinstance(empty_fig, go.Figure), 'Empty dotplot must still be Plotly Figure'
  print('Empty enrichment dotplot OK')
  
  # Test deconvolution plot
  mock_props = pd.DataFrame({
      'Keratinocyte': [0.5, 0.6, 0.4],
      'Fibroblast': [0.3, 0.2, 0.35],
      'Immune': [0.2, 0.2, 0.25]
  }, index=['S1', 'S2', 'S3'])
  fig2 = create_deconvolution_plot(mock_props)
  assert isinstance(fig2, go.Figure), 'Deconv plot must be Plotly Figure'
  print('Deconvolution plot OK')
  "
  ```

  **Commit**: YES
  - Message: `feat(viz): add Plotly enrichment dotplot and deconvolution stacked bar chart`
  - Files: `visualizations.py`
  - Pre-commit: `python -c "from visualizations import create_enrichment_dotplot, create_deconvolution_plot; print('OK')"`

---

- [ ] 6. Extend `EnrichmentResult` dataclass in `export_engine.py`

  **What to do**:
  - Add optional fields to the `EnrichmentResult` dataclass (lines 25-32 of `export_engine.py`) for the new enrichment databases and GSEA:
    ```python
    go_mf_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    go_cc_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    reactome_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    wikipathway_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    hallmark_results: pd.DataFrame = field(default_factory=pd.DataFrame)
    gsea_results: Dict[str, pd.DataFrame] = field(default_factory=dict)
    ```
  - All new fields use `field(default_factory=...)` so existing code creating `EnrichmentResult(go_results=..., kegg_results=..., genes_used=..., selection_note=...)` continues to work without modification
  - Add `from dataclasses import field` import if not already present

  **Must NOT do**:
  - Do NOT change existing field names or types
  - Do NOT change existing field order (would break positional construction at line 441 of main app)
  - Do NOT modify ExportEngine methods — new fields are optional and additive

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: YES
  - **Parallel Group**: Wave 2 (with Task 5)
  - **Blocks**: Task 7
  - **Blocked By**: Task 1

  **References**:
  - `export_engine.py:24-32` — Current `EnrichmentResult` dataclass definition
  - `rnaseq_analysis_platform.py:441-455` — Where `EnrichmentResult` is constructed (uses keyword args — safe for adding new optional fields)
  - `export_engine.py:1-20` — Import section (check if `field` is imported from dataclasses)

  **Acceptance Criteria**:
  ```bash
  python -c "
  from export_engine import EnrichmentResult
  import pandas as pd
  
  # Existing construction still works (backward compatible)
  old_style = EnrichmentResult(
      go_results=pd.DataFrame(),
      kegg_results=pd.DataFrame(),
      genes_used=['A', 'B'],
      selection_note='test'
  )
  assert old_style.error is None, 'Default error should be None'
  print('Backward compatible: OK')
  
  # New fields accessible with defaults
  assert isinstance(old_style.gsea_results, dict), 'GSEA results should default to empty dict'
  assert old_style.reactome_results.empty, 'Reactome should default to empty DataFrame'
  assert old_style.hallmark_results.empty, 'Hallmark should default to empty DataFrame'
  print('New fields with defaults: OK')
  
  # New construction with all fields
  full = EnrichmentResult(
      go_results=pd.DataFrame({'Term': ['A']}),
      kegg_results=pd.DataFrame({'Term': ['B']}),
      genes_used=['X'],
      selection_note='full test',
      reactome_results=pd.DataFrame({'Term': ['C']}),
      gsea_results={'Hallmark': pd.DataFrame({'Term': ['D']})}
  )
  assert len(full.reactome_results) == 1
  assert 'Hallmark' in full.gsea_results
  print('Full construction: OK')
  "
  ```

  **Commit**: YES
  - Message: `feat(export): extend EnrichmentResult with GSEA and additional database fields`
  - Files: `export_engine.py`
  - Pre-commit: `pytest tests/ -x -q`

---

- [ ] 7. Integrate all features into `rnaseq_analysis_platform.py`

  **What to do**:

  **A. New imports** (add after existing imports, ~line 27):
  ```python
  from interpretation_engine import InterpretationEngine
  from cell_deconvolution import CellTypeDeconvolution
  from advanced_qc import AdvancedQC
  from visualizations import create_enrichment_dotplot, create_deconvolution_plot
  ```

  **B. New session_state keys** (add to initialization block, ~line 46-67):
  ```python
  "gsea_results": {},
  "interpretations": [],
  "deconvolution_results": None,
  "outlier_info": None,
  "batch_effects": None,
  "pca_coordinates": None,
  ```

  **C. Extend analysis pipeline in Step 3** (after existing enrichment block at line 456, before visualizations at line 458):
  
  1. **Additional enrichment databases**: After the existing GO/KEGG enrichment loop, add a second pass that runs each additional database (`go_mf`, `go_cc`, `reactome`, `wikipathway`, `hallmark`) and stores results in the corresponding new `EnrichmentResult` fields.
  
  2. **GSEA**: After enrichment, for each comparison:
     - Create ranking metric from `de_res.results_df` using `PathwayEnrichment.create_ranking_metric`
     - Run GSEA against `MSigDB_Hallmark_2020` (default)
     - Store in `enrichment_results[comp].gsea_results`
     - Wrap in try/except to not block pipeline on GSEA failure
  
  3. **Interpretation**: After enrichment:
     - Create `InterpretationEngine()`
     - For each comparison, call `interpret_de_results` and `interpret_enrichment`
     - Store in `st.session_state["interpretations"]`
  
  4. **Cell deconvolution**: After visualization generation:
     - Create `CellTypeDeconvolution()`
     - Get default skin signature
     - Run NNLS on `norm_counts.T` (transpose from samples×genes to genes×samples)
     - Store result in `st.session_state["deconvolution_results"]`
     - Wrap in try/except
  
  5. **Outlier detection + batch effects**: During or after PCA computation:
     - Extract PCA coordinates (recompute from `create_pca_plot` or compute separately using sklearn PCA on `norm_counts`)
     - Store `pca_coordinates` in session_state
     - Run `AdvancedQC.detect_outliers` and store result
     - If metadata has a "batch" column, run `assess_batch_effects`

  **D. Expand Results tabs** (modify line 532):
  Change from 9 tabs to 13 tabs:
  ```python
  tab_summary, tab_de, tab_volcano, tab_ma, tab_heatmap, tab_pca, tab_enrichment, tab_gsea, tab_interpret, tab_deconv, tab_adv_qc, tab_panels, tab_export = st.tabs([
      "📊 Summary", "📋 DE Results", "🌋 Volcano", "📈 MA Plot",
      "🗺️ Heatmap", "🔬 PCA", "🧬 Enrichment", "🔍 GSEA",
      "🧠 Interpretation", "🧫 Deconvolution", "⚠️ QC Advanced",
      "🎯 Gene Panels", "📥 Export"
  ])
  ```

  **E. New tab contents**:

  **🧬 Enrichment tab** (enhance existing, lines 590-608):
  - Keep existing GO_BP and KEGG display
  - Add expandable sections for each additional database
  - Add enrichment dot plot visualization for GO and KEGG using `create_enrichment_dotplot`

  **🔍 GSEA tab** (new):
  - Show GSEA results DataFrame if available
  - Allow database selection via `st.selectbox` from `AVAILABLE_DATABASES`
  - Display NES (Normalized Enrichment Score) ranking
  - Show enrichment dotplot for GSEA results (adapt column names)

  **🧠 Interpretation tab** (new):
  - Display each `Interpretation` as an expander with icon based on level (🔴 critical, 🟡 important, 🔵 informative)
  - Show evidence dict as formatted details
  - Show recommendations as bullet list
  - Add "Generate Methods Text" button that generates and displays methods section
  - Add copy-to-clipboard via `st.code` block

  **🧫 Deconvolution tab** (new):
  - If `deconvolution_results` exists, show stacked bar chart via `create_deconvolution_plot`
  - Show proportions DataFrame below the chart
  - Add download button for proportions CSV
  - If deconvolution failed, show info message

  **⚠️ QC Advanced tab** (new):
  - If outliers detected, show warning with outlier sample names
  - Show outlier detection plot via `create_outlier_plot`
  - If batch effects computed, show batch effect bar chart via `create_batch_effect_plot`
  - Show distance table for all samples

  **Must NOT do**:
  - Do NOT modify existing tab content for Summary, DE Results, Volcano, MA, Heatmap, PCA, Gene Panels, Export
  - Do NOT change the 4-step wizard navigation
  - Do NOT change how existing modules are called
  - Do NOT let any new analysis step failure block the entire pipeline — always try/except with user-friendly error messages
  - Do NOT use matplotlib `st.pyplot` — all plots via `st.plotly_chart`

  **Recommended Agent Profile**:
  - **Category**: `unspecified-high`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (sequential)
  - **Blocks**: Task 8
  - **Blocked By**: Tasks 1-6

  **References**:
  - `rnaseq_analysis_platform.py:1-30` — Import section
  - `rnaseq_analysis_platform.py:46-67` — Session state initialization
  - `rnaseq_analysis_platform.py:411-498` — Analysis pipeline in Step 3 (insert new analyses here)
  - `rnaseq_analysis_platform.py:509-667` — Step 4 Results section (add new tabs here)
  - `rnaseq_analysis_platform.py:532-534` — Tab definition line (expand this)
  - `rnaseq_analysis_platform.py:590-608` — Existing enrichment tab (enhance with dotplots)
  - `rnaseq_analysis_platform.py:438-455` — How EnrichmentResult is constructed (extend for new DBs)
  - `pathway_enrichment.py` (Task 1 output) — New methods: `run_gsea`, `create_ranking_metric`, `get_go_mf_enrichment`, etc.
  - `interpretation_engine.py` (Task 2 output) — `InterpretationEngine` API
  - `cell_deconvolution.py` (Task 3 output) — `CellTypeDeconvolution` API
  - `advanced_qc.py` (Task 4 output) — `AdvancedQC` API
  - `visualizations.py` (Task 5 output) — `create_enrichment_dotplot`, `create_deconvolution_plot`

  **Acceptance Criteria**:
  ```bash
  # Verify imports work
  python -c "
  import rnaseq_analysis_platform
  print('Main app imports OK')
  "
  
  # Verify new session state keys
  python -c "
  # Check that the file contains the new session state keys
  with open('rnaseq_analysis_platform.py') as f:
      content = f.read()
  for key in ['gsea_results', 'interpretations', 'deconvolution_results', 'outlier_info']:
      assert key in content, f'Missing session_state key: {key}'
  print('Session state keys present')
  "
  
  # Verify tab count
  python -c "
  with open('rnaseq_analysis_platform.py') as f:
      content = f.read()
  assert 'GSEA' in content, 'GSEA tab missing'
  assert 'Interpretation' in content, 'Interpretation tab missing'
  assert 'Deconvolution' in content, 'Deconvolution tab missing'
  assert 'QC Advanced' in content, 'QC Advanced tab missing'
  print('All new tabs present')
  "
  
  # Verify no matplotlib usage
  python -c "
  with open('rnaseq_analysis_platform.py') as f:
      content = f.read()
  assert 'st.pyplot' not in content, 'Found st.pyplot — must use st.plotly_chart only'
  assert 'import matplotlib' not in content, 'Found matplotlib import in main app'
  print('No matplotlib in main app: OK')
  "
  ```

  **Commit**: YES
  - Message: `feat: integrate GSEA, interpretation, deconvolution, and advanced QC into Streamlit UI`
  - Files: `rnaseq_analysis_platform.py`
  - Pre-commit: `python -c "import rnaseq_analysis_platform; print('OK')"`

---

- [ ] 8. Write tests for all new modules

  **What to do**:
  - Create `tests/test_new_features.py` with pytest tests covering:
  
  **PathwayEnrichment extensions**:
  - `test_create_ranking_metric` — verify ranking metric computation with known values
  - `test_available_databases` — verify all 7 databases in `AVAILABLE_DATABASES`
  - `test_gsea_method_exists` — verify `run_gsea` is callable
  - `test_new_enrichment_methods_exist` — verify all 5 new getter methods exist
  
  **InterpretationEngine**:
  - `test_interpret_de_results` — verify returns list of Interpretation objects with correct fields
  - `test_interpret_enrichment_empty` — verify handles empty DataFrame
  - `test_interpret_enrichment_immune` — verify detects immune keywords
  - `test_generate_methods_text` — verify generates string containing tool name and parameters
  - `test_modest_effect_size_detection` — verify detects when median |LFC| < 0.5
  
  **CellTypeDeconvolution**:
  - `test_default_skin_signature` — verify shape and content
  - `test_nnls_basic` — verify returns DeconvolutionResult with correct dimensions
  - `test_nnls_proportions_sum_to_one` — verify normalization
  - `test_nnls_insufficient_overlap` — verify raises ValueError when < 50 genes overlap
  
  **AdvancedQC**:
  - `test_detect_outliers_basic` — verify detects obvious outlier
  - `test_detect_outliers_none` — verify returns empty list when no outliers
  - `test_batch_effects` — verify returns dict with PC keys and float values 0-1
  - `test_outlier_plot_returns_plotly` — verify returns go.Figure
  - `test_batch_plot_returns_plotly` — verify returns go.Figure
  
  **Visualization extensions**:
  - `test_enrichment_dotplot` — verify returns go.Figure
  - `test_enrichment_dotplot_empty` — verify handles empty DataFrame
  - `test_deconvolution_plot` — verify returns go.Figure
  
  **EnrichmentResult extension**:
  - `test_enrichment_result_backward_compatible` — verify old construction still works
  - `test_enrichment_result_new_fields` — verify new fields have correct defaults

  **Must NOT do**:
  - Do NOT require network access for tests (mock Enrichr/GSEA API calls)
  - Do NOT import streamlit in tests
  - Do NOT modify existing test files

  **Recommended Agent Profile**:
  - **Category**: `quick`
  - **Skills**: [`git-master`]

  **Parallelization**:
  - **Can Run In Parallel**: NO
  - **Parallel Group**: Wave 3 (after Task 7)
  - **Blocks**: None
  - **Blocked By**: Tasks 1-7

  **References**:
  - `tests/test_enhanced_viz.py` — Existing test patterns for visualization functions
  - `tests/test_de_bugfix.py` — Test patterns for DE analysis
  - `tests/test_demo_data.py` — Simple module test patterns
  - `conftest.py` — Existing test fixtures
  - All new module files (Tasks 1-6 outputs)

  **Acceptance Criteria**:
  ```bash
  # Run all tests
  pytest tests/test_new_features.py -v
  # Assert: All tests pass (0 failures)
  
  # Run full test suite to verify nothing broken
  pytest tests/ -v
  # Assert: All tests pass (existing + new)
  ```

  **Commit**: YES
  - Message: `test: add comprehensive tests for GSEA, interpretation, deconvolution, and advanced QC`
  - Files: `tests/test_new_features.py`
  - Pre-commit: `pytest tests/test_new_features.py -v`

---

## Commit Strategy

| After Task | Message | Files | Verification |
|------------|---------|-------|--------------|
| 1 | `feat(enrichment): add GSEA and 5 additional enrichment databases` | `pathway_enrichment.py` | `python -c "from pathway_enrichment import PathwayEnrichment"` |
| 2 | `feat: add interpretation engine with automated insights and methods text generation` | `interpretation_engine.py` | `python -c "from interpretation_engine import InterpretationEngine"` |
| 3 | `feat: add NNLS cell type deconvolution with skin tissue signature` | `cell_deconvolution.py` | `python -c "from cell_deconvolution import CellTypeDeconvolution"` |
| 4 | `feat: add advanced QC with outlier detection and batch effect assessment` | `advanced_qc.py` | `python -c "from advanced_qc import AdvancedQC"` |
| 5 | `feat(viz): add Plotly enrichment dotplot and deconvolution stacked bar chart` | `visualizations.py` | `python -c "from visualizations import create_enrichment_dotplot"` |
| 6 | `feat(export): extend EnrichmentResult with GSEA and additional database fields` | `export_engine.py` | `pytest tests/ -x -q` |
| 7 | `feat: integrate GSEA, interpretation, deconvolution, and advanced QC into Streamlit UI` | `rnaseq_analysis_platform.py` | `python -c "import rnaseq_analysis_platform"` |
| 8 | `test: add comprehensive tests for new features` | `tests/test_new_features.py` | `pytest tests/ -v` |

---

## Success Criteria

### Verification Commands
```bash
# All tests pass
pytest tests/ -v

# All imports work
python -c "
from pathway_enrichment import PathwayEnrichment
from interpretation_engine import InterpretationEngine
from cell_deconvolution import CellTypeDeconvolution
from advanced_qc import AdvancedQC
from visualizations import create_enrichment_dotplot, create_deconvolution_plot
from export_engine import EnrichmentResult
print('All imports OK')
"

# No matplotlib in UI code
grep -r "st.pyplot\|import matplotlib" rnaseq_analysis_platform.py && echo "FAIL: matplotlib found" || echo "PASS: no matplotlib"

# App starts without error
timeout 10 streamlit run rnaseq_analysis_platform.py --server.headless true 2>&1 | head -5
```

### Final Checklist
- [ ] All 8 original tabs work identically to before
- [ ] 4 new tabs appear: GSEA, Interpretation, Deconvolution, QC Advanced
- [ ] Enrichment tab enhanced with dot plots
- [ ] All visualizations use Plotly (zero matplotlib in UI)
- [ ] Full test suite passes: `pytest tests/ -v`
- [ ] Backward compatible: existing EnrichmentResult construction unbroken
- [ ] Pipeline errors in new features don't block existing analysis
