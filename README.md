# RNA-seq Analysis Platform

A user-friendly Streamlit-based platform for RNA-seq data analysis, designed for wet lab researchers in cosmetics and dermatological science. Eliminates the technical complexity of RNA-seq analysis while providing powerful differential expression, pathway enrichment, and visualization capabilities.

## Features

### 🧬 Data Loading & Validation
- **Multi-format support**: CSV, TSV, Excel (.xlsx)
- **Auto-detection**: Automatically identifies data type (raw counts, normalized, pre-analyzed)
- **Smart parsing**: Detects gene columns, sample orientation, and matrix structure
- **Validation**: Ensures data quality before analysis

### 📊 Differential Expression Analysis
- **PyDESeq2 integration**: Industry-standard statistical analysis
- **Multi-comparison support**: Analyze multiple treatment groups simultaneously
- **Adaptive thresholds**: Automatic fallback for low gene counts
- **Graceful error handling**: Failed comparisons don't block others

### 🔬 Pathway Enrichment
- **GO Biological Process**: Query Gene Ontology via Enrichr API
- **KEGG Pathways**: Identify enriched metabolic and signaling pathways
- **Adaptive gene selection**: Smart thresholding (padj < 0.05, fallback to 0.1)
- **Top 20 pathways**: Focused results for biological interpretation

### 📈 Interactive Visualizations
- **Volcano Plot**: Interactive scatter with significance coloring
- **Clustered Heatmap**: Hierarchical clustering with condition grouping
- **PCA Plot**: Sample clustering and quality control
- **Gene Panels**: Curated dermatology-specific gene sets

### 🎨 Dermatology Gene Panels
Six curated panels for skin biology research:
- **Anti-aging**: Collagen synthesis, ECM remodeling (COL1A1, ELN, MMPs)
- **Skin Barrier**: Stratum corneum, tight junctions (FLG, LOR, IVL)
- **Anti-inflammation**: NF-κB pathway cytokines (IL1A, IL6, TNFA)
- **Whitening/Melanogenesis**: Melanin synthesis (MITF, TYR, TYRP1)
- **Sebum Regulation**: Lipogenesis factors (SREBF1, PPARG, FASN)
- **Wound Healing**: Growth factors, keratins (VEGFA, FGF2, TGFB1)

### 💾 Export Capabilities
- **Excel Workbooks**: Multi-sheet with DE results, enrichment, settings
- **High-Resolution Images**: PNG (300 DPI), SVG, PDF (vector)
- **PDF Reports**: Comprehensive analysis reports with embedded figures
- **Jupyter Notebooks**: Templates for advanced customization

## Quick Start

### Installation

```bash
# Clone or download the repository
cd rnaseq-analysis-platform

# Create virtual environment (Python 3.10 or 3.11 required)
python3.10 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install Playwright for browser testing (optional)
playwright install chromium
```

### Running the Application

```bash
# Launch Streamlit app
streamlit run rnaseq_analysis_platform.py

# Or use Jupyter notebooks for advanced analysis
jupyter lab
# Open notebooks/rnaseq_analysis_template.ipynb
```

### Basic Workflow

1. **Upload Data**: CSV/Excel file with gene counts
2. **Assign Metadata**: Map samples to experimental conditions
3. **Select Comparisons**: Choose test vs reference groups
4. **Run Analysis**: Differential expression with PyDESeq2
5. **Explore Results**: Interactive visualizations and pathway enrichment
6. **Export**: Download Excel, PDF, or high-res images

## Supported Data Formats

### Input Files

**Raw Count Matrix** (required for DE analysis):
```
Gene      Sample_1  Sample_2  Sample_3  Sample_4
COL1A1    1523      1834      2145      1923
IL6       234       456       789       567
TYR       89        123       156       134
```

**Normalized Data** (visualization only):
```
Gene      Sample_1  Sample_2  Sample_3  Sample_4
COL1A1    12.5      13.2      14.1      13.8
IL6       8.3       9.1       10.2      9.5
TYR       6.7       7.2       7.8       7.4
```

**Pre-Analyzed Results** (from external tools):
```
Gene      log2FoldChange  padj        baseMean
COL1A1    2.34           0.0001      1500.2
IL6       -1.89          0.0023      450.8
TYR       1.45           0.0456      120.3
```

### Supported Formats
- **CSV**: Comma-separated (`.csv`)
- **TSV**: Tab-separated (`.tsv`, `.txt`)
- **Excel**: Multi-sheet workbooks (`.xlsx`)

### Matrix Orientations
The parser automatically detects and converts:
- **Genes as rows, samples as columns** → Canonical format
- **Samples as rows, genes as columns** → Auto-transposed

## Requirements

### System Requirements
- **Python**: 3.10 or 3.11 (PyDESeq2 compatibility)
- **RAM**: 4GB minimum, 8GB recommended
- **Storage**: 500MB for dependencies

### Python Dependencies
All dependencies listed in `requirements.txt`:
- **Analysis**: `pydeseq2==0.4.11`, `gseapy`, `scikit-learn`
- **Visualization**: `plotly`, `kaleido`
- **Export**: `openpyxl`, `reportlab`, `pypdf`
- **UI**: `streamlit`, `pandas`, `numpy`

## Project Structure

```
rnaseq-analysis-platform/
├── rnaseq_parser.py           # Multi-format data parser
├── de_analysis.py             # PyDESeq2 DE engine
├── pathway_enrichment.py      # GSEApy enrichment
├── visualizations.py          # Plotly visualizations
├── gene_panels.py             # Gene panel analyzer
├── export_engine.py           # Excel/PDF/image export
├── rnaseq_analysis_platform.py # Main Streamlit app
├── config/
│   └── gene_panels.yaml       # Dermatology gene panels
├── notebooks/
│   ├── rnaseq_analysis_template.ipynb
│   └── advanced_exploration.ipynb
├── tests/
│   └── data/                  # Test datasets
└── requirements.txt
```

## Usage Examples

### Python API (Jupyter Notebooks)

```python
from rnaseq_parser import RNASeqParser
from de_analysis import DEAnalysisEngine
from visualizations import create_volcano_plot

# Load data
parser = RNASeqParser()
result = parser.parse("data/counts.csv")

# Run DE analysis
engine = DEAnalysisEngine()
de_results = engine.run_all_comparisons(
    result.expression_df,
    metadata_df,
    comparisons=[("Treatment", "Control")],
    design_factor="condition"
)

# Visualize
de_result = de_results[("Treatment", "Control")]
fig = create_volcano_plot(de_result.results_df)
fig.show()
```

### Command Line

```bash
# Run analysis pipeline
python -c "
from rnaseq_parser import RNASeqParser
parser = RNASeqParser()
result = parser.parse('tests/data/sample_counts.csv')
print(f'Loaded {result.expression_df.shape[0]} samples, {result.expression_df.shape[1]} genes')
"
```

## Advanced Features

### Custom Gene Sets

Define your own gene panels in `config/gene_panels.yaml`:

```yaml
panels:
  My_Custom_Panel:
    description: "Custom pathway genes"
    genes:
      - GENE1
      - GENE2
      - GENE3
```

### Multi-Comparison Analysis

Analyze multiple treatment groups:

```python
comparisons = [
    ("Treatment_A", "Control"),
    ("Treatment_B", "Control"),
    ("Treatment_B", "Treatment_A")
]

de_results = engine.run_all_comparisons(
    counts_df, metadata_df, 
    comparisons=comparisons,
    design_factor="condition"
)
```

### Publication-Quality Figures

Customize Plotly figures for publication:

```python
fig = create_volcano_plot(de_result.results_df)
fig.update_layout(
    font={'family': 'Arial', 'size': 14},
    width=800,
    height=600
)

# Export high-resolution
from export_engine import ExportEngine
engine = ExportEngine()
engine.export_figure(fig, "volcano.png", format="png", scale=3)  # 300 DPI
```

## Troubleshooting

### Common Issues

**Import Error: PyDESeq2**
```bash
# Ensure Python 3.10 or 3.11
python --version

# Reinstall in clean environment
pip install --force-reinstall pydeseq2==0.4.11
```

**Enrichment API Timeout**
```
# Enrichr requires internet connection
# If offline, enrichment will return empty results with error message
```

**Memory Error with Large Datasets**
```python
# For >50,000 genes, filter low-expression genes first
counts_filtered = counts_df[counts_df.mean(axis=0) > 10]
```

### Getting Help

1. Check Jupyter notebook templates for examples
2. Review test data in `tests/data/` for format reference
3. Read module docstrings: `help(RNASeqParser)`
4. Check `.sisyphus/notepads/rnaseq-analysis-platform/learnings.md` for implementation notes

## Development

### Running Tests

```bash
# Install test dependencies
pip install pytest pytest-cov

# Run all tests
pytest tests/ -v

# With coverage
pytest tests/ --cov=. --cov-report=term-missing
```

### Code Quality

```bash
# Type checking
python -m py_compile *.py

# Import verification
python -c "from de_analysis import DEAnalysisEngine; print('OK')"
```

## Citation

If you use this platform in your research, please cite:

```
RNA-seq Analysis Platform for Dermatological Research
https://github.com/yourusername/rnaseq-analysis-platform
```

## License

This project is provided as-is for research purposes.

## Acknowledgments

- **PyDESeq2**: Python implementation of DESeq2
- **GSEApy**: Enrichr API integration
- **Plotly**: Interactive visualizations
- **Streamlit**: Web application framework

## Version History

### v1.0.0 (2026-01-26)
- Initial release
- Complete analysis pipeline (parser, DE, enrichment, visualization)
- 6 curated dermatology gene panels
- Multi-format export (Excel, PDF, images)
- Jupyter notebook templates

## Contact

For questions or issues, please open an issue on GitHub or contact the development team.

---

**Built for wet lab researchers by computational biologists** 🧬
