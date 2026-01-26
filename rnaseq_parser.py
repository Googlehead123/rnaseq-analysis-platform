"""
RNA-seq data parser module.

Handles multi-format RNA-seq data ingestion (CSV, TSV, Excel) with automatic:
- Data type detection (RAW_COUNTS, NORMALIZED, PRE_ANALYZED)
- Gene column detection
- Orientation detection and transposition to canonical shape (samples × genes)

Canonical output: samples × genes DataFrame (sample names as index, gene symbols as columns)
"""

from dataclasses import dataclass
from enum import Enum
from typing import Optional, List
import pandas as pd
import numpy as np
import re


# Gene column detection constants
KNOWN_SAMPLE_HEADERS = [
    "sample",
    "sample_id",
    "SampleID",
    "Sample_ID",
    "samplename",
    "Sample_Name",
    "samples",
    "id",
    "ID",
]

KNOWN_GENE_HEADERS = [
    "gene",
    "Gene",
    "GENE",
    "GeneSymbol",
    "gene_symbol",
    "gene_id",
    "SYMBOL",
    "GeneName",
    "gene_name",
]


class DataType(Enum):
    """RNA-seq data type classification."""

    RAW_COUNTS = "raw_counts"
    NORMALIZED = "normalized"
    PRE_ANALYZED = "pre_analyzed"


class ParserValidationError(Exception):
    """Raised when parsed data fails validation checks."""

    def __init__(self, message: str, details: Optional[dict] = None):
        self.message = message
        self.details = details or {}
        super().__init__(self.message)


@dataclass
class GeneColumnDetectionResult:
    """Result from detect_gene_column with explicit flags."""

    gene_column: Optional[str]  # Detected gene column name (None if not found/needed)
    candidates: List[str]  # Candidate columns if ambiguous
    first_col_is_sample_id: bool  # True if first col detected as sample ID (not gene)
    source: str  # How detection was made
    # source values: "known_gene_header:{col}", "heuristic_gene_col:{col}",
    #                "sample_id_detected", "ambiguous", "no_candidates"


@dataclass
class ParseResult:
    """Result from parsing RNA-seq data file.

    Field population rules by DataType:

    RAW_COUNTS:
        - expression_df: populated (samples × genes)
        - de_results_df: None
        - can_run_de: True
        - gene_column_source: how gene column was detected
        - needs_user_input: True if ambiguous gene column
        - gene_column_candidates: list if ambiguous, else []

    NORMALIZED:
        - expression_df: populated (samples × genes)
        - de_results_df: None
        - can_run_de: False (normalized data cannot be used for DE)
        - gene_column_source: how gene column was detected
        - needs_user_input: True if ambiguous gene column
        - gene_column_candidates: list if ambiguous, else []

    PRE_ANALYZED:
        - expression_df: None
        - de_results_df: populated (DE results table with canonical columns)
        - can_run_de: False (already analyzed)
        - gene_column_source: "pre_analyzed_gene_col"
        - needs_user_input: False
        - gene_column_candidates: []
    """

    # Core data (varies by DataType)
    expression_df: Optional[
        pd.DataFrame
    ]  # samples × genes matrix (RAW_COUNTS/NORMALIZED only)
    de_results_df: Optional[pd.DataFrame]  # DE results table (PRE_ANALYZED only)

    # Metadata
    data_type: DataType  # Detected data type
    can_run_de: bool  # True only if RAW_COUNTS
    warnings: List[str]  # Any warnings (e.g., "5 duplicate genes summed")
    dropped_columns: List[str]  # Non-numeric columns that were dropped
    gene_column_source: str  # How gene column was identified
    needs_user_input: bool  # True if parser couldn't auto-determine gene column
    gene_column_candidates: List[str]  # Candidate columns when needs_user_input=True


def looks_like_sample_names(values: List[str]) -> bool:
    """Check if values look like sample identifiers (not gene names)."""
    # Sample patterns: Sample_1, S1, Ctrl_1, Treatment_A, etc.
    # Characteristics: often contain underscores/numbers, sequential patterns
    sample_patterns = [
        r"^[Ss]ample[_-]?\d+$",  # Sample_1, Sample1
        r"^[Ss]\d+$",  # S1, S2
        r"^(Ctrl|Control|Treat|Treatment)[_-]?\d*",  # Ctrl_1, Treatment
        r"^[A-Za-z]+[_-]\d+$",  # Prefix_Number pattern
    ]
    match_count = sum(1 for v in values if any(re.match(p, v) for p in sample_patterns))
    return match_count >= len(values) * 0.5  # >50% match sample patterns


def detect_gene_column(df: pd.DataFrame) -> GeneColumnDetectionResult:
    """
    Detect gene column with confidence or return candidates.

    CRITICAL: Returns explicit flags to avoid ambiguity between "sample ID detected"
    and "no gene column candidates" cases.

    Returns:
        GeneColumnDetectionResult with gene_column, candidates, first_col_is_sample_id, source
    """
    first_col = df.columns[0]

    # Priority 0: Explicitly EXCLUDE known sample ID column headers
    if first_col.lower() in [h.lower() for h in KNOWN_SAMPLE_HEADERS]:
        return GeneColumnDetectionResult(
            gene_column=None,
            candidates=[],
            first_col_is_sample_id=True,
            source="sample_id_detected",
        )

    # Priority 1: Exact gene header match
    for col in df.columns:
        if col in KNOWN_GENE_HEADERS:
            return GeneColumnDetectionResult(
                gene_column=col,
                candidates=[],
                first_col_is_sample_id=False,
                source=f"known_gene_header:{col}",
            )

    # Priority 2: First column if non-numeric AND looks like gene names (NOT sample names)
    if df[first_col].dtype == "object":
        sample_values = df[first_col].dropna().head(10).astype(str).tolist()

        # Check if values look like sample names - if so, NOT a gene column
        if looks_like_sample_names(sample_values):
            return GeneColumnDetectionResult(
                gene_column=None,
                candidates=[],
                first_col_is_sample_id=True,
                source="sample_id_detected",
            )

        # Heuristic: gene names are typically 3-20 chars, alphanumeric
        looks_like_genes = all(
            3 <= len(v) <= 20 and v[0].isalpha() for v in sample_values
        )
        if looks_like_genes:
            return GeneColumnDetectionResult(
                gene_column=first_col,
                candidates=[],
                first_col_is_sample_id=False,
                source=f"heuristic_gene_col:{first_col}",
            )

    # Priority 3: Ambiguous - return all non-numeric column candidates
    non_numeric_cols = df.select_dtypes(exclude=[np.number]).columns.tolist()
    if non_numeric_cols:
        return GeneColumnDetectionResult(
            gene_column=None,
            candidates=non_numeric_cols,
            first_col_is_sample_id=False,
            source="ambiguous",
        )

    # No candidates at all - data might already be samples × genes with numeric-only columns
    return GeneColumnDetectionResult(
        gene_column=None,
        candidates=[],
        first_col_is_sample_id=False,
        source="no_candidates",
    )


def detect_orientation(
    df: pd.DataFrame, gene_column: Optional[str], first_col_is_sample_id: bool
) -> str:
    """
    Detect whether data is genes × samples or samples × genes.

    CRITICAL: Use gene_column detection result to inform orientation.

    Heuristics (in order):
    1. If first_col_is_sample_id=True (from detect_gene_column), data is ALREADY samples × genes
    2. If gene_column is provided/detected (non-None), genes are ROWS (genes × samples) → needs transpose
    3. If row count >> column count (e.g., 20000 rows × 10 cols), genes are ROWS → needs transpose
    4. If column count >> row count (e.g., 10 rows × 20000 cols), genes are COLUMNS → already canonical
    5. Default for small square-ish matrices: assume samples × genes (no transpose)

    Returns:
        "samples_as_rows" (already canonical, no transpose needed)
        "genes_as_rows" (needs transpose to become samples × genes)
    """
    # Heuristic 1: First column is sample ID → data is already samples × genes
    if first_col_is_sample_id:
        return "samples_as_rows"  # Already canonical, no transpose needed

    # Heuristic 2: Gene column detected → genes are rows, needs transpose
    if gene_column:
        return "genes_as_rows"

    # Heuristic 3: Shape ratio (typical RNA-seq: ~20k genes × 10-50 samples)
    n_rows, n_cols = df.shape
    if n_rows > n_cols * 10:  # Many more rows than columns → genes as rows
        return "genes_as_rows"
    elif n_cols > n_rows * 10:  # Many more columns than rows → already samples × genes
        return "samples_as_rows"

    # Default for small/square matrices: assume already samples × genes (no transpose)
    # This matches the common "samples as rows, genes as columns" export format
    return "samples_as_rows"


def is_close_to_integer(value: float, tolerance: float = 0.001) -> bool:
    """Check if value is within tolerance of nearest integer."""
    return abs(value - round(value)) <= tolerance


def are_counts_integer_like(df: pd.DataFrame) -> bool:
    """Check if all numeric values are close to integers (raw counts)."""
    numeric_values = df.select_dtypes(include=[np.number]).values.flatten()
    return all(is_close_to_integer(v) for v in numeric_values if not np.isnan(v))


def detect_data_type(df: pd.DataFrame) -> DataType:
    """
    Detect whether data is RAW_COUNTS, NORMALIZED, or PRE_ANALYZED.

    Detection logic:
    - RAW_COUNTS: All numeric values are integers (or close to integers within tolerance)
    - NORMALIZED: Contains floats with decimal places, typically 0-1000 range (TPM/FPKM)
    - PRE_ANALYZED: Contains columns like 'log2FoldChange', 'padj', 'pvalue'

    Returns:
        DataType enum value
    """
    # Check for PRE_ANALYZED first (has specific DE result columns)
    de_columns = ["log2FoldChange", "log2FC", "logFC", "padj", "FDR", "pvalue"]
    has_de_columns = any(col in df.columns for col in de_columns)
    if has_de_columns:
        return DataType.PRE_ANALYZED

    # Check if values are integer-like (RAW_COUNTS)
    if are_counts_integer_like(df):
        return DataType.RAW_COUNTS

    # Otherwise assume NORMALIZED (floats, TPM/FPKM/CPM)
    return DataType.NORMALIZED


def convert_to_canonical_shape(
    df: pd.DataFrame,
    gene_column: Optional[str],
    orientation: str,
    first_col_is_sample_id: bool,
) -> pd.DataFrame:
    """
    Convert input DataFrame to canonical shape: samples × genes.

    Canonical shape for DEResult/viz:
    - Index: sample names (e.g., ["Sample_1", "Sample_2", ...])
    - Columns: gene symbols (e.g., ["ACTB", "GAPDH", ...])

    Args:
        df: Raw parsed DataFrame
        gene_column: Column name containing gene identifiers (or None if not detected)
        orientation: "samples_as_rows" or "genes_as_rows"
        first_col_is_sample_id: True if first column is sample ID (from detect_gene_column)

    Returns:
        DataFrame in samples × genes format
    """
    if orientation == "samples_as_rows":
        # Data is ALREADY samples × genes (most common modern format)
        result = df.copy()

        # If first column is sample ID, set it as index
        if first_col_is_sample_id:
            first_col = df.columns[0]
            result = result.set_index(first_col)

        # Keep only numeric columns (drop Description, Length, etc.)
        numeric_cols = result.select_dtypes(include=[np.number]).columns
        dropped = [c for c in result.columns if c not in numeric_cols]
        result = result[numeric_cols]

        return result

    elif orientation == "genes_as_rows":
        # genes × samples format → needs transpose
        # Step 1: Set gene column as index
        if gene_column:
            result = df.set_index(gene_column)
        else:
            result = df.copy()  # Use existing index

        # Step 2: Keep only numeric columns (drop Description, Length, etc.)
        numeric_cols = result.select_dtypes(include=[np.number]).columns
        dropped = [c for c in result.columns if c not in numeric_cols]
        result = result[numeric_cols]

        # Step 3: Handle duplicate gene symbols (sum counts)
        if result.index.duplicated().any():
            n_dups = result.index.duplicated().sum()
            result = result.groupby(level=0).sum()
            # Warning added to ParseResult.warnings

        # Step 4: Transpose to samples × genes
        result = result.T

        return result

    else:
        raise ValueError(f"Unknown orientation: {orientation}")


def validate_for_de(df: pd.DataFrame) -> None:
    """
    Validate DataFrame is suitable for PyDESeq2 differential expression.

    Checks:
    - No negative values (count matrices cannot be negative)

    Raises:
        ParserValidationError: If validation fails
    """
    # Check for negative values
    numeric_cols = df.select_dtypes(include=[np.number])
    if (numeric_cols < 0).any().any():
        negative_count = (numeric_cols < 0).sum().sum()
        raise ParserValidationError(
            "Count matrices cannot contain negative values",
            details={"negative_count": int(negative_count)},
        )


def parse_csv(file_path: str) -> pd.DataFrame:
    """
    Parse CSV or TSV file with automatic delimiter detection.

    Args:
        file_path: Path to CSV/TSV file

    Returns:
        Parsed DataFrame (raw, before processing)

    Raises:
        ParserValidationError: If file cannot be parsed
    """
    try:
        # Try comma delimiter first
        df = pd.read_csv(file_path, sep=",")

        # If only one column, try tab delimiter
        if len(df.columns) == 1:
            df = pd.read_csv(file_path, sep="\t")

        # Validate minimum requirements
        if df.empty:
            raise ParserValidationError("File is empty")

        if len(df.columns) < 2:
            raise ParserValidationError(
                "File must have at least 2 columns (genes + samples)",
                details={"columns": len(df.columns)},
            )

        return df

    except pd.errors.EmptyDataError:
        raise ParserValidationError("File is empty")
    except pd.errors.ParserError as e:
        raise ParserValidationError(f"Failed to parse CSV: {str(e)}")
    except FileNotFoundError:
        raise ParserValidationError(f"File not found: {file_path}")
    except Exception as e:
        # Don't re-raise ParserValidationError
        if isinstance(e, ParserValidationError):
            raise
        raise ParserValidationError(f"Error reading file: {str(e)}")
