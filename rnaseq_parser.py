# pyright: reportMissingTypeArgument=false
from __future__ import annotations

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
from os import PathLike
from typing import Optional, List, Tuple, Any, Union
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

# PRE_ANALYZED column aliases
COLUMN_ALIASES = {
    "gene": ["gene", "Gene", "GENE", "GeneSymbol", "gene_symbol", "gene_id", "SYMBOL"],
    "log2FoldChange": ["log2FoldChange", "log2FC", "logFC", "log2_fold_change", "lfc"],
    "padj": ["padj", "FDR", "fdr", "adj.P.Val", "adjusted_pvalue", "q_value", "qvalue"],
    "baseMean": ["baseMean", "AveExpr", "mean_expression", "average_expression"],
    "pvalue": ["pvalue", "P.Value", "PValue", "raw_pvalue", "p_value"],
}

# Pattern-based DE column detection for non-standard formats
DE_COLUMN_PATTERNS = {
    "log2FoldChange": [r"^.+\.fc$", r"^.+\.logFC$", r"^.+\.log2FC$"],
    "padj": [r"^.+\.bh\.pval$", r"^.+\.adj\.pval$", r"^.+\.FDR$", r"^.+\.qvalue$"],
    "pvalue": [
        r"^.+\.raw\.pval$",
        r"^.+\.PValue$",
    ],  # Removed generic .pval$ to avoid matching .bh.pval
    "baseMean": [r"^.+\.baseMean$", r"^.+\.AveExpr$"],
}

# Sample column detection patterns
COUNT_COLUMN_PATTERN = r"^(.+)_(.+)_(.+)_Read_Count$"
FPKM_COLUMN_PATTERN = r"^(.+)_(.+)_(.+)_FPKM$"
TPM_COLUMN_PATTERN = r"^(.+)_(.+)_(.+)_TPM$"


class DataType(Enum):
    """RNA-seq data type classification."""

    RAW_COUNTS = "raw_counts"
    NORMALIZED = "normalized"
    PRE_ANALYZED = "pre_analyzed"


class ParserValidationError(Exception):
    """Raised when parsed data fails validation checks."""

    def __init__(self, message: str, details: Optional[dict[str, Any]] = None):
        self.message: str = message
        self.details: dict[str, Any] = details or {}
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
class SampleColumnInfo:
    """Result from detect_sample_columns with sample/condition information."""

    count_columns: List[str]
    fpkm_columns: List[str]
    tpm_columns: List[str]
    conditions_detected: set[str]
    sample_to_condition: dict[str, str]


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
    normalized_df: Optional[
        pd.DataFrame
    ]  # samples × genes with TPM/FPKM (NORMALIZED or multi-format)
    de_results_df: Optional[pd.DataFrame]  # DE results table (PRE_ANALYZED only)

    # Metadata
    data_type: DataType  # Detected data type
    can_run_de: bool  # True only if RAW_COUNTS
    warnings: List[str]  # Any warnings (e.g., "5 duplicate genes summed")
    dropped_columns: List[str]  # Non-numeric columns that were dropped
    gene_column_source: str  # How gene column was identified
    needs_user_input: bool  # True if parser couldn't auto-determine gene column
    gene_column_candidates: List[str]  # Candidate columns when needs_user_input=True
    data_types_detected: List[
        DataType
    ]  # All data types found in file (empty for single-type)


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
    first_col_label = df.columns[0]
    first_col = str(first_col_label)

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
        col_str = str(col)
        if col_str in KNOWN_GENE_HEADERS:
            return GeneColumnDetectionResult(
                gene_column=col_str,
                candidates=[],
                first_col_is_sample_id=False,
                source=f"known_gene_header:{col_str}",
            )

    # Priority 2: First column if non-numeric AND looks like gene names (NOT sample names)
    if df[first_col_label].dtype == "object":
        sample_values = df[first_col_label].dropna().head(10).astype(str).tolist()

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


def detect_de_columns(df: pd.DataFrame) -> Tuple[dict[str, str], Optional[str]]:
    """
    Detect DE result columns using pattern matching or standard aliases.

    Tries pattern matching first (for non-standard formats like Bt10U/none.fc),
    then falls back to COLUMN_ALIASES for standard column names.

    Returns:
        Tuple of (column_mapping, comparison_name) where:
        - column_mapping: Dict mapping canonical names to actual column names
          e.g., {"log2FoldChange": "Bt10U/none.fc", "padj": "Bt10U/none.bh.pval"}
        - comparison_name: Extracted comparison name (e.g., "Bt10U/none") or None
    """
    column_mapping = {}
    comparison_name = None

    for canonical_name, patterns in DE_COLUMN_PATTERNS.items():
        for col in df.columns:
            for pattern in patterns:
                if re.match(pattern, col, re.IGNORECASE):
                    column_mapping[canonical_name] = col

                    if comparison_name is None:
                        if canonical_name == "log2FoldChange":
                            comparison_name = re.sub(
                                r"\.fc$|\.logFC$|\.log2FC$",
                                "",
                                col,
                                flags=re.IGNORECASE,
                            )
                        elif canonical_name == "padj":
                            comparison_name = re.sub(
                                r"\.bh\.pval$|\.adj\.pval$|\.FDR$|\.qvalue$",
                                "",
                                col,
                                flags=re.IGNORECASE,
                            )
                        elif canonical_name == "pvalue":
                            comparison_name = re.sub(
                                r"\.raw\.pval$|\.pval$|\.PValue$",
                                "",
                                col,
                                flags=re.IGNORECASE,
                            )
                    break

    if column_mapping:
        return column_mapping, comparison_name

    for canonical_name, aliases in COLUMN_ALIASES.items():
        for col in df.columns:
            if col in aliases:
                column_mapping[canonical_name] = col

    return column_mapping, None


def detect_sample_columns(df: pd.DataFrame) -> SampleColumnInfo:
    """
    Detect sample columns and extract conditions from column names.

    Matches columns against patterns: {date}_{condition}_{timepoint}_{metric}
    - Read_Count: raw count data
    - FPKM: normalized expression
    - TPM: normalized expression

    Returns:
        SampleColumnInfo with detected columns, conditions, and mappings
    """
    count_columns = []
    fpkm_columns = []
    tpm_columns = []
    conditions_detected = set()
    sample_to_condition = {}

    for col in df.columns:
        match = re.match(COUNT_COLUMN_PATTERN, col)
        if match:
            count_columns.append(col)
            condition = match.group(2)
            conditions_detected.add(condition)
            sample_to_condition[col] = condition
            continue

        match = re.match(FPKM_COLUMN_PATTERN, col)
        if match:
            fpkm_columns.append(col)
            condition = match.group(2)
            conditions_detected.add(condition)
            sample_to_condition[col] = condition
            continue

        match = re.match(TPM_COLUMN_PATTERN, col)
        if match:
            tpm_columns.append(col)
            condition = match.group(2)
            conditions_detected.add(condition)
            sample_to_condition[col] = condition

    return SampleColumnInfo(
        count_columns=count_columns,
        fpkm_columns=fpkm_columns,
        tpm_columns=tpm_columns,
        conditions_detected=conditions_detected,
        sample_to_condition=sample_to_condition,
    )


def _is_multiformat_excel(df: pd.DataFrame) -> bool:
    """Check if Excel file contains multiple data types."""
    has_de = len(detect_de_columns(df)[0]) > 0
    sample_info = detect_sample_columns(df)
    has_samples = (
        len(sample_info.count_columns) > 0
        or len(sample_info.fpkm_columns) > 0
        or len(sample_info.tpm_columns) > 0
    )
    return has_de and has_samples


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
        numeric_cols = list(result.select_dtypes(include=[np.number]).columns)
        dropped = [c for c in result.columns if c not in numeric_cols]
        result = result.loc[:, numeric_cols]

        return result

    elif orientation == "genes_as_rows":
        # genes × samples format → needs transpose
        # Step 1: Set gene column as index
        if gene_column:
            result = df.set_index(gene_column)
        else:
            result = df.copy()  # Use existing index

        # Step 2: Keep only numeric columns (drop Description, Length, etc.)
        numeric_cols = list(result.select_dtypes(include=[np.number]).columns)
        dropped = [c for c in result.columns if c not in numeric_cols]
        result = result.loc[:, numeric_cols]

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
            f"Count matrices cannot contain negative values. Found {int(negative_count)} negative values. "
            f"Suggestion: Check if your data has been log-transformed or normalized. "
            f"PyDESeq2 requires raw count data (non-negative integers). "
            f"If you have normalized data, use the 'normalized' data type instead.",
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
            raise ParserValidationError(
                "File is empty. Ensure your file contains data rows with gene names and sample columns."
            )

        if len(df.columns) < 2:
            raise ParserValidationError(
                f"File must have at least 2 columns (genes + samples), but found {len(df.columns)}. "
                f"Suggestion: Check if your file uses the correct delimiter (comma for CSV, tab for TSV). "
                f"If using Excel, export as CSV first.",
                details={"columns": len(df.columns)},
            )

        return df

    except pd.errors.EmptyDataError:
        raise ParserValidationError(
            "File is empty or contains no readable data. "
            "Ensure your file has at least one row of data with gene names and sample values."
        )
    except pd.errors.ParserError as e:
        raise ParserValidationError(
            f"Failed to parse CSV file: {str(e)}. "
            f"Suggestion: Verify the file format is valid CSV/TSV. "
            f"Check for special characters or encoding issues. Try opening in a text editor to inspect."
        )
    except FileNotFoundError:
        raise ParserValidationError(
            f"File not found: {file_path}. "
            f"Suggestion: Check the file path is correct and the file exists in the specified location."
        )
    except Exception as e:
        # Don't re-raise ParserValidationError
        if isinstance(e, ParserValidationError):
            raise
        raise ParserValidationError(
            f"Error reading file: {str(e)}. "
            f"Suggestion: Ensure the file is not corrupted and is in a supported format (CSV, TSV, or Excel)."
        )


def parse_excel(file_path: str, sheet_name: Optional[str] = None) -> pd.DataFrame:
    """Parse Excel file (.xlsx)."""
    try:
        excel_file = pd.ExcelFile(file_path)

        if sheet_name:
            if sheet_name not in excel_file.sheet_names:
                available = ", ".join(excel_file.sheet_names[:5])
                raise ParserValidationError(
                    f"Sheet '{sheet_name}' not found. Available sheets: {available}{'...' if len(excel_file.sheet_names) > 5 else ''}. "
                    f"Suggestion: Check the sheet name spelling or select a different sheet.",
                    details={"available_sheets": excel_file.sheet_names},
                )
            df = pd.read_excel(file_path, sheet_name=sheet_name)
        else:
            df = pd.read_excel(file_path, sheet_name=0)

        if df.empty:
            raise ParserValidationError(
                "Sheet is empty. Ensure the selected sheet contains data rows with gene names and sample columns."
            )

        if len(df.columns) < 2:
            raise ParserValidationError(
                f"Sheet must have at least 2 columns (genes + samples), but found {len(df.columns)}. "
                f"Suggestion: Verify you selected the correct sheet with count data, not a metadata sheet.",
                details={"columns": len(df.columns)},
            )

        # Check if looks like metadata sheet
        if len(df) <= 20:
            non_numeric_pct = len(df.select_dtypes(exclude=[np.number]).columns) / len(
                df.columns
            )
            if non_numeric_pct > 0.5:
                raise ParserValidationError(
                    "First sheet appears to be metadata (mostly text), not count data. "
                    "Suggestion: Select a different sheet containing the count matrix with gene names and numeric sample columns."
                )

        return df

    except FileNotFoundError:
        raise ParserValidationError(
            f"File not found: {file_path}. "
            f"Suggestion: Check the file path is correct and the file exists."
        )
    except Exception as e:
        if isinstance(e, ParserValidationError):
            raise
        raise ParserValidationError(
            f"Error reading Excel file: {str(e)}. "
            f"Suggestion: Ensure the file is a valid Excel (.xlsx) file and not corrupted. "
            f"Try opening it in Excel to verify it's readable."
        )


def normalize_de_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Rename columns to canonical names for downstream compatibility."""
    df = df.copy()
    for canonical, aliases in COLUMN_ALIASES.items():
        for alias in aliases:
            if alias in df.columns and canonical not in df.columns:
                df = df.rename(columns={alias: canonical})
                break

    required = ["gene", "log2FoldChange", "padj"]
    missing = [col for col in required if col not in df.columns]
    if missing:
        available_str = ", ".join(df.columns.tolist()[:10])
        extra = f"... ({len(df.columns) - 10} more)" if len(df.columns) > 10 else ""
        raise ParserValidationError(
            f"Pre-analyzed file missing required columns: {', '.join(missing)}. "
            f"Found columns: {available_str}{extra}. "
            f"Suggestion: Ensure your file contains 'gene' (gene names), 'log2FoldChange' (or 'log2FC', 'logFC'), "
            f"and 'padj' (or 'FDR', 'adjusted_pvalue'). Check column names for typos or alternative naming conventions.",
            details={"missing": missing, "available": df.columns.tolist()},
        )

    return df


class RNASeqParser:
    """RNA-seq data parser with automatic format detection."""

    def _is_multiformat_excel(self, df: pd.DataFrame) -> bool:
        """Check if DataFrame contains multiple data types."""
        de_mapping, _ = detect_de_columns(df)
        sample_info = detect_sample_columns(df)
        has_de = len(de_mapping) > 0
        has_samples = (
            len(sample_info.count_columns) > 0
            or len(sample_info.fpkm_columns) > 0
            or len(sample_info.tpm_columns) > 0
        )
        return has_de and has_samples

    def _find_gene_column(self, df: pd.DataFrame) -> Optional[str]:
        """Find gene column in DataFrame."""
        gene_candidates = [
            "Gene_Symbol",
            "Gene",
            "gene",
            "GeneSymbol",
            "gene_symbol",
        ]
        for col in gene_candidates:
            if col in df.columns:
                return col
        return None

    def parse_multiformat(self, file_path: Union[str, PathLike[str]]) -> ParseResult:
        """Parse multi-format Excel with DE results + counts + normalized."""
        file_path = str(file_path)
        df = parse_excel(file_path)

        de_mapping, comparison_name = detect_de_columns(df)
        sample_info = detect_sample_columns(df)

        data_types = []
        warnings = []

        gene_col = self._find_gene_column(df)
        if gene_col is None:
            warnings.append("Gene column not found for multiformat extraction")

        de_results_df = None
        if de_mapping and gene_col:
            de_cols = [gene_col] + list(de_mapping.values())
            de_df: pd.DataFrame = pd.DataFrame(df.loc[:, de_cols]).copy()
            de_df = de_df.rename(columns=dict({v: k for k, v in de_mapping.items()}))
            de_df = de_df.rename(columns={gene_col: "gene"})
            de_results_df = de_df
            data_types.append(DataType.PRE_ANALYZED)

        expression_df = None
        if sample_info.count_columns and gene_col:
            count_df = pd.DataFrame(
                df.loc[:, [gene_col] + sample_info.count_columns]
            ).copy()
            count_df = count_df.set_index(gene_col)
            expression_df = count_df.T
            data_types.append(DataType.RAW_COUNTS)

        normalized_df = None
        if sample_info.tpm_columns:
            norm_cols = sample_info.tpm_columns
        elif sample_info.fpkm_columns:
            norm_cols = sample_info.fpkm_columns
        else:
            norm_cols = []

        if norm_cols and gene_col:
            norm_df_raw = pd.DataFrame(df.loc[:, [gene_col] + norm_cols]).copy()
            norm_df_raw = norm_df_raw.set_index(gene_col)
            normalized_df = norm_df_raw.T
            data_types.append(DataType.NORMALIZED)

        return ParseResult(
            expression_df=expression_df,
            normalized_df=normalized_df,
            de_results_df=de_results_df,
            data_type=data_types[0] if data_types else DataType.RAW_COUNTS,
            can_run_de=DataType.RAW_COUNTS in data_types,
            warnings=warnings,
            dropped_columns=[],
            gene_column_source="multiformat_detection",
            needs_user_input=False,
            gene_column_candidates=[],
            data_types_detected=data_types,
        )

    def parse(self, file_path: str, gene_column: Optional[str] = None) -> ParseResult:
        """Parse RNA-seq data file."""
        # Convert Path to string if needed
        file_path = str(file_path)

        # Detect format and parse
        if file_path.endswith(".xlsx") or file_path.endswith(".xls"):
            df = parse_excel(file_path)
            if self._is_multiformat_excel(df):
                return self.parse_multiformat(file_path)
        else:
            df = parse_csv(file_path)

        # Detect data type
        data_type = detect_data_type(df)

        # Handle PRE_ANALYZED
        if data_type == DataType.PRE_ANALYZED:
            df = normalize_de_columns(df)
            return ParseResult(
                expression_df=None,
                normalized_df=None,
                de_results_df=df,
                data_type=DataType.PRE_ANALYZED,
                can_run_de=False,
                warnings=[],
                dropped_columns=[],
                gene_column_source="pre_analyzed",
                needs_user_input=False,
                gene_column_candidates=[],
                data_types_detected=[],
            )

        # Detect gene column
        if gene_column:
            gene_det = GeneColumnDetectionResult(
                gene_column=gene_column,
                candidates=[],
                first_col_is_sample_id=False,
                source=f"user_specified:{gene_column}",
            )
        else:
            gene_det = detect_gene_column(df)

        # Check if ambiguous
        if gene_det.candidates and not gene_column:
            return ParseResult(
                expression_df=None,
                normalized_df=None,
                de_results_df=None,
                data_type=data_type,
                can_run_de=False,
                warnings=[],
                dropped_columns=[],
                gene_column_source=gene_det.source,
                needs_user_input=True,
                gene_column_candidates=gene_det.candidates,
                data_types_detected=[],
            )

        # Detect orientation and convert
        orientation = detect_orientation(
            df, gene_det.gene_column, gene_det.first_col_is_sample_id
        )
        expression_df = convert_to_canonical_shape(
            df, gene_det.gene_column, orientation, gene_det.first_col_is_sample_id
        )

        # Validate if RAW_COUNTS
        warnings = []
        if data_type == DataType.RAW_COUNTS:
            validate_for_de(expression_df)

        return ParseResult(
            expression_df=expression_df,
            normalized_df=None,
            de_results_df=None,
            data_type=data_type,
            can_run_de=(data_type == DataType.RAW_COUNTS),
            warnings=warnings,
            dropped_columns=[],
            gene_column_source=gene_det.source,
            needs_user_input=False,
            gene_column_candidates=[],
            data_types_detected=[],
        )
