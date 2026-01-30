"""Utility for merging multiple RNA-seq count files into a single expression matrix."""

import pandas as pd
import re
from typing import List, Dict, Optional
from pathlib import Path
from rnaseq_parser import ParseResult, DataType


def validate_merge_compatibility(results: List[ParseResult]) -> Dict:
    """Validate that multiple ParseResults can be merged.
    
    Returns dict with keys: compatible (bool), issues (list), report (dict with stats).
    """
    report = {
        "n_files": len(results),
        "data_types": [],
        "samples_per_file": [],
        "genes_per_file": [],
        "gene_overlap_pct": 0.0,
    }
    issues = []
    
    valid_results = [r for r in results if r.expression_df is not None]
    if len(valid_results) < 2:
        return {"compatible": len(valid_results) == 1, "issues": ["Need at least 2 files with expression data to merge"], "report": report}
    
    # Check DataType consistency
    data_types = set(r.data_type for r in valid_results)
    report["data_types"] = [dt.value for dt in data_types]
    if len(data_types) > 1:
        issues.append(f"Mixed data types: {', '.join(dt.value for dt in data_types)}. All files must be the same type.")
    
    # Gene overlap
    gene_sets = [set(r.expression_df.columns) for r in valid_results]
    all_genes = set.union(*gene_sets)
    shared_genes = set.intersection(*gene_sets)
    overlap_pct = len(shared_genes) / len(all_genes) * 100 if all_genes else 0
    report["gene_overlap_pct"] = round(overlap_pct, 1)
    report["shared_genes"] = len(shared_genes)
    report["total_genes"] = len(all_genes)
    
    if len(shared_genes) == 0:
        issues.append("No shared genes across files. Cannot merge.")
    elif overlap_pct < 50:
        issues.append(f"Low gene overlap ({overlap_pct:.1f}%). Merge will keep only {len(shared_genes)} of {len(all_genes)} genes.")
    
    for i, r in enumerate(valid_results):
        report["samples_per_file"].append(r.expression_df.shape[0])
        report["genes_per_file"].append(r.expression_df.shape[1])
    
    return {"compatible": len(issues) == 0 or (len(data_types) == 1 and len(shared_genes) > 0), "issues": issues, "report": report}


def _sanitize_sample_name(name: str) -> str:
    """Sanitize sample name: replace special chars with underscore."""
    return re.sub(r'[^a-zA-Z0-9_]', '_', str(name)).strip('_')


def _resolve_duplicate_samples(dfs: List[pd.DataFrame], filenames: List[str]) -> List[pd.DataFrame]:
    """Handle duplicate sample names across files by suffixing with filename stem."""
    all_samples = []
    for df in dfs:
        all_samples.extend(df.index.tolist())
    
    # Check for duplicates
    seen = set()
    has_duplicates = False
    for s in all_samples:
        if s in seen:
            has_duplicates = True
            break
        seen.add(s)
    
    if not has_duplicates:
        return dfs
    
    # Suffix all sample names with filename stem
    new_dfs = []
    for df, fname in zip(dfs, filenames):
        stem = Path(fname).stem
        stem = _sanitize_sample_name(stem)
        new_index = [f"{_sanitize_sample_name(s)}_{stem}" for s in df.index]
        new_df = df.copy()
        new_df.index = new_index
        new_dfs.append(new_df)
    return new_dfs


def merge_parse_results(results: List[ParseResult], filenames: List[str]) -> ParseResult:
    """Merge multiple ParseResults into a single combined result.
    
    Args:
        results: List of ParseResult objects from individual file parsing
        filenames: Original filenames (used for sample name deduplication)
    
    Returns:
        Single merged ParseResult with combined expression_df (inner join on genes)
    
    Raises:
        ValueError: If data types are mixed or no gene overlap exists
    """
    # Filter valid results
    valid_pairs = [(r, f) for r, f in zip(results, filenames) if r.expression_df is not None]
    
    if len(valid_pairs) == 0:
        raise ValueError("No files contain expression data")
    
    if len(valid_pairs) == 1:
        return valid_pairs[0][0]
    
    valid_results = [p[0] for p in valid_pairs]
    valid_filenames = [p[1] for p in valid_pairs]
    
    # Validate compatibility
    validation = validate_merge_compatibility(valid_results)
    
    # Check data type consistency
    data_types = set(r.data_type for r in valid_results)
    if len(data_types) > 1:
        raise ValueError(f"Cannot merge files with different data types: {', '.join(dt.value for dt in data_types)}")
    
    data_type = valid_results[0].data_type
    
    # Get expression DataFrames
    dfs = [r.expression_df for r in valid_results]
    
    # Resolve duplicate sample names
    dfs = _resolve_duplicate_samples(dfs, valid_filenames)
    
    # Inner join on genes (columns)
    shared_genes = set.intersection(*[set(df.columns) for df in dfs])
    if len(shared_genes) == 0:
        raise ValueError("No shared genes across files. Cannot merge.")
    
    shared_genes_sorted = sorted(shared_genes)
    trimmed_dfs = [df[shared_genes_sorted] for df in dfs]
    
    # Concatenate samples (axis=0)
    merged_df = pd.concat(trimmed_dfs, axis=0)
    
    # Build warnings
    all_warnings = []
    for r in valid_results:
        all_warnings.extend(r.warnings)
    
    total_genes = len(set.union(*[set(df.columns) for df in [r.expression_df for r in valid_results]]))
    genes_dropped = total_genes - len(shared_genes)
    if genes_dropped > 0:
        overlap_pct = len(shared_genes) / total_genes * 100
        all_warnings.append(
            f"Merged {len(valid_results)} files: kept {len(shared_genes)} shared genes "
            f"({overlap_pct:.1f}% overlap), dropped {genes_dropped} genes not shared across all files."
        )
    
    all_warnings.append(
        f"Merged expression matrix: {merged_df.shape[0]} samples × {merged_df.shape[1]} genes "
        f"from {len(valid_results)} files."
    )
    
    return ParseResult(
        expression_df=merged_df,
        normalized_df=None,
        de_results_df=None,
        data_type=data_type,
        can_run_de=(data_type == DataType.RAW_COUNTS),
        warnings=all_warnings,
        dropped_columns=[],
        gene_column_source="merged",
        needs_user_input=False,
        gene_column_candidates=[],
        data_types_detected=[data_type],
    )
