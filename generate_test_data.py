"""
Generate deterministic test data files for RNA-seq Analysis Platform.
This script creates sample CSV files with reproducible data using seed=42.
"""

import csv
import random

# Set seed for reproducibility
random.seed(42)

# Generate sample_counts.csv (10 samples x 100 genes, integers)
# Using negative binomial distribution approximation
print("Generating sample_counts.csv...")
with open("tests/data/sample_counts.csv", "w", newline="") as f:
    writer = csv.writer(f)

    # Header
    genes = [f"gene_{i + 1}" for i in range(100)]
    writer.writerow(["sample_id"] + genes)

    # Data rows
    for i in range(10):
        sample_id = f"sample_{i + 1}"
        # Approximate negative binomial with geometric distribution
        counts = [random.randint(50, 150) for _ in range(100)]
        writer.writerow([sample_id] + counts)

# Generate sample_metadata.csv (10 samples, sample_id and condition columns)
print("Generating sample_metadata.csv...")
with open("tests/data/sample_metadata.csv", "w", newline="") as f:
    writer = csv.writer(f)

    # Header
    writer.writerow(["sample_id", "condition"])

    # Data rows
    for i in range(10):
        sample_id = f"sample_{i + 1}"
        condition = "control" if i < 5 else "treatment"
        writer.writerow([sample_id, condition])

# Generate normalized_data.csv (10 samples x 100 genes, floats)
print("Generating normalized_data.csv...")
with open("tests/data/normalized_data.csv", "w", newline="") as f:
    writer = csv.writer(f)

    # Header
    genes = [f"gene_{i + 1}" for i in range(100)]
    writer.writerow(["sample_id"] + genes)

    # Data rows
    for i in range(10):
        sample_id = f"sample_{i + 1}"
        # Uniform distribution between 0 and 100
        values = [round(random.uniform(0, 100), 6) for _ in range(100)]
        writer.writerow([sample_id] + values)

print("Test data files generated successfully!")
