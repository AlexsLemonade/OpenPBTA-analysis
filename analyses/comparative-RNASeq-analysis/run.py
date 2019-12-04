"""Wrapper that imports and runs all comparative rnaseq analyses scripts"""
import os
import sys
import correlation_matrix as cm

# This script should always run as if it were being called from
# the directory it lives in.
os.chdir(sys.path[0])

# (input filename, desired prefix)
datasets_to_process = [
    ("pbta-gene-expression-rsem-tpm.polya.rds", "rsem-tpm-polya-"),
    ("pbta-gene-expression-rsem-tpm.stranded.rds", "rsem-tpm-stranded-")
]

for filename, prefix in datasets_to_process:
    # Create correlation matrix
    expression = cm.prepare_expression(filename, prefix=prefix)
    correlation_matrix = cm.calculate_correlation(expression, prefix=prefix)

    # Add future pipelines here!
