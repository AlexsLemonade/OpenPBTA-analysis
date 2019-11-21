
"""Wrapper that imports and runs all comparative rnaseq analyses scripts"""

import os
import correlation_matrix

data="../../data"

# Create correlation matrix
correlation_matrix.run(os.path.join(data, "pbta-gene-expression-rsem-tpm.polya.rds"))
