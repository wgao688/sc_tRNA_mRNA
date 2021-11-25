# This script runs through how to extract a scRNA-seq count matrix using Seurat and SCTransform.
# Replace all "foo" file names with actual file names. 
# The file paths that need to be supplied are indicated at the beginning of the script.

# Load libraries
library(Seurat)

#### IMPORTANT: supply the following file paths

# Insert path to count matrix
count_matrix_file = "foo.txt"
# Insert path to file with names of protein-coding genes, which should be named in the same way as the count matrix. 
# These will be the genes that are kept in the scRNA-seq count matrix for further analysis.
coding_gene_names <- read.table("foo_protein.txt")
# Insert path to output file
output <- "foo.out"

####

# load count matrix
scRNA_seq_matrix <- readRDS(count_matrix)

# trim scRNA-seq data only to protein-coding genes, which are provided by the coding_gene_names input file
scRNA_seq_matrix <- scRNA_seq_matrix[coding_gene_names, ]

# create Seurat object and perform SCTransform
scRNA_seq_transform <- CreateSeuratObject(counts = scRNA_seq_matrix)
scRNA_seq_transform <- SCTransform(object = scRNA_seq_transform, assay = "RNA", ncells = 5000, conserve.memory = TRUE)

# save SCTransform corrected count matrix as RDS
corrected_count_matrix <- GetAssayData(object = scRNA_seq_transform, assay = "SCT", slot = "count")
saveRDS(corrected_count_matrix, output)