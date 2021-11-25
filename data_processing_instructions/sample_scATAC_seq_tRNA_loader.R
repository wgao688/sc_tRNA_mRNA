# This script runs through how to extract a scATAC-seq count matrix for tRNA genes using Signac. 
# Replace all "foo" file names with actual file names. 
# The file paths that need to be supplied are indicated at the beginning of the script.

# Load libraries
library(Biostrings)
library(dplyr)
library(tibble)
library(tidyr)
library(Signac)
library(Seurat)
library(tRNAscanImport)

#### IMPORTANT: supply the following file paths

# insert path to count matrix RDS file
counts <- readRDS("foo.RDS")
# insert path to fragment file (note: you can list also list a vector of different fragment files rather than just a single one using c() function)
fragment_file = "foo.fragments.txt.gz"
# insert path to tRNA annotation. Use tRNA annotations, downloaded from GtRNAdb: http://gtrnadb.ucsc.edu/index.html. This should be the file that ends in .ss from the zip folder
tRNA_granges <- import.tRNAscanAsGRanges("foo.ss")
# insert path to output file
output_file <- "foo.out"

####

# Create Chromatin Assay in Signac
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c("-", "-"),
  fragments = fragment_file,
  min.cells = 10,
  min.features = 0
)

# Create Seurat object
chrom_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)

# Create tRNA expression matrix, using the barcodes of the cells as the column names
cell_barcodes <- colnames(chrom_obj)

# Extract metadata from tRNA gene annotation object
tRNA_granges <- tRNA_granges[which(tRNA_granges$tRNA_anticodon != "NNN")] # remove tRNAs which have unknown/undefined anticodon
ranges(tRNA_granges) <- ranges(tRNA_granges) + 100 # include the flanking 100 nt of each tRNA
# create cut data matrix for tRNA regions of genome
names(tRNA_granges) <- 1:length(tRNA_granges)
tRNA_granges$gene_name <- paste(names(tRNA_granges), tRNA_granges$tRNA_type, tRNA_granges$tRNA_anticodon, sep='_')
tRNA_granges$gene_biotype <- 'tRNA'

tRNA_expression_matrix <- FeatureMatrix(
  fragments = Fragments(chrom_obj[['peaks']]), 
  features = tRNA_granges, 
  cells = cell_barcodes,
  sep = c("-", "-"),
  verbose = TRUE
)

# save tRNA expression matrix as RDS
saveRDS(tRNA_expression_matrix, file = output_file)