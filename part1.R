# Part 1
#Install Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db", update=FALSE)

library(org.Mm.eg.db)
library(ggplot2)

#Load Data
counts <- read.delim("SRP082327.tsv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim("metadata_SRP082327.tsv", header = TRUE, row.names = 1)

#Check Matrix Size
cat("Matrix Size:", dim(counts)[1], "genes x", dim(counts)[2], "samples\n")

#Convert Ensembl IDs to Gene Names
clean_ensembl_ids <- gsub("\\..*", "", rownames(counts))

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = clean_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

gene_symbols <- ifelse(is.na(gene_map), rownames(counts), gene_map)
rownames(counts) <- make.unique(gene_symbols)

#Calculation
log_counts <- log2(counts + 1)
gene_medians <- apply(log_counts, 1, median)

df_density <- data.frame(MedianExpression = gene_medians)

#Density Plot
ggplot(df_density, aes(x = MedianExpression)) +
  geom_density(fill = "#69b3a2", color = "#e9ecef", alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Per-Gene Median Expression Distribution (Log2)",
    x = "Log2(Expression + 1) Median",
    y = "Density"
  )