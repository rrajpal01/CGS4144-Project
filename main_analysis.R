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

# Part 2
#Install Packages
if (!requireNamespace("Rtsne", quietly = TRUE)) install.packages("Rtsne")
if (!requireNamespace("umap", quietly = TRUE)) install.packages("umap")

library(Rtsne)
library(umap)
library(ggplot2)

#Divide Data to Two Age Groups
samples_2groups <- metadata$refinebio_age %in% c("2", "12")
log_counts_2g <- log_counts[, samples_2groups]
meta_2g <- metadata[samples_2groups, ]

group_2g <- as.factor(meta_2g$refinebio_age)

#Use All Genes
gene_vars <- apply(log_counts_2g, 1, var)
non_zero <- gene_vars > 0
expr_sub <- t(log_counts_2g[non_zero, ])

#PCA Plot
pca_res <- prcomp(expr_sub, scale. = TRUE)

# Calculate Percentage of Variance
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_var <- round(var_explained[1] * 100, 2)
pc2_var <- round(var_explained[2] * 100, 2)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  Age = group_2g
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Age)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "PCA Plot (2 Years vs 12 Weeks)",
    x = paste0("PC1: ", pc1_var, "% variance"),
    y = paste0("PC2: ", pc2_var, "% variance"),
    color = "Age Group"
  )
print(p_pca)

#t-SNE Plot
set.seed(42)
tsne_res <- Rtsne(expr_sub, perplexity = 30, check_duplicates = FALSE)
tsne_df <- data.frame(
  tSNE1 = tsne_res$Y[, 1],
  tSNE2 = tsne_res$Y[, 2],
  Age = group_2g
)

p_tsne <- ggplot(tsne_df, aes(x = tSNE1, y = tSNE2, color = Age)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "t-SNE Plot (2 Years vs 12 Weeks)",
    x = "t-SNE Dimension 1",
    y = "t-SNE Dimension 2",
    color = "Age Group"
  )
print(p_tsne)

#UMAP Plot
umap_res <- umap(expr_sub)
umap_df <- data.frame(
  UMAP1 = umap_res$layout[, 1],
  UMAP2 = umap_res$layout[, 2],
  Age = group_2g
)

p_umap <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Age)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "UMAP Plot (2 Years vs 12 Weeks)",
    x = "UMAP Dimension 1",
    y = "UMAP Dimension 2",
    color = "Age Group"
  )
print(p_umap)