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

#Part 3
library(limma)
library(ggplot2)
library(org.Mm.eg.db)

#Two Age Groups
target_samples <- metadata$refinebio_age %in% c("2", "12")
log_counts_sub <- log_counts[, target_samples]
meta_sub <- metadata[target_samples, ]
group <- factor(meta_sub$refinebio_age, levels = c("2", "12"))
design <- model.matrix(~ group)

#Fit Linear Model with Limma
fit <- lmFit(log_counts_sub, design)
fit <- eBayes(fit)

#Differential Expression
res_df <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")

#Store Feature Names
res_df$gene_symbol <- rownames(res_df)

#Gene Name Column
res_df$gene_name <- NA

#Map
symbols_in_df <- res_df$gene_symbol[!grepl("^ENSMUSG", res_df$gene_symbol)]
valid_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
matching_symbols <- intersect(symbols_in_df, valid_symbols)

if (length(matching_symbols) > 0) {
  name_map_from_symbol <- mapIds(
    org.Mm.eg.db,
    keys = matching_symbols,
    column = "GENENAME",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  res_df$gene_name <- name_map_from_symbol[res_df$gene_symbol]
}

ensembl_in_df <- sub("\\..*", "", res_df$gene_symbol[grepl("^ENSMUSG", res_df$gene_symbol)])
valid_ensembl <- keys(org.Mm.eg.db, keytype = "ENSEMBL")
matching_ensembl <- intersect(ensembl_in_df, valid_ensembl)

if (length(matching_ensembl) > 0) {
  name_map_from_ensembl <- mapIds(
    org.Mm.eg.db,
    keys = matching_ensembl,
    column = "GENENAME",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  symbol_map_from_ensembl <- mapIds(
    org.Mm.eg.db,
    keys = matching_ensembl,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  unmapped_mask <- grepl("^ENSMUSG", res_df$gene_symbol)
  clean_unmapped <- sub("\\..*", "", res_df$gene_symbol[unmapped_mask])
  
  updated_names <- name_map_from_ensembl[clean_unmapped]
  updated_symbols <- symbol_map_from_ensembl[clean_unmapped]
  
  res_df$gene_name[unmapped_mask] <- ifelse(is.na(updated_names), res_df$gene_name[unmapped_mask], updated_names)
  res_df$gene_symbol[unmapped_mask] <- ifelse(is.na(updated_symbols), res_df$gene_symbol[unmapped_mask], updated_symbols)
}

#N/A Gene Names
res_df$gene_name[is.na(res_df$gene_name)] <- "Unknown / Unannotated"

# Rename columns to standard DEG format
colnames(res_df)[colnames(res_df) == "logFC"] <- "log2FoldChange"
colnames(res_df)[colnames(res_df) == "P.Value"] <- "pvalue"
colnames(res_df)[colnames(res_df) == "adj.P.Val"] <- "padj"

#Sort
res_df <- res_df[order(res_df$padj), ]

#Clean
res_df <- res_df[, c("gene_symbol", "gene_name", "log2FoldChange", "AveExpr", "t", "pvalue", "padj", "B")]

#Save CSV
if (!dir.exists("results")) dir.create("results")

write.csv(res_df, file = "results/full_differential_expression_results.csv", row.names = FALSE)

top50_degs <- head(res_df, 50)
write.csv(top50_degs, file = "results/top50_differentially_expressed_genes.csv", row.names = FALSE)

#Volcano Plot
res_df$diffexpressed <- "NO"
res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "UP"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "DOWN"

volcano_p <- ggplot(data = res_df, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "#2b5c8f", "NO" = "grey", "UP" = "#e76f51")) +
  geom_vline(xintercept = c(-1, 1), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "black", linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot: 2 Years vs 12 Weeks",
    subtitle = "Differential Expression Analysis (limma)",
    x = "Log2 Fold Change",
    y = "-Log10 P-Value",
    color = "Expression Status"
  )
print(volcano_p)

#Part 4
library(pheatmap)

#Filter DEGs
sig_genes_df <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

#Sort by p-value
sig_genes_df <- sig_genes_df[order(sig_genes_df$padj), ]

#Take Top 50 DEGs
top_deg_count <- min(50, nrow(sig_genes_df))
sig_genes_top <- sig_genes_df[1:top_deg_count, ]

#Extract Expression Matrix
heatmap_matrix <- log_counts_sub[rownames(log_counts_sub) %in% rownames(sig_genes_top), ]

#Gene Symbols as Row Names
if ("gene_symbol" %in% colnames(sig_genes_top)) {
  matched_symbols <- sig_genes_top$gene_symbol[match(rownames(heatmap_matrix), rownames(sig_genes_top))]
  rownames(heatmap_matrix) <- make.unique(as.character(matched_symbols))
}

#Sample Annotation Sidebar
annotation_col <- data.frame(
  Age_Group = factor(meta_sub$refinebio_age, levels = c("2", "12"))
)
rownames(annotation_col) <- colnames(heatmap_matrix)

#Color
ann_colors <- list(
  Age_Group = c("2" = "#2b5c8f", "12" = "#e76f51")
)

#Heatmap
if (!dir.exists("results")) dir.create("results")

#Save PDF
pdf("results/significant_genes_heatmap.pdf", width = 8, height = 10)
pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = FALSE,             
  show_rownames = TRUE,              
  fontsize_row = 8,                 
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = paste("Heatmap of Top", top_deg_count, "Significant DEGs (2 Years vs 12 Weeks)")
)
dev.off()

#Save PNG
png("results/significant_genes_heatmap.png", width = 900, height = 1100, res = 130)
pheatmap(
  heatmap_matrix,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = paste("Heatmap of Top", top_deg_count, "Significant DEGs (2 Years vs 12 Weeks)")
)
dev.off()

# Part 5
library(topGO)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)

if (!dir.exists("results")) dir.create("results")

#topGO & Gene Ontology

p_vals <- res_df$padj
names(p_vals) <- gsub("\\..*", "", res_df$gene_symbol)

valid_idx <- !is.na(p_vals) & names(p_vals) != ""
p_vals <- p_vals[valid_idx]

entrez_map <- mapIds(org.Mm.eg.db, keys = names(p_vals), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

valid_entrez <- !is.na(entrez_map)
gene_universe <- p_vals[valid_entrez]
names(gene_universe) <- entrez_map[valid_entrez]

#Deduplicate
gene_universe <- tapply(gene_universe, names(gene_universe), min)

top_diff_genes <- function(all_score) {
  return(all_score < 0.05)
}

topgo_data <- new(
  "topGOdata",
  description = "GO Enrichment using topGO",
  ontology    = "BP",
  allGenes    = gene_universe,
  geneSel     = top_diff_genes,
  nodeSize    = 10,
  annot       = annFUN.org,
  mapping     = "org.Mm.eg.db",
  ID          = "entrez"
)

result_classic <- runTest(topgo_data, algorithm = "classic", statistic = "fisher")

top_terms_count <- min(200, length(score(result_classic)))
topgo_results <- GenTable(
  topgo_data,
  classicFisher = result_classic,
  orderBy = "classicFisher",
  topNodes = top_terms_count
)

topgo_clean <- topgo_results %>%
  transmute(
    Method = "topGO",
    Ontology = "Gene Ontology (BP)",
    Term_ID = GO.ID,
    Term_Description = Term,
    Annotated = Annotated,
    Significant = Significant,
    Expected = Expected,
    p_value = as.numeric(classicFisher)
  ) %>%
  arrange(p_value)

write.csv(topgo_clean, "results/enrichment_topGO_GO_BP.csv", row.names = FALSE)

#clusterProfiler & Gene Ontology
sig_degs <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

deg_genes <- sig_degs$gene_symbol
universe_genes <- res_df$gene_symbol

#Convert Gene Symbols
deg_mapped <- bitr(deg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
if (nrow(deg_mapped) < 5) {
  deg_mapped <- bitr(gsub("\\..*", "", deg_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}

universe_mapped <- bitr(universe_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
if (nrow(universe_mapped) < 10) {
  universe_mapped <- bitr(gsub("\\..*", "", universe_genes), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
}

deg_entrez <- deg_mapped$ENTREZID
universe_entrez <- universe_mapped$ENTREZID

#GO Biological Process Enrichment
ego_bp <- enrichGO(
  gene          = deg_entrez,
  universe      = universe_entrez,
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

cp_res <- as.data.frame(ego_bp)

if (nrow(cp_res) > 0) {
  cp_clean <- cp_res %>%
    mutate(
      k = as.numeric(sub("/.*", "", GeneRatio)),
      n = as.numeric(sub(".*/", "", GeneRatio)),
      M = as.numeric(sub("/.*", "", BgRatio)),
      N = as.numeric(sub(".*/", "", BgRatio)),
      Exp = round(n * (M / N), 2)
    ) %>%
    transmute(
      Method           = "clusterProfiler",
      Ontology         = "Gene Ontology (BP)",
      Term_ID          = ID,
      Term_Description = Description,
      Annotated        = M,
      Significant      = k,
      Expected         = Exp,
      p_value          = pvalue,
      p.adjust         = p.adjust,
      qvalue           = qvalue,
      geneID           = geneID,
      Count            = Count
    ) %>%
    arrange(p_value)
  
  #Save CSV
  write.csv(cp_clean, "results/enrichment_clusterProfiler_GO_BP.csv", row.names = FALSE)
} else {
  message("No significantly enriched terms found at the specified cutoff.")
}

if (nrow(cp_res) > 0) {
  p_dot <- dotplot(ego_bp, showCategory = 15) + 
    ggtitle("clusterProfiler: Top 15 GO Biological Processes")
  print(p_dot)
  
  p_bar <- barplot(ego_bp, showCategory = 15) + 
    ggtitle("clusterProfiler: Top 15 GO Biological Processes")
  print(p_bar)
  
  ggsave("results/clusterProfiler_dotplot.png", plot = p_dot, width = 8, height = 6)
  ggsave("results/clusterProfiler_barplot.png", plot = p_bar, width = 8, height = 6)
}

#gProfiler2 & Gene Ontology

if (!requireNamespace("gprofiler2", quietly = TRUE)) install.packages("gprofiler2")

library(gprofiler2)
library(dplyr)
library(org.Mm.eg.db)

if (!dir.exists("results")) dir.create("results")

sig_degs <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

deg_genes <- sig_degs$gene_symbol
universe_genes <- res_df$gene_symbol

#Run gProfiler2
gost_res <- gost(
  query = deg_genes,
  organism = "mmusculus",
  ordered_query = FALSE,
  multi_query = FALSE,
  significant = TRUE,
  exclude_iea = FALSE,
  measure_underrepresentation = FALSE,
  evcodes = TRUE,
  user_threshold = 0.05,
  correction_method = "g_SCS",
  domain_scope = "custom",
  custom_bg = universe_genes,
  sources = c("GO:BP")
)

if (!is.null(gost_res$result) && nrow(gost_res$result) > 0) {
  gp_df <- gost_res$result
  
  gp_clean <- gp_df %>%
    transmute(
      Method           = "gProfiler2",
      Ontology         = "Gene Ontology (BP)",
      Term_ID          = term_id,
      Term_Description = term_name,
      Annotated        = term_size,
      Significant      = intersection_size,
      Expected         = round(query_size * (term_size / effective_domain_size), 2),
      p_value          = p_value
    ) %>%
    arrange(p_value)
  
  write.csv(gp_clean, "results/enrichment_gProfiler2_GO_BP.csv", row.names = FALSE)
} else {
  message("No significantly enriched terms found using gProfiler2.")
}

#Visualization
if (!is.null(gost_res$result)) {
  p_gost <- gostplot(gost_res, capped = TRUE, interactive = FALSE)
  ggsave("results/gprofiler2_gostplot.png", plot = p_gost, width = 9, height = 6)
}

# Part 6
# Read Generated Enrichment Results
files_to_read <- list.files("results", pattern = "^enrichment_.*\\.csv$", full.names = TRUE)

all_methods_list <- list()

for (f in files_to_read) {
  df_in <- read.csv(f, stringsAsFactors = FALSE)
  if (nrow(df_in) > 0) {
    all_methods_list[[df_in$Method[1]]] <- df_in
  }
}

#Total Number of Methods Tested
total_methods_count <- length(all_methods_list)

#Unique Term IDs
all_terms_df <- do.call(rbind, lapply(all_methods_list, function(x) {
  x[, c("Term_ID", "Term_Description")]
})) %>% distinct(Term_ID, .keep_all = TRUE)

combined_wide <- all_terms_df

for (method_name in names(all_methods_list)) {
  m_df <- all_methods_list[[method_name]]
  m_sub <- m_df %>%
    select(Term_ID, p_value, Significant) %>%
    rename(
      !!paste0(method_name, "_pvalue") := p_value,
      !!paste0(method_name, "_sig_count") := Significant
    )
  
  combined_wide <- left_join(combined_wide, m_sub, by = "Term_ID")
}

#Summary Columns
p_col_names <- colnames(combined_wide)[grepl("_pvalue$", colnames(combined_wide))]

combined_wide <- combined_wide %>%
  rowwise() %>%
  mutate(
    methods_included = sum(!is.na(c_across(all_of(p_col_names)))),
    methods_significant = sum(c_across(all_of(p_col_names)) < 0.05, na.rm = TRUE),
    mean_p_value = mean(c_across(all_of(p_col_names)), na.rm = TRUE)
  ) %>%
  ungroup() %>%
  filter(methods_included > 0) %>% 
  arrange(desc(methods_significant), desc(methods_included), mean_p_value)

#Save CSV
write.csv(combined_wide, "results/joint_enrichment_results_wide.csv", row.names = FALSE)

# Part 7
library(dplyr)

#Table from Part 6
joint_df <- read.csv("results/joint_enrichment_results_wide.csv", stringsAsFactors = FALSE)

#Extract Top 10 Terms
top10_terms <- joint_df %>%
  arrange(desc(methods_significant), desc(methods_included), mean_p_value) %>%
  head(10)

# Save Top 10 Table
write.csv(top10_terms, "results/top10_enrichment_results_combined.csv", row.names = FALSE)
