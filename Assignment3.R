# Setup
required_packages <- c("ggplot2", "pheatmap", "cluster", "stats", "dplyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}
if (!requireNamespace("kernlab", quietly = TRUE)) install.packages("kernlab")

library(org.Mm.eg.db)
library(ggplot2)
library(pheatmap)
library(cluster)
library(dplyr)
library(kernlab)
library(RColorBrewer)

if (!dir.exists("results")) dir.create("results")

# Part 1
# Load Data
counts <- read.delim("SRP082327.tsv", header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim("metadata_SRP082327.tsv", header = TRUE, row.names = 1)

# Clean Ensembl IDs
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

log_counts <- log2(counts + 1)

# Subset metadata
samples_2groups <- metadata$refinebio_age %in% c("2", "12")
log_counts_sub <- log_counts[, samples_2groups]
meta_sub <- metadata[samples_2groups, ]

# Create age group factor
age_group <- factor(meta_sub$refinebio_age, levels = c("2", "12"))

expr_matrix <- log_counts_sub
gene_vars <- apply(expr_matrix, 1, var)
sorted_genes <- names(sort(gene_vars, decreasing = TRUE))
expr_matrix <- expr_matrix[sorted_genes, ]

cat("Step 1 Complete: Data loaded and preprocessed for", ncol(log_counts_sub), "samples.\n")

# Extract Top 5000 Most Variable Genes
top5000_genes <- sorted_genes[1:min(5000, length(sorted_genes))]
top5000_matrix <- expr_matrix[top5000_genes, ]
expr_5k <- t(top5000_matrix)

k_values <- 2:6

# Compute PCA
pca_res <- prcomp(expr_5k, scale. = TRUE)
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)

# Part 2
# K-Means
set.seed(42)
kmeans_results <- lapply(k_values, function(k) {
  kmeans(expr_5k, centers = k, nstart = 25)
})
names(kmeans_results) <- paste0("k", k_values)

k2_clusters <- factor(kmeans_results$k2$cluster)
k3_clusters <- factor(kmeans_results$k3$cluster)
k4_clusters <- factor(kmeans_results$k4$cluster)
k5_clusters <- factor(kmeans_results$k5$cluster)
k6_clusters <- factor(kmeans_results$k6$cluster)

# Elbow Plot
tot_withinss <- sapply(kmeans_results, function(x) x$tot.withinss)
elbow_df <- data.frame(k = k_values, Tot_WithinSS = tot_withinss)

p_elbow <- ggplot(elbow_df, aes(x = k, y = Tot_WithinSS)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "red", size = 3) +
  theme_minimal() +
  labs(
    title = "K-Means Elbow Method (Top 5000 Variable Genes)",
    x = "Number of Clusters (k)",
    y = "Total Within-Cluster Sum of Squares"
  )
print(p_elbow)
ggsave("results/kmeans_elbow_plot.png", plot = p_elbow, width = 7, height = 5)

# PCA Plot
pca_df_kmeans <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  KMeans_Cluster = k2_clusters,
  True_Age_Group = age_group
)

p_kmeans_pca <- ggplot(pca_df_kmeans, aes(x = PC1, y = PC2, color = KMeans_Cluster, shape = True_Age_Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  theme_minimal() +
  labs(
    title = "K-Means Clustering Analysis (Top 5000 Variable Genes, k=2)",
    x = paste0("PC1: ", round(var_explained[1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(var_explained[2] * 100, 2), "% variance"),
    color = "K-Means Cluster",
    shape = "True Age Group"
  )
print(p_kmeans_pca)
ggsave("results/kmeans_pca_5000genes.png", plot = p_kmeans_pca, width = 8, height = 6)

# PAM Clustering
set.seed(42)
pam_results <- lapply(k_values, function(k) {
  pam(expr_5k, k = k)
})
names(pam_results) <- paste0("k", k_values)

pam_k2_clusters <- factor(pam_results$k2$clustering)
pam_k3_clusters <- factor(pam_results$k3$clustering)
pam_k4_clusters <- factor(pam_results$k4$clustering)
pam_k5_clusters <- factor(pam_results$k5$clustering)
pam_k6_clusters <- factor(pam_results$k6$clustering)

# Average Silhouette Width Plot
avg_sil_widths <- sapply(pam_results, function(x) x$silinfo$avg.width)
pam_sil_df <- data.frame(k = k_values, Avg_Silhouette = avg_sil_widths)

p_pam_sil <- ggplot(pam_sil_df, aes(x = k, y = Avg_Silhouette)) +
  geom_line(color = "forestgreen", linewidth = 1) +
  geom_point(color = "darkgreen", size = 3) +
  theme_minimal() +
  labs(
    title = "PAM Average Silhouette Width (Top 5000 Variable Genes)",
    x = "Number of Clusters (k)",
    y = "Average Silhouette Width"
  )
print(p_pam_sil)
ggsave("results/pam_silhouette_plot.png", plot = p_pam_sil, width = 7, height = 5)

# PCA Plot
pca_df_pam <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PAM_Cluster = pam_k2_clusters,
  True_Age_Group = age_group
)

p_pam_pca <- ggplot(pca_df_pam, aes(x = PC1, y = PC2, color = PAM_Cluster, shape = True_Age_Group)) +
  geom_point(size = 3.5, alpha = 0.85) +
  theme_minimal() +
  labs(
    title = "PAM Clustering Analysis (Top 5000 Variable Genes, k=2)",
    x = paste0("PC1: ", round(var_explained[1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(var_explained[2] * 100, 2), "% variance"),
    color = "PAM Cluster",
    shape = "True Age Group"
  )
print(p_pam_pca)
ggsave("results/pam_pca_5000genes.png", plot = p_pam_pca, width = 8, height = 6)

# Spectral Clustering
spectral_results <- list()
for (k in k_values) {
  set.seed(42)
  spectral_results[[paste0("k", k)]] <- specc(expr_5k, centers = k, kernel = "rbfdot")
}

spec_k2_clusters <- factor(as.vector(spectral_results$k2))
spec_k3_clusters <- factor(as.vector(spectral_results$k3))
spec_k4_clusters <- factor(as.vector(spectral_results$k4))
spec_k5_clusters <- factor(as.vector(spectral_results$k5))
spec_k6_clusters <- factor(as.vector(spectral_results$k6))

# Spectral PCA Facet Plot
spec_k_df_list <- lapply(k_values, function(k) {
  data.frame(
    PC1 = pca_res$x[, 1],
    PC2 = pca_res$x[, 2],
    Cluster = factor(as.vector(spectral_results[[paste0("k", k)]])),
    k_val = paste0("k = ", k),
    True_Age_Group = age_group
  )
})

spec_k_df <- do.call(rbind, spec_k_df_list)

p_spec_k_compare <- ggplot(spec_k_df, aes(x = PC1, y = PC2, color = Cluster, shape = True_Age_Group)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ k_val, ncol = 3) +
  theme_minimal() +
  labs(
    title = "Spectral Clustering Membership Comparison Across k = 2 to 6",
    x = paste0("PC1: ", round(var_explained[1] * 100, 2), "% variance"),
    y = paste0("PC2: ", round(var_explained[2] * 100, 2), "% variance"),
    color = "Cluster ID",
    shape = "True Age Group"
  )

print(p_spec_k_compare)
ggsave("results/spectral_k_comparison_pca.png", plot = p_spec_k_compare, width = 12, height = 8, dpi = 300)


# Cluster Membership Transition
membership_matrix <- data.frame(
  Sample = rownames(expr_5k),
  `k = 2` = spec_k2_clusters,
  `k = 3` = spec_k3_clusters,
  `k = 4` = spec_k4_clusters,
  `k = 5` = spec_k5_clusters,
  `k = 6` = spec_k6_clusters,
  True_Age = age_group,
  check.names = FALSE
)

library(reshape2)
membership_long <- melt(membership_matrix, id.vars = c("Sample", "True_Age"), 
                        variable.name = "k_level", value.name = "Cluster")

p_spec_bars <- ggplot(membership_long, aes(x = k_level, fill = Cluster)) +
  geom_bar(position = "fill", color = "black", linewidth = 0.2) +
  facet_wrap(~ True_Age) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal() +
  labs(
    title = "Spectral Cluster Membership Proportion across k (Split by Age Group)",
    x = "Cluster Resolution (k)",
    y = "Percentage of Samples",
    fill = "Cluster ID"
  )

print(p_spec_bars)
ggsave("results/spectral_membership_bar_comparison.png", plot = p_spec_bars, width = 10, height = 6, dpi = 300)

# Cluster Membership across k = 2..6
get_k_breakdown <- function(cluster_factor, age_factor) {
  clusters <- levels(cluster_factor)
  parts <- sapply(clusters, function(cl) {
    sub_ages <- age_factor[cluster_factor == cl]
    c2  <- sum(sub_ages == "2")
    c12 <- sum(sub_ages == "12")
    tot <- length(sub_ages)
    sprintf("%d (%d 2Y / %d 12W)", tot, c2, c12)
  })
  paste(parts, collapse = ", ")
}

k_table_rows <- list()
for (k in 2:6) {
  k_key <- paste0("k", k)
  
  km_cl   <- kmeans_results[[k_key]]$cluster
  pam_cl  <- pam_results[[k_key]]$clustering
  spec_cl <- as.vector(spectral_results[[k_key]])
  
  k_table_rows[[length(k_table_rows) + 1]] <- data.frame(
    `k` = paste0("k = ", k),
    `K-Means Clusters` = get_k_breakdown(factor(km_cl), age_group),
    `PAM Clusters`     = get_k_breakdown(factor(pam_cl), age_group),
    `Spectral Clusters`= get_k_breakdown(factor(spec_cl), age_group),
    check.names = FALSE
  )
}

table1_k_membership <- do.call(rbind, k_table_rows)
print(table1_k_membership)
write.csv(table1_k_membership, "results/k_membership_comparison.csv", row.names = FALSE)


# Different Number of Genes
run_chisq <- function(cluster_vec, age_vec) {
  contingency_tbl <- table(cluster_vec, age_vec)
  
  res <- suppressWarnings(chisq.test(contingency_tbl, correct = FALSE))
  stat <- as.numeric(res$statistic)
  
  df_val <- (nrow(contingency_tbl) - 1) * (ncol(contingency_tbl) - 1)
  raw_p <- pchisq(stat, df = df_val, lower.tail = FALSE)
  
  return(sprintf("stat = %.4f, p = %.2e", stat, raw_p))
}

get_two_cluster_str <- function(cluster_vec, age_vec) {
  cl_fac <- factor(cluster_vec)
  levels(cl_fac) <- c("1", "2")
  
  sub1 <- age_vec[cl_fac == "1"]
  c1_tot <- length(sub1)
  c1_a2  <- sum(sub1 == "2")
  c1_a12 <- sum(sub1 == "12")
  
  sub2 <- age_vec[cl_fac == "2"]
  c2_tot <- length(sub2)
  c2_a2  <- sum(sub2 == "2")
  c2_a12 <- sum(sub2 == "12")
  
  sprintf("C1: %d (%d 2Y / %d 12W) | C2: %d (%d 2Y / %d 12W)", 
          c1_tot, c1_a2, c1_a12, c2_tot, c2_a2, c2_a12)
}

gene_subsets <- c(10, 100, 1000, 10000)

kmeans_gene_results <- list()
pam_gene_results    <- list()
spec_gene_results   <- list()

gene_summary_rows <- list()

for (n in gene_subsets) {
  n_str <- as.character(n)
  
  sub_expr <- expr_matrix[1:min(n, nrow(expr_matrix)), ]
  sub_data <- t(sub_expr)
  
  set.seed(42)
  km_res <- kmeans(sub_data, centers = 2, nstart = 25)$cluster
  
  set.seed(42)
  pam_res <- cluster::pam(sub_data, k = 2)$clustering
  
  set.seed(42)
  spec_res <- kernlab::specc(sub_data, centers = 2)
  spec_cl <- as.vector(spec_res)
  
  kmeans_gene_results[[n_str]] <- km_res
  pam_gene_results[[n_str]]    <- pam_res
  spec_gene_results[[n_str]]   <- spec_cl
  
  km_str  <- get_two_cluster_str(km_res, age_group)
  pam_str <- get_two_cluster_str(pam_res, age_group)
  sp_str  <- get_two_cluster_str(spec_cl, age_group)
  
  km_chisq  <- run_chisq(km_res, age_group)
  pam_chisq <- run_chisq(pam_res, age_group)
  sp_chisq  <- run_chisq(spec_cl, age_group)
  
  gene_summary_rows[[length(gene_summary_rows) + 1]] <- data.frame(
    `Gene Subset`              = n_str,
    `K-Means Clusters`        = km_str,
    `K-Means Chi-Sq`           = km_chisq,
    `PAM Clusters`            = pam_str,
    `PAM Chi-Sq`               = pam_chisq,
    `Spectral Clusters`       = sp_str,
    `Spectral Chi-Sq`          = sp_chisq,
    check.names = FALSE
  )
}

gene_summary_table <- do.call(rbind, gene_summary_rows)

# Part 3
# Heatmap
graphics.off()

if (!dir.exists("results")) {
  dir.create("results")
}

# Format Annotations
annotation_col <- data.frame(
  True_Age_Group   = as.character(age_group),
  KMeans_Cluster   = as.character(k2_clusters),
  PAM_Cluster      = as.character(pam_k2_clusters),
  Spectral_Cluster = as.character(spec_k2_clusters),
  stringsAsFactors = FALSE
)
rownames(annotation_col) <- colnames(top5000_matrix)

ann_colors <- list(
  True_Age_Group   = c("2" = "#2b5c8f", "12" = "#e76f51"),
  KMeans_Cluster   = c("1" = "#E41A1C", "2" = "#377EB8"),
  PAM_Cluster      = c("1" = "#4DAF4A", "2" = "#984EA3"),
  Spectral_Cluster = c("1" = "#FF7F00", "2" = "#FFFF33")
)

png("results/heatmap_5000genes.png", width = 12, height = 10, units = "in", res = 300)

# Generate Heatmap
pheatmap(
  top5000_matrix,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,  # Set to FALSE because 811 sample names will crash graphics
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  legend = TRUE,
  annotation_legend = TRUE,
  silent = TRUE,         # Prevents rendering lag in RStudio plot pane
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Part 3: Gene Expression Heatmap (Top 5,000 Genes) with Cluster Annotations"
)

dev.off()

# Part 4
# Chi-squared
run_chisq_test <- function(vec1, vec2, test_label, method_name, k_val, gene_count) {
  contingency_tbl <- table(vec1, vec2)
  
  df_val <- (nrow(contingency_tbl) - 1) * (ncol(contingency_tbl) - 1)

  res <- suppressWarnings(chisq.test(contingency_tbl, correct = FALSE))
  stat <- as.numeric(res$statistic)
  
  # Exact upper-tail Chi-squared p-value
  raw_p <- pchisq(stat, df = df_val, lower.tail = FALSE)
  
  return(data.frame(
    Method = method_name,
    Number_of_Genes = as.character(gene_count),
    k_Clusters = as.character(k_val),
    Comparison_Label = test_label,
    Chi2_Statistic = round(stat, 4),
    P_Value = raw_p,
    stringsAsFactors = FALSE
  ))
}

# True Age Group vs Clusters
all_stats_df <- list()
# Groups vs. Clusters across k = 2..6 (5000 Genes)
kmeans_all <- list(k2 = k2_clusters, k3 = k3_clusters, k4 = k4_clusters, k5 = k5_clusters, k6 = k6_clusters)
pam_all    <- list(k2 = pam_k2_clusters, k3 = pam_k3_clusters, k4 = pam_k4_clusters, k5 = pam_k5_clusters, k6 = pam_k6_clusters)
spec_all   <- list(k2 = spec_k2_clusters, k3 = spec_k3_clusters, k4 = spec_k4_clusters, k5 = spec_k5_clusters, k6 = spec_k6_clusters)

for (k_name in names(kmeans_all)) {
  k_val <- as.numeric(gsub("k", "", k_name))
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(kmeans_all[[k_name]], age_group, paste0("K-Means (5000 genes, k=", k_val, ") vs Assignment 1 Age Group"), "K-Means", k_val, 5000)
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(pam_all[[k_name]], age_group, paste0("PAM (5000 genes, k=", k_val, ") vs Assignment 1 Age Group"), "PAM", k_val, 5000)
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(spec_all[[k_name]], age_group, paste0("Spectral (5000 genes, k=", k_val, ") vs Assignment 1 Age Group"), "Spectral", k_val, 5000)
}

# Groups vs. Gene Subset Clusters (k = 2)
for (n_genes in gene_subsets) {
  n_str <- as.character(n_genes)
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(kmeans_gene_results[[n_str]], age_group, paste0("K-Means (", n_genes, " genes, k=2) vs Assignment 1 Age Group"), "K-Means", 2, n_genes)
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(pam_gene_results[[n_str]], age_group, paste0("PAM (", n_genes, " genes, k=2) vs Assignment 1 Age Group"), "PAM", 2, n_genes)
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(spec_gene_results[[n_str]], age_group, paste0("Spectral (", n_genes, " genes, k=2) vs Assignment 1 Age Group"), "Spectral", 2, n_genes)
}

# Different Clustering Methods
# Pairwise Comparisons Across k (k=2 vs k=3, k=3 vs k=4, etc.)
all_k_lists <- list("K-Means" = kmeans_all, "PAM" = pam_all, "Spectral" = spec_all)
for (m_name in names(all_k_lists)) {
  m_k_list <- all_k_lists[[m_name]]
  for (i in 2:5) {
    k_curr <- paste0("k", i)
    k_next <- paste0("k", i + 1)
    label_str <- paste0(m_name, " (5000 genes) Resolution Comparison: ", k_curr, " vs ", k_next)
    all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(m_k_list[[k_curr]], m_k_list[[k_next]], label_str, m_name, paste0(i, "_vs_", i + 1), 5000)
  }
}

# Pairwise Feature Subset Comparisons (10 vs 100, 100 vs 1000, 1000 vs 10000 genes)
methods_list <- list("K-Means" = kmeans_gene_results, "PAM" = pam_gene_results, "Spectral" = spec_gene_results)
for (m_name in names(methods_list)) {
  m_res <- methods_list[[m_name]]
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(m_res[["10"]], m_res[["100"]], paste(m_name, "10 vs 100 genes"), m_name, 2, "10_vs_100")
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(m_res[["100"]], m_res[["1000"]], paste(m_name, "100 vs 1000 genes"), m_name, 2, "100_vs_1000")
  all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(m_res[["1000"]], m_res[["10000"]], paste(m_name, "1000 vs 10000 genes"), m_name, 2, "1000_vs_10000")
}

# Cross-Method Agreement Comparisons (K-Means vs PAM vs Spectral at k=2, 5000 genes)
all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(k2_clusters, pam_k2_clusters,  "Inter-Algorithm Comparison: K-Means vs PAM (5000 genes, k=2)",     "K-Means_vs_PAM",     2, 5000)
all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(k2_clusters, spec_k2_clusters, "Inter-Algorithm Comparison: K-Means vs Spectral (5000 genes, k=2)", "K-Means_vs_Spectral", 2, 5000)
all_stats_df[[length(all_stats_df) + 1]] <- run_chisq_test(pam_k2_clusters, spec_k2_clusters, "Inter-Algorithm Comparison: PAM vs Spectral (5000 genes, k=2)",     "PAM_vs_Spectral",     2, 5000)

# Combination and Output
master_stats_df <- do.call(rbind, all_stats_df);

master_stats_df$Adjusted_P_Value <- p.adjust(master_stats_df$P_Value, method = "BH");
master_stats_df$Significant_FDR_0.05 <- ifelse(master_stats_df$Adjusted_P_Value < 0.05, "Yes", "No");

# Append Adjusted P-values for all three algorithms
for (m_name in c("K-Means", "PAM", "Spectral")) {
  m_indices <- match(
    paste(m_name, gene_subsets, 2),
    paste(master_stats_df$Method, master_stats_df$Number_of_Genes, master_stats_df$k_Clusters)
  )
  
  col_name <- paste0(m_name, " Chi2 Adj P")
  gene_summary_table[[col_name]] <- ifelse(
    is.na(m_indices),
    "N/A",
    sprintf("%.2E", master_stats_df$Adjusted_P_Value[m_indices])
  )
}

print(gene_summary_table)
write.csv(gene_summary_table, "results/gene_subsets_cluster_breakdown.csv", row.names = FALSE)
print(master_stats_df, row.names = FALSE);
write.csv(master_stats_df, "results/chisq_statistics.csv", row.names = FALSE);