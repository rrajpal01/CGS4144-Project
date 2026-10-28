# Load required packages
required_sup_pkgs <- c("e1071", "class", "glmnet", "pROC", "ggplot2", "pheatmap", "dplyr", "org.Mm.eg.db", "kernlab", "cluster")
for (pkg in required_sup_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}

library(e1071)
library(class)
library(glmnet)
library(pROC)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(org.Mm.eg.db)
library(kernlab)
library(cluster)

if (!dir.exists("results")) dir.create("results", showWarnings = FALSE)

# Fixes the connection error
if (file.exists("SRP082327.tsv")) {
  data_path <- "SRP082327.tsv"
  meta_path <- "metadata_SRP082327.tsv"
} else if (file.exists("SRP082327/SRP082327.tsv")) {
  data_path <- "SRP082327/SRP082327.tsv"
  meta_path <- "SRP082327/metadata_SRP082327.tsv"
} else if (file.exists("../SRP082327.tsv")) {
  data_path <- "../SRP082327.tsv"
  meta_path <- "../metadata_SRP082327.tsv"
} else if (file.exists("../SRP082327/SRP082327.tsv")) {
  data_path <- "../SRP082327/SRP082327.tsv"
  meta_path <- "../SRP082327/metadata_SRP082327.tsv"
} else {
  stop("Error: Could not find 'SRP082327.tsv'. Please check your working directory.")
}

# Load Assignment 2 Data
counts <- read.delim(data_path, header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim(meta_path, header = TRUE, row.names = 1)

# Clean Ensembl IDs & Log Transform
clean_ensembl_ids <- gsub("\\..*", "", rownames(counts))
gene_map <- mapIds(org.Mm.eg.db, keys = clean_ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(counts) <- make.unique(ifelse(is.na(gene_map), rownames(counts), gene_map))
log_counts <- log2(counts + 1)

# Subset to 2 groups (2 years vs 12 weeks)
samples_2groups <- metadata$refinebio_age %in% c("2", "12")
log_counts_sub <- log_counts[, samples_2groups]
meta_sub <- metadata[samples_2groups, ]
age_group <- factor(meta_sub$refinebio_age, levels = c("2", "12"))

# Top 5000 Most Variable Genes
gene_vars <- apply(log_counts_sub, 1, var)
sorted_genes <- names(sort(gene_vars, decreasing = TRUE))
top5000_matrix <- log_counts_sub[sorted_genes[1:min(5000, length(sorted_genes))], ]
expr_5k <- t(top5000_matrix)

# Recreate Assignment 3 Clusters
set.seed(42)
k2_clusters <- factor(kmeans(expr_5k, centers = 2, nstart = 25)$cluster)
pam_k2_clusters <- factor(pam(expr_5k, k = 2)$clustering)
spec_k2_clusters <- factor(as.vector(specc(expr_5k, centers = 2, kernel = "rbfdot")))

# Train models to predict Assignment 1 Age Groups
s_dat_a1 <- as.data.frame(expr_5k)
s_dat_a1$Age <- age_group

# SVM for Age Group
set.seed(42)
svm_age_model <- svm(Age ~ ., data = s_dat_a1, probability = TRUE, kernel = "radial")
svm_age_preds <- predict(svm_age_model, s_dat_a1)

# Logistic Regression for Age Group
pca_a1 <- prcomp(expr_5k, scale. = TRUE)
glm_a1_dat <- as.data.frame(pca_a1$x[, 1:50])
glm_a1_dat$Age_bin <- ifelse(age_group == "12", 1, 0)
glm_age_model <- glm(Age_bin ~ ., data = glm_a1_dat, family = binomial)
glm_age_probs <- predict(glm_age_model, type = "response")
glm_age_preds <- ifelse(glm_age_probs > 0.5, "12", "2")

# KNN for Age Group
set.seed(42)
knn_age_preds <- knn.cv(train = expr_5k, cl = age_group, k = 5)

# SVM predicting K-Means Cluster 2
svm_c_dat <- as.data.frame(expr_5k)
svm_c_dat$Cluster <- k2_clusters
svm_clust_model <- svm(Cluster ~ ., data = svm_c_dat, probability = TRUE, kernel = "radial")
svm_clust_preds <- predict(svm_clust_model, svm_c_dat)

# Logistic Regression predicting PAM Cluster 2
glm_c_dat <- as.data.frame(pca_a1$x[, 1:50])
glm_c_dat$Cluster_bin <- ifelse(pam_k2_clusters == "2", 1, 0)
glm_clust_model <- glm(Cluster_bin ~ ., data = glm_c_dat, family = binomial)
glm_clust_probs <- predict(glm_clust_model, type = "response")
glm_clust_preds <- ifelse(glm_clust_probs > 0.5, "2", "1")

# KNN predicting Spectral Cluster 2
knn_clust_preds <- knn.cv(train = expr_5k, cl = spec_k2_clusters, k = 5)

# Build comprehensive sample-by-model matrix
sample_model_matrix <- data.frame(
  Sample = rownames(expr_5k),
  True_Age = as.character(age_group),
  SVM_Age_Pred = as.character(svm_age_preds),
  LogReg_Age_Pred = as.character(glm_age_preds),
  KNN_Age_Pred = as.character(knn_age_preds),
  
  KMeans_Cluster = as.character(k2_clusters),
  SVM_Clust_Pred = as.character(svm_clust_preds),
  PAM_Cluster = as.character(pam_k2_clusters),
  LogReg_Clust_Pred = as.character(glm_clust_preds),
  Spectral_Cluster = as.character(spec_k2_clusters),
  KNN_Clust_Pred = as.character(knn_clust_preds),
  stringsAsFactors = FALSE
)

# Calculate consensus per sample
sample_model_matrix$Models_Predicting_Age_12 <- rowSums(sample_model_matrix[, c("SVM_Age_Pred", "LogReg_Age_Pred", "KNN_Age_Pred")] == "12")
sample_model_matrix$Models_Predicting_Cluster_2 <- rowSums(sample_model_matrix[, c("SVM_Clust_Pred", "LogReg_Clust_Pred", "KNN_Clust_Pred")] == "2")

write.csv(sample_model_matrix, "results/sample_model_predictions_matrix.csv", row.names = FALSE)

stability_cor_test <- cor.test(
  sample_model_matrix$Models_Predicting_Age_12, 
  sample_model_matrix$Models_Predicting_Cluster_2, 
  method = "spearman"
)

cat(sprintf("Stability Correlation Test - rho: %.4f, raw p-value: %.2e\n", 
            stability_cor_test$estimate, stability_cor_test$p.value))

# Save statistical test result
stat_res_df <- data.frame(
  Comparison = "Age Label Stability vs Cluster Prediction Stability",
  Method = "Spearman",
  Rho_Statistic = stability_cor_test$estimate,
  Raw_P_Value = stability_cor_test$p.value,
  Adjusted_P_Value = p.adjust(stability_cor_test$p.value, method = "BH")
)
write.csv(stat_res_df, "results/stability_correlation_statistics.csv", row.names = FALSE)

# Different Gene Counts
gene_subsets <- c(10, 100, 1000, 10000)
perf_results <- list()

for (n in gene_subsets) {
  actual_n <- min(n, nrow(expr_5k))
  sub_expr <- expr_5k[, 1:actual_n]
  
  # SVM
  m_svm <- svm(sub_expr, age_group, probability = TRUE, kernel = "radial")
  pr_svm <- attr(predict(m_svm, sub_expr, probability = TRUE), "probabilities")[, "12"]
  auc_svm <- as.numeric(auc(roc(as.numeric(age_group == "12"), pr_svm, quiet = TRUE)))
  
  # LogReg
  pca_sub <- prcomp(sub_expr, scale. = TRUE)
  lr_dat <- as.data.frame(pca_sub$x[, 1:min(5, ncol(pca_sub$x))])
  lr_dat$Cl <- ifelse(pam_k2_clusters == "2", 1, 0)
  m_glm <- glm(Cl ~ ., data = lr_dat, family = binomial)
  pr_glm <- predict(m_glm, type = "response")
  auc_glm <- as.numeric(auc(roc(as.numeric(pam_k2_clusters == "2"), pr_glm, quiet = TRUE)))
  
  # KNN
  pr_knn_cv <- knn.cv(train = sub_expr, cl = spec_k2_clusters, k = 5)
  auc_knn <- as.numeric(auc(roc(as.numeric(spec_k2_clusters == "2"), as.numeric(pr_knn_cv == "2"), quiet = TRUE)))
  
  perf_results[[length(perf_results) + 1]] <- data.frame(NumGenes = actual_n, Model = "SVM (Age Group)", AUC = auc_svm)
  perf_results[[length(perf_results) + 1]] <- data.frame(NumGenes = actual_n, Model = "Logistic Regression (PAM)", AUC = auc_glm)
  perf_results[[length(perf_results) + 1]] <- data.frame(NumGenes = actual_n, Model = "KNN (Spectral)", AUC = auc_knn)
}

auc_gene_df <- do.call(rbind, perf_results)
write.csv(auc_gene_df, "results/model_auc_by_gene_count_assignment4.csv", row.names = FALSE)

p_auc_gene <- ggplot(auc_gene_df, aes(x = factor(NumGenes), y = AUC, color = Model, group = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  ylim(0, 1.05) +
  theme_minimal() +
  labs(
    title = "Supervised Model Performance (AUC) Across Gene Set Sizes",
    x = "Number of Genes Included",
    y = "Area Under ROC Curve (AUC)",
    color = "Supervised Model"
  )
ggsave("results/supervised_auc_vs_genecount.png", plot = p_auc_gene, width = 9, height = 6)

# Extract signature genes
set.seed(42)
svm_linear_model <- svm(expr_5k, k2_clusters, kernel = "linear", scale = TRUE)
svm_weights <- drop(t(svm_linear_model$coefs) %*% svm_linear_model$SV)
svm_top_genes <- names(sort(abs(svm_weights), decreasing = TRUE))[1:25]

cv_lasso <- cv.glmnet(expr_5k, as.numeric(pam_k2_clusters) - 1, family = "binomial", alpha = 1)
lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
lasso_top_genes <- rownames(expr_5k)[which(lasso_coefs[-1] != 0)]

predictive_signature_genes <- unique(c(svm_top_genes, head(lasso_top_genes, 25)))
if(length(predictive_signature_genes) < 10) predictive_signature_genes <- head(rownames(expr_5k), 30)

heatmap_sig_matrix <- t(expr_5k[, colnames(expr_5k) %in% predictive_signature_genes])

# Map to Gene Symbols
clean_ids <- gsub("\\..*", "", rownames(heatmap_sig_matrix))
gene_map_syms <- mapIds(org.Mm.eg.db, keys = clean_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
symbols_final <- ifelse(is.na(gene_map_syms), rownames(heatmap_sig_matrix), gene_map_syms)
rownames(heatmap_sig_matrix) <- make.unique(as.character(symbols_final))

annotation_col_a4 <- data.frame(
  True_Age = as.character(age_group),
  KMeans_Cl = as.character(k2_clusters),
  PAM_Cl = as.character(pam_k2_clusters),
  Spectral_Cl = as.character(spec_k2_clusters),
  stringsAsFactors = FALSE
)
rownames(annotation_col_a4) <- colnames(heatmap_sig_matrix)

ann_colors_a4 <- list(
  True_Age = c("2" = "#2b5c8f", "12" = "#e76f51"),
  KMeans_Cl = c("1" = "#E41A1C", "2" = "#377EB8"),
  PAM_Cl = c("1" = "#4DAF4A", "2" = "#984EA3"),
  Spectral_Cl = c("1" = "#FF7F00", "2" = "#FFFF33")
)

# Heatmap
png("results/supervised_predictive_genes_heatmap.png", width = 12, height = 10, units = "in", res = 300)
pheatmap(
  heatmap_sig_matrix,
  scale = "row",
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "ward.D2",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 8,
  annotation_col = annotation_col_a4,
  annotation_colors = ann_colors_a4,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
  main = "Heatmap of Predictive Signature Genes with Cluster & Age Annotations"
)
dev.off()