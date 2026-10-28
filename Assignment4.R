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
stopifnot(all(colnames(counts) %in% rownames(metadata)))
metadata <- metadata[colnames(counts), , drop = FALSE]

# Clean Ensembl IDs & Log Transform
clean_ensembl_ids <- gsub("\\..*", "", rownames(counts))
gene_map <- mapIds(org.Mm.eg.db, keys = clean_ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(counts) <- make.unique(ifelse(is.na(gene_map), rownames(counts), gene_map))
log_counts <- log2(counts + 1)

# Subset to 2 age groups (2 months vs 12 months)
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

# Use the same held-out samples for every supervised model
set.seed(42)
train_idx <- unlist(lapply(split(seq_along(age_group), age_group), function(idx) {
  sample(idx, floor(0.8 * length(idx)))
}), use.names = FALSE)
test_idx <- setdiff(seq_along(age_group), train_idx)

# Train all six models on training samples only
x_train <- expr_5k[train_idx, , drop = FALSE]
x_test <- expr_5k[test_idx, , drop = FALSE]

# SVM for Age Group
set.seed(42)
svm_age_model <- svm(x_train, age_group[train_idx], probability = TRUE, kernel = "radial")
svm_age_preds <- predict(svm_age_model, x_test)

# Logistic Regression for Age Group (PCA fitted on training data)
# Retain only non-constant genes for scaled PCA
pca_features <- colnames(x_train)[apply(x_train, 2, sd) > 0]
pca_a1 <- prcomp(x_train[, pca_features, drop = FALSE], scale. = TRUE)
num_pcs <- min(50, ncol(pca_a1$x), nrow(x_train) - 2)
glm_a1_dat <- as.data.frame(pca_a1$x[, seq_len(num_pcs), drop = FALSE])
glm_a1_dat$Age_bin <- as.numeric(age_group[train_idx] == "12")
glm_age_model <- glm(Age_bin ~ ., data = glm_a1_dat, family = binomial)
glm_age_test <- as.data.frame(predict(pca_a1, newdata = x_test[, pca_features, drop = FALSE])[, seq_len(num_pcs), drop = FALSE])
glm_age_probs <- predict(glm_age_model, newdata = glm_age_test, type = "response")
glm_age_preds <- ifelse(glm_age_probs > 0.5, "12", "2")

# KNN for Age Group
set.seed(42)
knn_age_preds <- knn(train = x_train, test = x_test, cl = age_group[train_idx], k = 5)

# SVM predicting K-Means Cluster
set.seed(42)
svm_clust_model <- svm(x_train, k2_clusters[train_idx], probability = TRUE, kernel = "radial")
svm_clust_preds <- predict(svm_clust_model, x_test)

# Logistic Regression predicting PAM Cluster
# Use the same training-only PCA coordinates
glm_c_dat <- as.data.frame(pca_a1$x[, seq_len(num_pcs), drop = FALSE])
glm_c_dat$Cluster_bin <- as.numeric(pam_k2_clusters[train_idx] == "2")
glm_clust_model <- glm(Cluster_bin ~ ., data = glm_c_dat, family = binomial)
glm_clust_probs <- predict(glm_clust_model, newdata = glm_age_test, type = "response")
glm_clust_preds <- ifelse(glm_clust_probs > 0.5, "2", "1")

# KNN predicting Spectral Cluster
set.seed(42)
knn_clust_preds <- knn(train = x_train, test = x_test, cl = spec_k2_clusters[train_idx], k = 5)

# Build sample-by-model matrix for held-out samples only
sample_model_matrix <- data.frame(
  Sample = rownames(expr_5k)[test_idx],
  True_Age = as.character(age_group[test_idx]),
  SVM_Age_Pred = as.character(svm_age_preds),
  LogReg_Age_Pred = as.character(glm_age_preds),
  KNN_Age_Pred = as.character(knn_age_preds),
  KMeans_Cluster = as.character(k2_clusters[test_idx]),
  SVM_Clust_Pred = as.character(svm_clust_preds),
  PAM_Cluster = as.character(pam_k2_clusters[test_idx]),
  LogReg_Clust_Pred = as.character(glm_clust_preds),
  Spectral_Cluster = as.character(spec_k2_clusters[test_idx]),
  KNN_Clust_Pred = as.character(knn_clust_preds),
  stringsAsFactors = FALSE
)

# Compare age predictions directly. Cluster numbers differ by method and
# cannot be counted as equivalent votes for "cluster 2".
sample_model_matrix$Models_Predicting_Age_12 <- rowSums(
  sample_model_matrix[, c("SVM_Age_Pred", "LogReg_Age_Pred", "KNN_Age_Pred")] == "12"
)
write.csv(sample_model_matrix, "results/sample_model_predictions_matrix.csv", row.names = FALSE)

# Paired significance tests of correctness on the SAME held-out samples
age_predictions <- list(SVM = svm_age_preds, LogReg = glm_age_preds, KNN = knn_age_preds)
age_correct <- lapply(age_predictions, function(pred) as.character(pred) == as.character(age_group[test_idx]))
model_pairs <- combn(names(age_correct), 2, simplify = FALSE)
stat_res_df <- do.call(rbind, lapply(model_pairs, function(pair) {
  a <- age_correct[[pair[1]]]
  b <- age_correct[[pair[2]]]
  discordant_ab <- sum(a & !b)
  discordant_ba <- sum(!a & b)
  # Exact McNemar (binomial) is valid even with few discordant samples
  p <- if (discordant_ab + discordant_ba == 0) 1 else
    binom.test(discordant_ab, discordant_ab + discordant_ba, p = 0.5)$p.value
  data.frame(Comparison = paste(pair, collapse = " vs "), Method = "Exact McNemar",
             Discordant_A_Correct = discordant_ab, Discordant_B_Correct = discordant_ba,
             Raw_P_Value = p)
}))
stat_res_df$Adjusted_P_Value <- p.adjust(stat_res_df$Raw_P_Value, method = "BH")
write.csv(stat_res_df, "results/model_pairwise_significance.csv", row.names = FALSE)

# Held-out accuracy for each prediction task; cluster targets differ by model
prediction_tasks <- list(
  list("SVM", "Age", svm_age_preds, age_group[test_idx]),
  list("LogReg", "Age", glm_age_preds, age_group[test_idx]),
  list("KNN", "Age", knn_age_preds, age_group[test_idx]),
  list("SVM", "KMeans", svm_clust_preds, k2_clusters[test_idx]),
  list("LogReg", "PAM", glm_clust_preds, pam_k2_clusters[test_idx]),
  list("KNN", "Spectral", knn_clust_preds, spec_k2_clusters[test_idx])
)
accuracy_df <- do.call(rbind, lapply(prediction_tasks, function(z) {
  data.frame(Model = z[[1]], Target = z[[2]],
             Test_Accuracy = mean(as.character(z[[3]]) == as.character(z[[4]])))
}))
write.csv(accuracy_df, "results/model_heldout_accuracy.csv", row.names = FALSE)

# Different Gene Counts
gene_subsets <- c(10, 100, 1000, 10000)
perf_results <- list()

for (n in gene_subsets) {
  # Use all available genes, not just the original 5,000
  actual_n <- min(n, length(sorted_genes))
  selected_genes <- head(sorted_genes, actual_n)
  sub_expr <- t(log_counts_sub[selected_genes, , drop = FALSE])
  sub_train <- sub_expr[train_idx, , drop = FALSE]
  sub_test <- sub_expr[test_idx, , drop = FALSE]
  cat("Training models with", actual_n, "genes\n")

  # SVM (held-out AUC)
  set.seed(42)
  m_svm <- svm(sub_train, age_group[train_idx], probability = TRUE, kernel = "radial")
  pr_svm <- attr(predict(m_svm, sub_test, probability = TRUE), "probabilities")[, "12"]
  auc_svm <- as.numeric(auc(roc(as.numeric(age_group[test_idx] == "12"), pr_svm, quiet = TRUE, direction = "<")))

  # Logistic Regression predicting PAM clusters (training-only PCA)
  nonconstant <- apply(sub_train, 2, sd) > 0
  pca_sub <- prcomp(sub_train[, nonconstant, drop = FALSE], scale. = TRUE)
  n_pcs <- min(5, ncol(pca_sub$x), nrow(sub_train) - 2)
  lr_dat <- as.data.frame(pca_sub$x[, seq_len(n_pcs), drop = FALSE])
  lr_dat$Cl <- as.numeric(pam_k2_clusters[train_idx] == "2")
  m_glm <- glm(Cl ~ ., data = lr_dat, family = binomial)
  lr_test <- as.data.frame(predict(pca_sub, newdata = sub_test[, nonconstant, drop = FALSE])[, seq_len(n_pcs), drop = FALSE])
  pr_glm <- predict(m_glm, newdata = lr_test, type = "response")
  auc_glm <- as.numeric(auc(roc(as.numeric(pam_k2_clusters[test_idx] == "2"), pr_glm, quiet = TRUE, direction = "<")))

  # KNN predicting spectral clusters (held-out predictions)
  set.seed(42)
  pr_knn <- knn(train = sub_train, test = sub_test, cl = spec_k2_clusters[train_idx], k = 5, prob = TRUE)
  # attr(prob) is confidence of predicted class; convert to probability of class 2
  pr_knn_2 <- ifelse(as.character(pr_knn) == "2", attr(pr_knn, "prob"), 1 - attr(pr_knn, "prob"))
  auc_knn <- as.numeric(auc(roc(as.numeric(spec_k2_clusters[test_idx] == "2"), pr_knn_2, quiet = TRUE, direction = "<")))

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

# Extract signature genes using correct gene names (column names)
set.seed(42)
svm_linear_model <- svm(expr_5k, k2_clusters, kernel = "linear", scale = TRUE)
svm_weights <- drop(t(svm_linear_model$coefs) %*% svm_linear_model$SV)
names(svm_weights) <- colnames(expr_5k)
svm_top_genes <- head(names(sort(abs(svm_weights), decreasing = TRUE)), 25)

set.seed(42)
cv_lasso <- cv.glmnet(expr_5k, as.numeric(pam_k2_clusters) - 1, family = "binomial", alpha = 1)
lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
lasso_vec <- lasso_coefs[-1, 1]
lasso_top_genes <- names(sort(abs(lasso_vec[lasso_vec != 0]), decreasing = TRUE))

# Report overlap between selected SVM and LASSO signatures
lasso_top25 <- head(lasso_top_genes, 25)
signature_overlap <- intersect(svm_top_genes, lasso_top25)
write.csv(data.frame(Gene = signature_overlap), "results/signature_gene_overlap.csv", row.names = FALSE)
write.csv(data.frame(Method = c("SVM", "LASSO", "Overlap"),
                     Gene_Count = c(length(svm_top_genes), length(lasso_top25), length(signature_overlap))),
          "results/signature_overlap_summary.csv", row.names = FALSE)


predictive_signature_genes <- unique(c(svm_top_genes, lasso_top25))
if (length(predictive_signature_genes) < 10) {
  predictive_signature_genes <- unique(c(predictive_signature_genes, head(colnames(expr_5k), 30)))
}
heatmap_sig_matrix <- t(expr_5k[, predictive_signature_genes, drop = FALSE])
# Gene symbols were already mapped before building expr_5k

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
