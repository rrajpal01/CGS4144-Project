# Load required packages
library(e1071)
library(org.Mm.eg.db)

# Part 1
# Load Data & Preprocessing
counts <- read.delim(file.choose(), header = TRUE, row.names = 1, check.names = FALSE)
metadata <- read.delim(file.choose(), header = TRUE, row.names = 1)

# Clean Ensembl IDs, map symbols, and log-transform
clean_ensembl_ids <- gsub("\\..*", "", rownames(counts))
gene_map <- mapIds(org.Mm.eg.db, keys = clean_ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(counts) <- make.unique(ifelse(is.na(gene_map), rownames(counts), gene_map))
log_counts <- log2(counts + 1)

# Subset to the two working sample groups
samples_2groups <- metadata$refinebio_age %in% c("2", "12")
log_counts_sub <- log_counts[, samples_2groups]
meta_sub <- metadata[samples_2groups, ]
age_group <- factor(meta_sub$refinebio_age, levels = c("2", "12"))

# Part 2
# Feature Selection (Top 5,000 Genes)
gene_vars <- apply(log_counts_sub, 1, var)
sorted_genes <- names(sort(gene_vars, decreasing = TRUE))
top5000_matrix <- log_counts_sub[sorted_genes[1:min(5000, length(sorted_genes))], ]
expr_5k <- t(top5000_matrix)

# SVM Predicting Assignment 1 Groups (Age)
svm_age_data <- as.data.frame(expr_5k)
svm_age_data$Age <- age_group

set.seed(42)
svm_age_model <- svm(Age ~ ., data = svm_age_data, probability = TRUE, kernel = "radial")
svm_age_preds <- predict(svm_age_model, svm_age_data)
print(table(True_Age = age_group, Predicted_SVM_Age = svm_age_preds))

# SVM Predicting Assignment 3 K-Means Clusters
set.seed(42)
k2_clusters <- factor(kmeans(expr_5k, centers = 2, nstart = 25)$cluster)

svm_cluster_data <- as.data.frame(expr_5k)
svm_cluster_data$Cluster <- k2_clusters

set.seed(42)
svm_cluster_model <- svm(Cluster ~ ., data = svm_cluster_data, probability = TRUE, kernel = "radial")
svm_cluster_preds <- predict(svm_cluster_model, svm_cluster_data)
print(table(True_KMeans = k2_clusters, Predicted_SVM_Cluster = svm_cluster_preds))