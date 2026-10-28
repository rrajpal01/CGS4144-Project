# Part 1
#Install Packages

library(org.Mm.eg.db)
library(ggplot2)
library(glmnet)
library(pROC)
library(cluster)

# make res folder
output_dir <- "results/log_reg"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load data
counts <- read.delim("data/SRP082327/SRP082327.tsv",
                     header = TRUE, row.names = 1,
                     check.names = FALSE)

metadata <- read.delim("data/SRP082327/metadata_SRP082327.tsv",
                       header = TRUE, row.names = 1)

# Convert Ensembl IDs to gene symbols
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

stopifnot(all(colnames(counts) %in% rownames(metadata)))
metadata <- metadata[colnames(counts), , drop = FALSE]

# Log transform
log_counts <- log2(as.matrix(counts) + 1)

samples_2groups <- !is.na(metadata$refinebio_age) &
  as.character(metadata$refinebio_age) %in% c("2", "12")

log_counts_sub <- log_counts[, samples_2groups, drop = FALSE]
meta_sub <- metadata[samples_2groups, , drop = FALSE]

age_group <- factor(as.character(meta_sub$refinebio_age),
                    levels = c("2", "12"))




set.seed(42)

train_idx <- unlist(lapply(
  split(seq_along(age_group), age_group),
  function(idx) sample(idx, floor(0.8 * length(idx)))
))


gene_vars <- apply(log_counts_sub[, train_idx, drop = FALSE], 1, var)
sorted_genes <- names(sort(gene_vars, decreasing = TRUE))

# do top 5k
top5000_matrix <- log_counts_sub[
  head(sorted_genes, 5000), , drop = FALSE
]

expr_5k <- t(top5000_matrix)

X_train <- expr_5k[train_idx, , drop = FALSE]
X_test <- expr_5k[-train_idx, , drop = FALSE]

y_train <- age_group[train_idx]
y_test <- age_group[-train_idx]


# do log reg

set.seed(42)

logistic_age_model <- cv.glmnet(
  x = X_train,
  y = y_train,
  family = "binomial",
  alpha = 0,
  nfolds = 5
)

logistic_age_probs <- as.vector(predict(
  logistic_age_model,
  newx = X_test,
  s = "lambda.min",
  type = "response"
))

logistic_age_preds <- factor(
  ifelse(logistic_age_probs >= 0.5, "12", "2"),
  levels = c("2", "12")
)

age_accuracy <- mean(logistic_age_preds == y_test)

age_roc <- roc(
  response = y_test,
  predictor = logistic_age_probs,
  levels = c("2", "12"),
  direction = "<",
  quiet = TRUE
)

age_auc <- as.numeric(auc(age_roc))


p_age <- ggroc(age_roc) +
  ggtitle("Logistic Regression ROC - Age Groups") +
  theme_minimal()

ggsave(file.path(output_dir, "logistic_age_roc.png"),
       plot = p_age, width = 7, height = 5)


# compare to pam

# do pam again
pam_gene_vars <- apply(log_counts_sub, 1, var)
pam_sorted_genes <- names(sort(pam_gene_vars, decreasing = TRUE))

pam_expr_5k <- t(log_counts_sub[
  head(pam_sorted_genes, 5000), , drop = FALSE
])

set.seed(42)
pam_result <- pam(pam_expr_5k, k = 2)
pam_clusters <- factor(pam_result$clustering)


y_pam_train <- pam_clusters[train_idx]
y_pam_test <- pam_clusters[-train_idx]


set.seed(42)

logistic_pam_model <- cv.glmnet(
  x = X_train,
  y = y_pam_train,
  family = "binomial",
  alpha = 0,
  nfolds = 5
)

logistic_pam_probs <- as.vector(predict(
  logistic_pam_model,
  newx = X_test,
  s = "lambda.min",
  type = "response"
))

logistic_pam_preds <- factor(
  ifelse(logistic_pam_probs >= 0.5, "2", "1"),
  levels = c("1", "2")
)

pam_accuracy <- mean(logistic_pam_preds == y_pam_test)

pam_roc <- roc(
  response = y_pam_test,
  predictor = logistic_pam_probs,
  levels = c("1", "2"),
  direction = "<",
  quiet = TRUE
)

pam_auc <- as.numeric(auc(pam_roc))


p_pam <- ggroc(pam_roc) +
  ggtitle("Logistic Regression ROC - PAM Clusters") +
  theme_minimal()

ggsave(file.path(output_dir, "logistic_pam_roc.png"),
       plot = p_pam, width = 7, height = 5)


# gene cnt compare

gene_counts <- c(10, 100, 1000, 10000)

results_age <- data.frame(Genes = integer(), AUC = numeric())
results_pam <- data.frame(Genes = integer(), AUC = numeric())

for (n in gene_counts) {


  selected_genes <- head(sorted_genes, n)

  X_train_n <- t(log_counts_sub[
    selected_genes, train_idx, drop = FALSE
  ])

  X_test_n <- t(log_counts_sub[
    selected_genes, -train_idx, drop = FALSE
  ])


  for (target in c("Age", "PAM")) {

    train_labels <- if (target == "Age") y_train else y_pam_train
    test_labels <- if (target == "Age") y_test else y_pam_test
    class_levels <- if (target == "Age") c("2", "12") else c("1", "2")

    set.seed(42)

    model_n <- cv.glmnet(
      x = X_train_n,
      y = train_labels,
      family = "binomial",
      alpha = 0,
      nfolds = 5
    )

    probs_n <- as.vector(predict(
      model_n,
      newx = X_test_n,
      s = "lambda.min",
      type = "response"
    ))

    roc_n <- roc(
      response = test_labels,
      predictor = probs_n,
      levels = class_levels,
      direction = "<",
      quiet = TRUE
    )

    result <- data.frame(
      Genes = n,
      AUC = as.numeric(auc(roc_n))
    )

    if (target == "Age") {
      results_age <- rbind(results_age, result)
    } else {
      results_pam <- rbind(results_pam, result)
    }
  }
}




comparison_df <- rbind(
  data.frame(
    Genes = c(results_age$Genes, 5000),
    AUC = c(results_age$AUC, age_auc),
    Model = "Age"
  ),
  data.frame(
    Genes = c(results_pam$Genes, 5000),
    AUC = c(results_pam$AUC, pam_auc),
    Model = "PAM"
  )
)

comparison_df <- comparison_df[
  order(comparison_df$Model, comparison_df$Genes), ]

p_auc <- ggplot(comparison_df,
                aes(x = Genes, y = AUC, color = Model)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_log10(breaks = c(10, 100, 1000, 5000, 10000)) +
  labs(
    title = "Logistic Regression AUC by Number of Genes",
    x = "Number of Genes",
    y = "ROC AUC",
    color = "Prediction Target"
  ) +
  theme_minimal()

ggsave(file.path(output_dir, "logistic_auc_comparison.png"),
       plot = p_auc, width = 8, height = 5)


performance_summary <- data.frame(
  Model = c("Age", "PAM"),
  Accuracy = c(age_accuracy, pam_accuracy),
  AUC = c(age_auc, pam_auc)
)

write.csv(
  performance_summary,
  file.path(output_dir, "logistic_performance.csv"),
  row.names = FALSE
)


write.csv(
  comparison_df,
  file.path(output_dir, "logistic_gene_comparison.csv"),
  row.names = FALSE
)


predictions_df <- data.frame(
  Sample = rownames(X_test),
  True_Age = as.character(y_test),
  Predicted_Age = as.character(logistic_age_preds),
  Age_Probability_12 = logistic_age_probs,
  True_PAM = as.character(y_pam_test),
  Predicted_PAM = as.character(logistic_pam_preds),
  PAM_Probability_2 = logistic_pam_probs
)

write.csv(
  predictions_df,
  file.path(output_dir, "logistic_predictions.csv"),
  row.names = FALSE
)
