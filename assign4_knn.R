library(class)
library(kernlab)
library(caret)
library(pROC)
library(org.Mm.eg.db)

# output directory
output_dir <- "results/knn"

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

#load data
counts <- read.delim(
  file.choose(),
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

metadata <- read.delim(
  file.choose(),
  header = TRUE,
  row.names = 1
)


#clean ids
clean_ids <- gsub("\\..*", "", rownames(counts))

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = clean_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rownames(counts) <- make.unique(
  ifelse(
    is.na(gene_map),
    rownames(counts),
    gene_map
  )
)

#log transformation
log_counts <- log2(counts + 1)

#grouping
samples_keep <- metadata$refinebio_age %in% c("2", "12")

log_counts_sub <- log_counts[, samples_keep]

meta_sub <- metadata[samples_keep, ]

age_group <- factor(
  meta_sub$refinebio_age,
  levels = c("2", "12")
)

#top variable genes
gene_vars <- apply(
  log_counts_sub,
  1,
  var
)

ordered_genes <- names(
  sort(
    gene_vars,
    decreasing = TRUE
  )
)

#knn function
run_knn <- function(n_genes){
  
  selected_genes <- ordered_genes[
    1:min(n_genes, length(ordered_genes))
  ]
  
  expr <- t(
    log_counts_sub[selected_genes, ]
  )
  
  expr <- scale(expr)
  
  set.seed(42)
  
  train_index <- createDataPartition(
    age_group,
    p = 0.70,
    list = FALSE
  )
  
  train_x <- expr[train_index, ]
  test_x <- expr[-train_index, ]
  
  train_y <- age_group[train_index]
  test_y <- age_group[-train_index]
  
  pred <- knn(
    train = train_x,
    test = test_x,
    cl = train_y,
    k = 5,
    prob = TRUE
  )
  
  vote_prob <- attr(pred, "prob")
  
  prob_age12 <- ifelse(
    pred == "12",
    vote_prob,
    1 - vote_prob
  )
  
  roc_obj <- roc(
    response = test_y,
    predictor = prob_age12,
    levels = c("2", "12"),
    quiet = TRUE
  )
  
  auc_value <- as.numeric(
    auc(roc_obj)
  )
  
  cm <- confusionMatrix(
    pred,
    test_y
  )
  
  write.csv(
    as.data.frame(cm$table),
    file.path(
      output_dir,
      paste0(
        "confusion_matrix_",
        n_genes,
        "_genes.csv"
      )
    ),
    row.names = FALSE
  )
  
  return(
    list(
      auc = auc_value,
      roc = roc_obj,
      genes = selected_genes
    )
  )
}

#run model on different number of genes
age_10 <- run_knn(10)
age_100 <- run_knn(100)
age_1000 <- run_knn(1000)
age_5000 <- run_knn(5000)
age_10000 <- run_knn(10000)


#AUC summary 
age_auc_summary <- data.frame(
  Genes = c(
    10,
    100,
    1000,
    5000,
    10000
  ),
  AUC = c(
    age_10$auc,
    age_100$auc,
    age_1000$auc,
    age_5000$auc,
    age_10000$auc
  )
)

print(age_auc_summary)

write.csv(
  age_auc_summary,
  file.path(
    output_dir,
    "age_auc_summary.csv"
  ),
  row.names = FALSE
)


#ROC plots
png(
  file.path(
    output_dir,
    "age_roc_curves.png"
  ),
  width = 1200,
  height = 900
)

plot(
  age_10$roc,
  col = "red",
  lwd = 2,
  main = "ROC Curves: Age Prediction"
)

plot(age_100$roc,
     add = TRUE,
     col = "blue",
     lwd = 2)

plot(age_1000$roc,
     add = TRUE,
     col = "green",
     lwd = 2)

plot(age_5000$roc,
     add = TRUE,
     col = "purple",
     lwd = 2)

plot(age_10000$roc,
     add = TRUE,
     col = "orange",
     lwd = 2)

abline(
  a = 0,
  b = 1,
  lty = 2
)

legend(
  "bottomright",
  legend = c(
    "10 Genes",
    "100 Genes",
    "1000 Genes",
    "5000 Genes",
    "10000 Genes"
  ),
  col = c(
    "red",
    "blue",
    "green",
    "purple",
    "orange"
  ),
  lwd = 2
)

dev.off()


#spectral clustering

expr_5000 <- t(
  log_counts_sub[
    ordered_genes[
      1:min(5000, length(ordered_genes))
    ],
  ]
)

set.seed(42)

spectral_clusters <- factor(
  specc(
    expr_5000,
    centers = 2
  )
)

cluster_results <- data.frame(
  Sample = colnames(log_counts_sub),
  Age_Group = age_group,
  Spectral_Cluster = spectral_clusters
)

print(
  table(
    Age_Group = age_group,
    Spectral_Cluster = spectral_clusters
  )
)

write.csv(
  cluster_results,
  file.path(
    output_dir,
    "spectral_clusters.csv"
  ),
  row.names = FALSE
)

#knn on spectral clusters
expr_cluster <- scale(expr_5000)

set.seed(42)

cluster_train <- createDataPartition(
  spectral_clusters,
  p = 0.70,
  list = FALSE
)

cluster_pred <- knn(
  train = expr_cluster[cluster_train, ],
  test = expr_cluster[-cluster_train, ],
  cl = spectral_clusters[cluster_train],
  k = 5
)

cluster_cm <- confusionMatrix(
  cluster_pred,
  spectral_clusters[-cluster_train]
)

print(cluster_cm)

write.csv(
  as.data.frame(cluster_cm$table),
  file.path(
    output_dir,
    "cluster_prediction_confusion_matrix.csv"
  ),
  row.names = FALSE
)

