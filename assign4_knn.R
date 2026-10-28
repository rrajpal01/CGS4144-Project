# =====================================================
# Required Packages
# =====================================================

library(class)
library(kernlab)
library(pROC)
library(caret)
library(org.Mm.eg.db)

# =====================================================
# Load Data
# =====================================================

counts <- read.delim(file.choose(),
                     header = TRUE,
                     row.names = 1,
                     check.names = FALSE)

metadata <- read.delim(file.choose(),
                       header = TRUE,
                       row.names = 1)

# =====================================================
# Gene Mapping
# =====================================================

clean_ensembl_ids <- gsub("\\..*", "", rownames(counts))

gene_map <- mapIds(
  org.Mm.eg.db,
  keys = clean_ensembl_ids,
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rownames(counts) <- make.unique(
  ifelse(is.na(gene_map),
         rownames(counts),
         gene_map)
)

# =====================================================
# Log Transform
# =====================================================

log_counts <- log2(counts + 1)

# =====================================================
# Assignment 1 groups
# =====================================================

samples_keep <- metadata$refinebio_age %in% c("2","12")

log_counts_sub <- log_counts[,samples_keep]

meta_sub <- metadata[samples_keep,]

age_group <- factor(
  meta_sub$refinebio_age,
  levels = c("2","12")
)

# =====================================================
# Top Variable Genes
# =====================================================

gene_vars <- apply(log_counts_sub,1,var)

ordered_genes <- names(
  sort(gene_vars,
       decreasing=TRUE)
)

# =====================================================
# Function for kNN Models
# =====================================================

run_knn_model <- function(num_genes){
  
  selected <- ordered_genes[
    1:min(num_genes,length(ordered_genes))
  ]
  
  expression_matrix <- t(
    log_counts_sub[selected,]
  )
  
  expression_matrix <- scale(expression_matrix)
  
  df <- as.data.frame(expression_matrix)
  
  set.seed(42)
  
  train_idx <- createDataPartition(
    age_group,
    p=.70,
    list=FALSE
  )
  
  train_x <- df[train_idx,]
  test_x <- df[-train_idx,]
  
  train_y <- age_group[train_idx]
  test_y <- age_group[-train_idx]
  
  pred <- knn(
    train=train_x,
    test=test_x,
    cl=train_y,
    k=5,
    prob=TRUE
  )
  
  probs <- attr(pred,"prob")
  
  auc <- roc(
    response=test_y,
    predictor=probs
  )
  
  cat("\n==========================\n")
  cat("Genes:",num_genes,"\n")
  cat("AUC:",auc$auc,"\n")
  
  print(
    confusionMatrix(
      pred,
      test_y
    )
  )
  
  return(
    list(
      auc=auc$auc,
      genes=selected,
      predictions=pred
    )
  )
}

# =====================================================
# Run Multiple Gene Sets
# =====================================================

results_10 <- run_knn_model(10)

results_100 <- run_knn_model(100)

results_1000 <- run_knn_model(1000)

results_5000 <- run_knn_model(5000)

results_10000 <- run_knn_model(10000)

# =====================================================
# Spectral Clustering
# Assignment 3
# =====================================================

expr_5000 <- t(
  log_counts_sub[
    ordered_genes[1:min(5000,length(ordered_genes))],
  ]
)

set.seed(42)

spectral_clusters <- factor(
  specc(
    expr_5000,
    centers=2
  )
)

table(
  Age=age_group,
  Spectral=spectral_clusters
)

#===============================================

library(pheatmap)

signature_genes <- results_100$genes

heatmap_matrix <- log_counts_sub[
  signature_genes,
]

annotation_df <- data.frame(
  Age = age_group
)

rownames(annotation_df) <- colnames(heatmap_matrix)

pheatmap(
  heatmap_matrix,
  scale="row",
  annotation_col=annotation_df,
  clustering_method="complete",
  fontsize_row=5,
  main="Predictive Gene Signature Heatmap"
)
