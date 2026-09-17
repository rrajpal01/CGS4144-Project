if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!("org.Mm.eg.db" %in% installed.packages())) {
  BiocManager::install("org.Mm.eg.db", update = FALSE)
}

#Attach the library
library(org.Mm.eg.db)

# We will need this so we can use the pipe: %>%
library(magrittr)

#Read in metadata TSV file
metadata <- readr::read_tsv("C:/Jasmine/UF Courses/FALL 2026/CGS4144/project/data/SRP082327/metadata_SRP082327.tsv")

# Read in data TSV file
expression_df <- readr::read_tsv("C:/Jasmine/UF Courses/FALL 2026/CGS4144/project/data/SRP082327/SRP082327.tsv") %>%
  # Tuck away the Gene ID column as row names
  tibble::column_to_rownames("Gene")

# Make the data in the order of the metadata
expression_df <- expression_df %>%
  dplyr::select(metadata$refinebio_accession_code)

# Check if this is in the same order
all.equal(colnames(expression_df), metadata$refinebio_accession_code)

# Bring back the "Gene" column in preparation for mapping
expression_df <- expression_df %>%
  tibble::rownames_to_column("Gene")

# Map Ensembl IDs to their associated Entrez IDs
mapped_list <- mapIds(
  org.Mm.eg.db, # Replace with annotation package for your organism
  keys = expression_df$Gene,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "list"
)

head(mapped_list)

# Let's make our list a bit more manageable by turning it into a data frame
mapped_df <- mapped_list %>%
  tibble::enframe(name = "Ensembl", value = "Symbol") %>%
  # enframe() makes a `list` column; we will simplify it with unnest()
  # This will result in one row of our data frame per list item
  tidyr::unnest(cols = Symbol)

head(mapped_df)

# Use the `summary()` function to show the distribution of Entrez values
# We need to use `as.factor()` here to get the count of unique values
# `maxsum = 10` limits the summary to 10 distinct values
summary(as.factor(mapped_df$Symbol), maxsum = 10)

final_mapped_df <- mapped_df %>%
  # Join to expression data by Ensembl ID
  dplyr::inner_join(expression_df, by = c("Ensembl" = "Gene")) %>%
  # Optional: drop the Ensembl column now that you have Symbol
  dplyr::select(-Ensembl)

head(final_mapped_df)

library(dplyr)

# Assuming final_mapped_df has: Symbol, and then all your SRR sample columns
# First, drop rows with no gene symbol (NA) - these can't be interpreted biologically
final_mapped_df <- final_mapped_df %>%
  dplyr::filter(!is.na(Symbol))

# Calculate a summary expression value per row (mean across all samples)
# This lets us rank duplicate Symbol rows by overall expression level
final_mapped_df <- final_mapped_df %>%
  dplyr::mutate(
    mean_expression = rowMeans(dplyr::select(., -Symbol), na.rm = TRUE)
  )

# For each Symbol, keep only the row with the highest mean expression
collapsed_df <- final_mapped_df %>%
  dplyr::group_by(Symbol) %>%
  dplyr::slice_max(order_by = mean_expression, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(-mean_expression)  # drop the helper column now that we're done

# Confirm no more duplicates
sum(duplicated(collapsed_df$Symbol))  # should be 0

# Now Symbol is unique - safe to use as row names for downstream tools like DESeq2
final_matrix <- collapsed_df %>%
  tibble::column_to_rownames("Symbol")

head(final_matrix)

final_matrix %>%
  tibble::rownames_to_column("Symbol") %>%
  readr::write_tsv(file.path("C:/Jasmine/UF Courses/FALL 2026/CGS4144/project/data/SRP082327", 
                             "SRP082327_gene_symbols_collapsed.tsv"))


