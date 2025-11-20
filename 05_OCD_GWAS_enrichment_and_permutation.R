

# Load required libraries
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(tidyr)

# Set options for better printing
options(scipen = 999)

################################################################################
# SECTION 1: HELPER FUNCTIONS
################################################################################

# Function to extract tissue name from filename
get_tissue_name <- function(filename) {
  tissue <- gsub("gene_tpm_v10_brain_", "", filename)
  tissue <- gsub(".gct.gz", "", tissue)
  tissue <- gsub("_", " ", tissue)
  tissue <- tools::toTitleCase(tissue)
  return(tissue)
}

# Function to read GTEx data with consistent filtering
read_gtex_data <- function(filename, apply_median_filter = TRUE, median_threshold = 0.1) {
  
  cat("Reading", filename, "...\n")
  
  # Read the data
  data <- fread(filename, skip = 2, showProgress = FALSE)
  gene_ids <- data[[1]]
  gene_names <- data[[2]]
  expr_matrix <- as.matrix(data[, 3:ncol(data)])
  
  # Get sample size
  n_samples <- ncol(expr_matrix)
  
  if (apply_median_filter) {
    # Apply median TPM filter
    median_tpm <- apply(expr_matrix, 1, median)
    keep_idx <- which(median_tpm >= median_threshold)
    
    # Always keep KNCN and MKNK1-AS1 if they exist
    kncn_idx <- which(gene_names == "KNCN" | grepl("KNCN", gene_ids))[1]
    mknk_idx <- which(gene_names == "MKNK1-AS1" | grepl("MKNK1-AS1", gene_ids))[1]
    
    if (!is.na(kncn_idx)) keep_idx <- union(keep_idx, kncn_idx)
    if (!is.na(mknk_idx)) keep_idx <- union(keep_idx, mknk_idx)
    
    # Apply filter
    expr_matrix <- expr_matrix[keep_idx, ]
    gene_names <- gene_names[keep_idx]
    gene_ids <- gene_ids[keep_idx]
    
    cat("  Retained", length(keep_idx), "genes after filtering (median TPM >=", median_threshold, ")\n")
  }
  
  return(list(
    gene_ids = gene_ids,
    gene_names = gene_names,
    expr_matrix = expr_matrix,
    n_samples = n_samples,
    tissue = get_tissue_name(filename)
  ))
}

################################################################################
# SECTION 2: CORRELATION ANALYSIS WITH VISUALIZATION
################################################################################

# Enhanced correlation plot function
create_correlation_plot <- function(data_obj, gene1 = "KNCN", gene2 = "MKNK1-AS1", panel_label = NULL) {
  
  # Extract data
  gene_names <- data_obj$gene_names
  gene_ids <- data_obj$gene_ids
  expr_matrix <- data_obj$expr_matrix
  n_samples <- data_obj$n_samples
  tissue <- data_obj$tissue
  
  # Find gene indices
  gene1_idx <- which(gene_names == gene1 | grepl(gene1, gene_ids))[1]
  gene2_idx <- which(gene_names == gene2 | grepl(gene2, gene_ids))[1]
  
  if (is.na(gene1_idx) || is.na(gene2_idx)) {
    cat("  Warning: Could not find", gene1, "or", gene2, "in", tissue, "\n")
    return(NULL)
  }
  
  # Extract and log-transform expression
  gene1_expr <- log2(as.numeric(expr_matrix[gene1_idx, ]) + 1)
  gene2_expr <- log2(as.numeric(expr_matrix[gene2_idx, ]) + 1)
  
  # Create data frame
  plot_data <- data.frame(gene1 = gene1_expr, gene2 = gene2_expr)
  
  # Pearson correlation
  cor_test <- cor.test(gene1_expr, gene2_expr, method = "pearson")
  
  # Create plot title
  plot_title <- if (!is.null(panel_label)) paste0("(", panel_label, ") ", tissue) else tissue
  
  # Create the plot
  p <- ggplot(plot_data, aes(x = gene1, y = gene2)) +
    geom_point(alpha = 0.5, size = 1, color = "darkblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", linewidth = 0.8) +
    labs(
      title = plot_title,
      x = expression(log[2]("KNCN TPM + 1")),
      y = expression(log[2]("MKNK1-AS1 TPM + 1"))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7),
      panel.grid.minor = element_blank()
    ) +
    annotate("text", 
             x = Inf, y = -Inf, 
             label = sprintf("r = %.3f\np = %.3e\nn = %d", 
                             cor_test$estimate, 
                             cor_test$p.value,
                             n_samples),
             hjust = 1.1, vjust = -0.5, size = 3)
  
  return(p)
}

################################################################################
# SECTION 3: HIGH CORRELATION GENE IDENTIFICATION
################################################################################

# Function to compute correlations with significance
compute_correlations_with_significance <- function(data_obj, 
                                                   cor_threshold = 0.7, 
                                                   p_threshold = 0.05,
                                                   return_all = FALSE) {
  
  gene_names <- data_obj$gene_names
  expr_matrix <- data_obj$expr_matrix
  
  # Find KNCN and MKNK1-AS1 indices
  kncn_idx <- which(gene_names == "KNCN")[1]
  mknk_idx <- which(gene_names == "MKNK1-AS1")[1]
  
  if (is.na(kncn_idx) || is.na(mknk_idx)) {
    return(NULL)
  }
  
  # Extract and log-transform expression
  kncn_expr <- log2(as.numeric(expr_matrix[kncn_idx, ]) + 1)
  mknk_expr <- log2(as.numeric(expr_matrix[mknk_idx, ]) + 1)
  
  # Initialize results
  gene_cor_table <- data.frame(
    Gene = gene_names,
    Cor_KNCN = NA,
    P_KNCN = NA,
    Cor_MKNK1AS1 = NA,
    P_MKNK1AS1 = NA
  )
  
  # Compute correlations
  cat("  Computing correlations for", length(gene_names), "genes...\n")
  
  for (i in 1:nrow(expr_matrix)) {
    gene_expr <- log2(as.numeric(expr_matrix[i, ]) + 1)
    
    # KNCN correlation
    if (sd(gene_expr) > 0 && sd(kncn_expr) > 0) {
      cor_test_kncn <- cor.test(kncn_expr, gene_expr, method = "pearson")
      gene_cor_table$Cor_KNCN[i] <- cor_test_kncn$estimate
      gene_cor_table$P_KNCN[i] <- cor_test_kncn$p.value
    }
    
    # MKNK1-AS1 correlation
    if (sd(gene_expr) > 0 && sd(mknk_expr) > 0) {
      cor_test_mknk <- cor.test(mknk_expr, gene_expr, method = "pearson")
      gene_cor_table$Cor_MKNK1AS1[i] <- cor_test_mknk$estimate
      gene_cor_table$P_MKNK1AS1[i] <- cor_test_mknk$p.value
    }
  }
  
  if (return_all) {
    return(gene_cor_table)
  }
  
  # Filter by correlation threshold and significance
  filtered_genes <- gene_cor_table %>%
    filter(!is.na(Cor_KNCN) & !is.na(Cor_MKNK1AS1)) %>%
    filter(abs(Cor_KNCN) >= cor_threshold & P_KNCN < p_threshold) %>%
    filter(abs(Cor_MKNK1AS1) >= cor_threshold & P_MKNK1AS1 < p_threshold)
  
  return(filtered_genes)
}

################################################################################
# SECTION 4: MAIN ANALYSIS PIPELINE
################################################################################

# Get all GTEx brain tissue files
files <- list.files(pattern = "gene_tpm_v10_brain_.*\\.gct\\.gz$")
cat("\nFound", length(files), "GTEx brain tissue files\n\n")

# 1. Check sample sizes across tissues
cat("=== SAMPLE SIZE CHECK ===\n")
sample_sizes <- data.frame()

for (file in files) {
  data <- fread(file, skip = 2, showProgress = FALSE)
  n_samples <- ncol(data) - 2  # Subtract gene ID and name columns
  tissue <- get_tissue_name(file)
  
  sample_sizes <- rbind(sample_sizes, data.frame(
    Tissue = tissue,
    N_samples = n_samples,
    File = file,
    stringsAsFactors = FALSE
  ))
}

print(sample_sizes[order(sample_sizes$N_samples), ])
write.csv(sample_sizes, "GTEx_brain_tissue_sample_sizes.csv", row.names = FALSE)

# 2. Create correlation plots for all tissues
cat("\n=== CREATING CORRELATION PLOTS ===\n")
plot_list <- list()
correlation_results <- list()
panel_labels <- letters[1:length(files)]

for (i in 1:length(files)) {
  # Read data with consistent filtering
  data_obj <- read_gtex_data(files[i], apply_median_filter = TRUE)
  
  # Create plot
  p <- create_correlation_plot(data_obj, panel_label = panel_labels[i])
  
  if (!is.null(p)) {
    plot_list[[i]] <- p
    
    # Store correlation results
    gene_names <- data_obj$gene_names
    expr_matrix <- data_obj$expr_matrix
    
    kncn_idx <- which(gene_names == "KNCN")[1]
    mknk_idx <- which(gene_names == "MKNK1-AS1")[1]
    
    if (!is.na(kncn_idx) && !is.na(mknk_idx)) {
      kncn_expr <- log2(as.numeric(expr_matrix[kncn_idx, ]) + 1)
      mknk_expr <- log2(as.numeric(expr_matrix[mknk_idx, ]) + 1)
      cor_test <- cor.test(kncn_expr, mknk_expr, method = "pearson")
      
      correlation_results[[i]] <- data.frame(
        Tissue = data_obj$tissue,
        Correlation = cor_test$estimate,
        P_value = cor_test$p.value,
        N_samples = data_obj$n_samples,
        stringsAsFactors = FALSE
      )
    }
  }
}

# Remove NULL plots and results
plot_list <- plot_list[!sapply(plot_list, is.null)]
correlation_results <- do.call(rbind, correlation_results)

# Create combined figure
n_plots <- length(plot_list)
n_cols <- 4
n_rows <- ceiling(n_plots / n_cols)

combined_plot <- plot_grid(
  plotlist = plot_list,
  ncol = n_cols,
  nrow = n_rows,
  align = "hv"
)

# Save correlation plots
ggsave("KNCN_MKNK1AS1_correlations_all_tissues.pdf", 
       combined_plot, 
       width = 12, 
       height = 10,
       dpi = 600)

png("KNCN_MKNK1AS1_correlations_all_tissues.png", 
    width = 12, 
    height = 10, 
    units = "in", 
    res = 600)
print(combined_plot)
dev.off()

# 3. Create correlation summary bar plot
cat("\n=== CREATING CORRELATION SUMMARY PLOT ===\n")

# Add significance indicators
correlation_results$Significant <- correlation_results$P_value < 0.05
correlation_results$Significance <- cut(correlation_results$P_value, 
                                        breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                                        labels = c("***", "**", "*", ""))

# Create bar plot
cor_bar_plot <- ggplot(correlation_results, 
                       aes(x = reorder(Tissue, Correlation), 
                           y = Correlation, 
                           fill = Correlation > 0.7)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Significance, 
                y = Correlation + ifelse(Correlation > 0, 0.02, -0.02)),
            size = 5) +
  # Add sample size labels
  geom_text(aes(label = paste0("n=", N_samples), 
                y = 0.1),
            size = 3, angle = 90, hjust = 0) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkred"),
                    labels = c("r ≤ 0.7", "r > 0.7")) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Brain Region",
       y = "Pearson Correlation Coefficient",
       fill = "Correlation Strength",
       caption = "Significance: * p < 0.05, ** p < 0.01, *** p < 0.001") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom")

ggsave("KNCN_MKNK1AS1_correlation_barplot.png", 
       cor_bar_plot,
       width = 8, 
       height = 6, 
       dpi = 600)

# Save correlation results
write.csv(correlation_results, "KNCN_MKNK1AS1_correlation_summary.csv", row.names = FALSE)

# 4. Identify high correlation genes for each tissue
cat("\n=== IDENTIFYING HIGH CORRELATION GENES ===\n")

for (file in files) {
  tissue <- get_tissue_name(file)
  cat("\nProcessing", tissue, "\n")
  
  # Read data with filtering
  data_obj <- read_gtex_data(file, apply_median_filter = TRUE)
  
  # Compute correlations
  high_cor_genes <- compute_correlations_with_significance(data_obj, 
                                                           cor_threshold = 0.7,
                                                           p_threshold = 0.05)
  
  if (!is.null(high_cor_genes) && nrow(high_cor_genes) > 0) {
    high_cor_genes$Tissue <- tissue
    output_filename <- paste0("high_cor_genes_", gsub(" ", "_", tissue), ".csv")
    write.csv(high_cor_genes, file = output_filename, row.names = FALSE)
    cat("  Found", nrow(high_cor_genes), "high correlation genes\n")
  } else {
    cat("  No genes passed the correlation threshold\n")
    # Write empty file
    write.csv(data.frame(Gene = character(0),
                         Cor_KNCN = numeric(0),
                         P_KNCN = numeric(0),
                         Cor_MKNK1AS1 = numeric(0),
                         P_MKNK1AS1 = numeric(0),
                         Tissue = character(0)),
              file = paste0("high_cor_genes_", gsub(" ", "_", tissue), "_EMPTY.csv"),
              row.names = FALSE)
  }
}

# 5. Analyze shared genes across tissues
cat("\n=== ANALYZING SHARED HIGH CORRELATION GENES ===\n")

# Read all high correlation gene files
high_cor_files <- list.files(pattern = "^high_cor_genes_.*\\.csv$")
high_cor_files <- high_cor_files[!grepl("_EMPTY\\.csv$", high_cor_files)]

all_high_cor <- list()
for (file in high_cor_files) {
  tissue_data <- read.csv(file, stringsAsFactors = FALSE)
  if (nrow(tissue_data) > 0) {
    all_high_cor[[file]] <- tissue_data
  }
}

# Find genes shared across multiple tissues
if (length(all_high_cor) > 0) {
  gene_tissue_matrix <- table(unlist(lapply(all_high_cor, function(x) x$Gene)))
  shared_genes <- names(gene_tissue_matrix)[gene_tissue_matrix >= 2]
  
  cat("\nGenes found in 2+ tissues:", length(shared_genes), "\n")
  
  # Create detailed shared gene report
  shared_gene_details <- data.frame()
  for (gene in shared_genes) {
    for (tissue_data in all_high_cor) {
      if (gene %in% tissue_data$Gene) {
        gene_row <- tissue_data[tissue_data$Gene == gene, ]
        shared_gene_details <- rbind(shared_gene_details, gene_row)
      }
    }
  }
  
  write.csv(shared_gene_details, "shared_high_correlation_genes_details.csv", row.names = FALSE)
}

################################################################################
# SECTION 5: OCD GENE ENRICHMENT ANALYSIS
################################################################################

cat("\n=== OCD GENE ENRICHMENT ANALYSIS ===\n")

library(readxl)

# Read OCD genes
ocd_genes_pre <- read_excel("../gwas/41588_2025_2189_MOESM3_ESM.xlsx", skip=1, sheet=15)
ocd_genes <- ocd_genes_pre$Gene
head(ocd_genes)

cat("Loaded", length(ocd_genes), "OCD-associated genes\n") ## 251 genes

# Perform enrichment analysis
enrichment_results <- data.frame()

for (file in files) {
  tissue <- get_tissue_name(file)
  cat("\nAnalyzing", tissue, "\n")
  
  # Read corresponding high correlation genes file
  high_cor_file <- paste0("high_cor_genes_", gsub(" ", "_", tissue), ".csv")
  
  if (!file.exists(high_cor_file) || grepl("_EMPTY\\.csv$", high_cor_file)) {
    cat("  No high correlation genes found\n")
    next
  }
  
  # Read high correlation genes
  high_cor_data <- read.csv(high_cor_file, stringsAsFactors = FALSE)
  if (nrow(high_cor_data) == 0) next
  
  high_cor_genes <- unique(high_cor_data$Gene)
  
  # Read original data to get background genes
  data_obj <- read_gtex_data(file, apply_median_filter = TRUE)
  background_genes <- data_obj$gene_names
  
  # Create contingency table
  a <- sum(high_cor_genes %in% ocd_genes)  # OCD genes in high cor
  b <- sum(!(high_cor_genes %in% ocd_genes))  # Non-OCD genes in high cor
  c <- sum((background_genes %in% ocd_genes) & !(background_genes %in% high_cor_genes))  # OCD not in high cor
  d <- sum(!(background_genes %in% ocd_genes) & !(background_genes %in% high_cor_genes))  # Non-OCD not in high cor
  
  cat("  Contingency table: a=", a, ", b=", b, ", c=", c, ", d=", d, "\n")
  
  # Fisher's exact test
  if (a + b > 0) {
    contingency_matrix <- matrix(c(a, b, c, d), nrow = 2)
    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
    
    enrichment_results <- rbind(enrichment_results, data.frame(
      Tissue = tissue,
      OCD_genes_in_high_cor = a,
      Non_OCD_genes_in_high_cor = b,
      OCD_genes_non_in_high_cor = c,
      Non_OCD_genes_non_in_high_cor = d,
      Total_high_cor_genes = length(high_cor_genes),
      Total_OCD_genes_in_background = sum(background_genes %in% ocd_genes),
      Total_genes_in_background = length(background_genes),
      Enrichment_ratio = (a / length(high_cor_genes)) / (sum(background_genes %in% ocd_genes) / length(background_genes)),
      P_value = fisher_result$p.value,
      Odds_ratio = fisher_result$estimate,
      stringsAsFactors = FALSE
    ))
  }
}

# Apply FDR correction
enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
enrichment_results$Significant <- enrichment_results$FDR < 0.05

# Sort by FDR
enrichment_results <- enrichment_results[order(enrichment_results$FDR), ]

# Save results
write.csv(enrichment_results, "OCD_gene_enrichment_results.csv", row.names = FALSE)
write.csv(enrichment_results, "OCD_gene_enrichment_results_251genes.csv", row.names = FALSE)

# Print summary
cat("\n=== ENRICHMENT ANALYSIS SUMMARY ===\n")
cat("Total tissues analyzed:", nrow(enrichment_results), "\n")
cat("Tissues with significant enrichment (FDR < 0.05):", sum(enrichment_results$Significant), "\n\n")

if (sum(enrichment_results$Significant) > 0) {
  cat("Significant tissues:\n")
  sig_results <- enrichment_results[enrichment_results$Significant, ]
  for (i in 1:nrow(sig_results)) {
    cat(sprintf("  %s: %d/%d OCD genes (%.1f%%), OR = %.2f, FDR = %.3e\n",
                sig_results$Tissue[i],
                sig_results$OCD_genes_in_high_cor[i],
                sig_results$Total_high_cor_genes[i],
                100 * sig_results$OCD_genes_in_high_cor[i] / sig_results$Total_high_cor_genes[i],
                sig_results$Odds_ratio[i],
                sig_results$FDR[i]))
  }
}

# Create enrichment visualization
enrichment_plot <- ggplot(enrichment_results, 
                          aes(x = reorder(Tissue, -log10(FDR)), 
                              y = -log10(FDR),
                              fill = Significant)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen")) +
  coord_flip() +
  labs(x = "Brain Region",
       y = "-log10(FDR)",
       title = "OCD Gene Enrichment Across Brain Regions",
       fill = "FDR < 0.05") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"))

ggsave("OCD_gene_enrichment_barplot.png",
       enrichment_plot,
       width = 8,
       height = 6,
       dpi = 600)




################################################################################
# SECTION 5: OCD GENE ENRICHMENT ANALYSIS (Fisher + Permutation)
################################################################################

cat("\n=== OCD GENE ENRICHMENT ANALYSIS (Fisher + Permutation) ===\n")

library(readxl)

## 1) Load OCD GWAS gene list
ocd_genes_pre <- read_excel("../gwas/41588_2025_2189_MOESM3_ESM.xlsx",
                            skip = 1, sheet = 15)
ocd_genes <- unique(ocd_genes_pre$Gene)
cat("Loaded", length(ocd_genes), "OCD-associated genes\n")

## 2) Permutation settings
n_perm <- 10000          # Number of permutations (adjust to e.g. 5000 or 20000 if needed)
set.seed(20241117)       # Fix seed for reproducibility

enrichment_results <- data.frame()
magma_geneset_list <- list()  # List for MAGMA gene-set file (used below in step 2)

tissue_idx <- 0

for (file in files) {
  tissue <- get_tissue_name(file)
  tissue_idx <- tissue_idx + 1
  cat("\nAnalyzing", tissue, "\n")
  
  # File name for high-correlation genes
  high_cor_file <- paste0("high_cor_genes_", gsub(" ", "_", tissue), ".csv")
  
  if (!file.exists(high_cor_file) || grepl("_EMPTY\\.csv$", high_cor_file)) {
    cat("  No high correlation genes found\n")
    next
  }
  
  # Read high-correlation genes
  high_cor_data <- read.csv(high_cor_file, stringsAsFactors = FALSE)
  if (nrow(high_cor_data) == 0) {
    cat("  high_cor_data has 0 rows\n")
    next
  }
  
  high_cor_genes <- unique(high_cor_data$Gene)
  
  # Background genes (all genes in this GTEx dataset)
  data_obj <- read_gtex_data(file, apply_median_filter = TRUE)
  background_genes <- unique(data_obj$gene_names)
  
  # Restrict OCD genes to genes present in the background
  ocd_in_background <- intersect(ocd_genes, background_genes)
  
  # Compute contingency table
  a <- sum(high_cor_genes %in% ocd_in_background)                          # OCD ∩ high_cor
  b <- sum(!(high_cor_genes %in% ocd_in_background))                       # non-OCD ∩ high_cor
  c <- sum((background_genes %in% ocd_in_background) &
             !(background_genes %in% high_cor_genes))                      # OCD ∩ !high_cor
  d <- sum(!(background_genes %in% ocd_in_background) &
             !(background_genes %in% high_cor_genes))                      # non-OCD ∩ !high_cor
  
  cat("  Contingency table: a =", a, ", b =", b, ", c =", c, ", d =", d, "\n")
  
  # Numbers in background: OCD genes, total genes, and high_cor set size
  K_bg <- length(ocd_in_background)        # OCD genes in background
  N_bg <- length(background_genes)         # Total genes in background
  n_high <- length(high_cor_genes)        # Size of high_cor set
  
  # Fisher's exact test (one-sided, enrichment)
  fisher_p <- NA
  fisher_or <- NA
  enrichment_ratio <- NA
  
  if ((a + b) > 0 && K_bg > 0) {
    contingency_matrix <- matrix(c(a, b, c, d), nrow = 2)
    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
    fisher_p <- fisher_result$p.value
    fisher_or <- unname(fisher_result$estimate)
    
    enrichment_ratio <- (a / n_high) / (K_bg / N_bg)
  }
  
  ## 3) Empirical P-value based on permutation
  empirical_p <- NA
  perm_mean   <- NA
  perm_sd     <- NA
  perm_z      <- NA
  
  if (n_high > 0 && K_bg > 0) {
    # From background_genes, sample gene sets of the same size as high_cor n_perm times
    perm_overlaps <- replicate(n_perm, {
      random_genes <- sample(background_genes, size = n_high, replace = FALSE)
      sum(random_genes %in% ocd_in_background)
    })
    
    perm_mean <- mean(perm_overlaps)
    perm_sd   <- sd(perm_overlaps)
    perm_z    <- if (!is.na(perm_sd) && perm_sd > 0) (a - perm_mean) / perm_sd else NA
    
    # One-sided empirical P: Pr(overlap_perm >= observed a)
    empirical_p <- (sum(perm_overlaps >= a) + 1) / (n_perm + 1)
    
    cat("  Permutation: mean overlap =", perm_mean,
        ", sd =", perm_sd,
        ", empirical P =", empirical_p, "\n")
  } else {
    cat("  Skip permutation (no OCD genes in background or no high_cor genes)\n")
  }
  
  ## 4) Data for MAGMA gene-set (later written out as a file)
  if (n_high > 0) {
    magma_geneset_list[[tissue]] <- data.frame(
      SET  = gsub(" ", "_", tissue),
      GENE = high_cor_genes,
      stringsAsFactors = FALSE
    )
  }
  
  ## 5) Store results
  enrichment_results <- rbind(enrichment_results, data.frame(
    Tissue = tissue,
    OCD_genes_in_high_cor = a,
    Non_OCD_genes_in_high_cor = b,
    OCD_genes_non_in_high_cor = c,
    Non_OCD_genes_non_in_high_cor = d,
    Total_high_cor_genes = n_high,
    Total_OCD_genes_in_background = K_bg,
    Total_genes_in_background = N_bg,
    Enrichment_ratio = enrichment_ratio,
    P_value_Fisher = fisher_p,
    Odds_ratio_Fisher = fisher_or,
    Empirical_P_perm = empirical_p,
    Perm_mean_overlap = perm_mean,
    Perm_sd_overlap   = perm_sd,
    Perm_Z_score      = perm_z,
    stringsAsFactors = FALSE
  ))
}

## 6) FDR correction
enrichment_results$FDR_Fisher   <- p.adjust(enrichment_results$P_value_Fisher,
                                            method = "BH")
enrichment_results$FDR_Empirical <- p.adjust(enrichment_results$Empirical_P_perm,
                                             method = "BH")

enrichment_results$Significant_Fisher   <- enrichment_results$FDR_Fisher < 0.05
enrichment_results$Significant_Empirical <- enrichment_results$FDR_Empirical < 0.05

## 7) Sort by FDR
enrichment_results <- enrichment_results[order(enrichment_results$FDR_Empirical), ]

## 8) Save results
write.csv(enrichment_results,
          "OCD_gene_enrichment_results_Fisher_perm.csv",
          row.names = FALSE)



################################################################################



cat("\n=== OCD GENE ENRICHMENT ANALYSIS ===\n")
library(readxl)
# Read OCD genes
ocd_genes_pre <- read_excel("../gwas/41588_2025_2189_MOESM3_ESM.xlsx", skip=1, sheet=15)
ocd_genes <- ocd_genes_pre$Gene
head(ocd_genes)
cat("Loaded", length(ocd_genes), "OCD-associated genes\n") ## 251 genes

# Create directory for overlapping high correlation genes with OCD risk genes
dir.create("high_cor_OCD_overlap", showWarnings = FALSE)

# Perform enrichment analysis
enrichment_results <- data.frame()
all_overlapping_genes <- list()  # Store all overlapping genes by tissue

for (file in files) {
  tissue <- get_tissue_name(file)
  cat("\nAnalyzing", tissue, "\n")
  
  # Read corresponding high correlation genes file
  high_cor_file <- paste0("high_cor_genes_", gsub(" ", "_", tissue), ".csv")
  
  if (!file.exists(high_cor_file) || grepl("_EMPTY\\.csv$", high_cor_file)) {
    cat("  No high correlation genes found\n")
    next
  }
  
  # Read high correlation genes
  high_cor_data <- read.csv(high_cor_file, stringsAsFactors = FALSE)
  if (nrow(high_cor_data) == 0) next
  
  high_cor_genes <- unique(high_cor_data$Gene)
  
  # Read original data to get background genes
  data_obj <- read_gtex_data(file, apply_median_filter = TRUE)
  background_genes <- data_obj$gene_names
  
  # Identify high correlation genes that overlap with OCD risk genes
  high_cor_OCD_overlap <- high_cor_genes[high_cor_genes %in% ocd_genes]
  
  # Save overlapping high correlation genes with OCD risk genes to individual tissue file
  if (length(high_cor_OCD_overlap) > 0) {
    tissue_filename <- paste0("high_cor_OCD_overlap/high_cor_OCD_overlap_", gsub(" ", "_", tissue), ".txt")
    writeLines(high_cor_OCD_overlap, tissue_filename)
    cat("  Saved", length(high_cor_OCD_overlap), "high correlation genes overlapping with OCD risk genes to", tissue_filename, "\n")
    
    # Store in list for combined file
    all_overlapping_genes[[tissue]] <- high_cor_OCD_overlap
  }
  
  # Create contingency table
  a <- sum(high_cor_genes %in% ocd_genes)  # OCD genes in high cor
  b <- sum(!(high_cor_genes %in% ocd_genes))  # Non-OCD genes in high cor
  c <- sum((background_genes %in% ocd_genes) & !(background_genes %in% high_cor_genes))  # OCD not in high cor
  d <- sum(!(background_genes %in% ocd_genes) & !(background_genes %in% high_cor_genes))  # Non-OCD not in high cor
  
  cat("  Contingency table: a=", a, ", b=", b, ", c=", c, ", d=", d, "\n")
  
  # Fisher's exact test
  if (a + b > 0) {
    contingency_matrix <- matrix(c(a, b, c, d), nrow = 2)
    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
    
    enrichment_results <- rbind(enrichment_results, data.frame(
      Tissue = tissue,
      OCD_genes_in_high_cor = a,
      Non_OCD_genes_in_high_cor = b,
      OCD_genes_non_in_high_cor = c,
      Non_OCD_genes_non_in_high_cor = d,
      Total_high_cor_genes = length(high_cor_genes),
      Total_OCD_genes_in_background = sum(background_genes %in% ocd_genes),
      Total_genes_in_background = length(background_genes),
      Enrichment_ratio = (a / length(high_cor_genes)) / (sum(background_genes %in% ocd_genes) / length(background_genes)),
      P_value = fisher_result$p.value,
      Odds_ratio = fisher_result$estimate,
      stringsAsFactors = FALSE
    ))
  }
}

# Save combined file of high correlation genes overlapping with OCD risk genes
if (length(all_overlapping_genes) > 0) {
  combined_file <- "high_cor_OCD_overlap/all_high_cor_OCD_overlap_by_tissue.txt"
  sink(combined_file)
  cat("High correlation genes overlapping with OCD risk genes\n")
  cat("Total OCD risk genes:", length(ocd_genes), "\n\n")
  
  for (tissue in names(all_overlapping_genes)) {
    cat("===", tissue, "===\n")
    cat("Number of overlapping genes:", length(all_overlapping_genes[[tissue]]), "\n")
    cat(all_overlapping_genes[[tissue]], sep = "\n")
    cat("\n")
  }
  sink()
  cat("\nSaved combined high correlation genes overlapping with OCD risk genes to", combined_file, "\n")
}

# Apply FDR correction
enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
enrichment_results$Significant <- enrichment_results$FDR < 0.05
# Sort by FDR
enrichment_results <- enrichment_results[order(enrichment_results$FDR), ]
# Save results
write.csv(enrichment_results, "OCD_gene_enrichment_results_251genes.csv", row.names = FALSE)



                                            
