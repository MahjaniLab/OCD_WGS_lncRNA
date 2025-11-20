

library(dplyr)
library(stringr)
library(data.table)
library(SKAT)


# =====================================================
# STEP 1: Load and Prepare Your Data
# =====================================================

load_genotype_and_gene_data <- function(genotype_file = "input_variant_GT_code_transverse.txt",
                                       variant_tag_file = "variant_tag.txt") {
  
  cat("Loading genotype matrix...\n")
  # Read genotype matrix (variants as columns, samples as rows)
  genotype_data <- read.table(genotype_file, 
                             header = TRUE, 
                             row.names = 1, 
                             check.names = FALSE,
                             stringsAsFactors = FALSE)
  
  # Force conversion to numeric matrix
  genotype_matrix <- apply(genotype_data, 2, as.numeric)
  rownames(genotype_matrix) <- rownames(genotype_data)
  
  cat("Loaded genotype matrix:", nrow(genotype_matrix), "samples x", 
      ncol(genotype_matrix), "variants\n")
  
  # Check data type
  cat("Data type:", class(genotype_matrix[1,1]), "\n")
  
  # Check for any remaining issues
  if (!is.numeric(genotype_matrix[1,1])) {
    stop("Failed to convert genotype matrix to numeric. Please check your input file format.")
  }
  
  # Read variant-gene mapping
  cat("Loading variant-gene mapping...\n")
  variant_gene_map <- read.table(variant_tag_file, 
                                header = TRUE, 
                                stringsAsFactors = FALSE,
                                col.names = c("variant", "ENSG"))
  
  cat("Loaded", nrow(variant_gene_map), "variant-gene mappings\n")
  
  # Check which variants from mapping are in genotype matrix
  common_variants <- intersect(variant_gene_map$variant, colnames(genotype_matrix))
  cat("Found", length(common_variants), "variants in both files\n")
  
  # Filter to common variants
  genotype_matrix <- genotype_matrix[, common_variants, drop = FALSE]
  variant_gene_map <- variant_gene_map[variant_gene_map$variant %in% common_variants, ]
  
  return(list(
    genotype = genotype_matrix,
    gene_map = variant_gene_map
  ))
}

# =====================================================
# STEP 2: Prepare Phenotype and Covariates
# =====================================================

prepare_phenotype_covariates <- function(sample_ids, phenotype_file = NULL, 
                                       covariate_file = NULL) {
  
  n_samples <- length(sample_ids)
  
  if (!is.null(covariate_file) && file.exists(covariate_file)) {
    # Load combined phenotype and covariate file
    cat("Loading phenotype and covariates from:", covariate_file, "\n")
    
    # Read the file
    data <- read.table(covariate_file, header = TRUE, stringsAsFactors = FALSE)
    
    # Match to sample IDs
    # The ID column should match the sample IDs from genotype matrix
    matched_idx <- match(sample_ids, data$ID)
    
    if (sum(!is.na(matched_idx)) == 0) {
      stop("No matching sample IDs found between genotype and covariate files!")
    }
    
    # Extract phenotype (assuming column name is 'pheno')
    if ("pheno" %in% colnames(data)) {
      phenotype <- data$pheno[matched_idx]
      names(phenotype) <- sample_ids
      cat("Extracted phenotype from 'pheno' column\n")
    } else {
      stop("No 'pheno' column found in covariate file!")
    }
    
    # Extract covariates (only X1-X5 for PC1-5)
    pc_cols <- paste0("X", 1:5)  # X1, X2, X3, X4, X5
    available_pcs <- pc_cols[pc_cols %in% colnames(data)]
    
    if (length(available_pcs) > 0) {
      covariates <- data[matched_idx, available_pcs, drop = FALSE]
      # Rename columns to be more descriptive
      colnames(covariates) <- paste0("PC", 1:length(available_pcs))
      rownames(covariates) <- sample_ids
      cat("Extracted", ncol(covariates), "principal components:", paste(colnames(covariates), collapse = ", "), "\n")
    } else {
      # If no X columns, use all numeric columns except ID, pheno, matches
      exclude_cols <- c("ID", "pheno", "matches")
      covariate_cols <- setdiff(colnames(data), exclude_cols)
      numeric_cols <- covariate_cols[sapply(data[covariate_cols], is.numeric)]
      
      if (length(numeric_cols) > 0) {
        covariates <- data[matched_idx, numeric_cols, drop = FALSE]
        rownames(covariates) <- sample_ids
        cat("Extracted", ncol(covariates), "covariates:", paste(numeric_cols, collapse = ", "), "\n")
      } else {
        cat("Warning: No covariate columns found. Using intercept-only model.\n")
        covariates <- data.frame(intercept = rep(1, n_samples), row.names = sample_ids)
      }
    }
    
    # Report matching statistics
    cat("Matched", sum(!is.na(matched_idx)), "out of", length(sample_ids), "samples\n")
    
  } else if (!is.null(phenotype_file) && file.exists(phenotype_file)) {
    # Load phenotype from separate file
    pheno_data <- read.table(phenotype_file, header = TRUE, stringsAsFactors = FALSE)
    # Match to sample IDs
    phenotype <- pheno_data$phenotype[match(sample_ids, pheno_data$sample_id)]
    
    # Generate example covariates
    cat("No covariate file provided. Generating example covariates.\n")
    covariates <- data.frame(
      PC1 = rnorm(n_samples),
      PC2 = rnorm(n_samples),
      PC3 = rnorm(n_samples),
      row.names = sample_ids
    )
  } else {
    # Generate example phenotype and covariates
    cat("No phenotype/covariate file provided. Generating example data.\n")
    phenotype <- sample(0:1, n_samples, replace = TRUE)
    names(phenotype) <- sample_ids
    
    covariates <- data.frame(
      PC1 = rnorm(n_samples),
      PC2 = rnorm(n_samples),
      PC3 = rnorm(n_samples),
      row.names = sample_ids
    )
  }
  
  # Remove samples with missing phenotype
  missing_pheno <- is.na(phenotype)
  if (sum(missing_pheno) > 0) {
    cat("Removing", sum(missing_pheno), "samples with missing phenotype\n")
    phenotype <- phenotype[!missing_pheno]
    covariates <- covariates[!missing_pheno, , drop = FALSE]
  }
  
  return(list(phenotype = phenotype, covariates = covariates))
}

# =====================================================
# STEP 3: Main SKAT-O Analysis Function
# =====================================================

run_SKAT_O_analysis <- function(genotype_file = "input_variant_GT_code_transverse.txt",
                               variant_tag_file = "variant_tag.txt",
                               phenotype_file = NULL,
                               covariate_file = NULL,
                               output_prefix = "SKAT_O_results",
                               maf_threshold = 0.01) {
  
  # Load data
  cat("=== Loading Data ===\n")
  data <- load_genotype_and_gene_data(genotype_file, variant_tag_file)
  genotype_matrix <- data$genotype
  variant_gene_map <- data$gene_map
  
  # Get sample IDs
  sample_ids <- rownames(genotype_matrix)
  
  # Prepare phenotype and covariates
  cat("\n=== Preparing Phenotype and Covariates ===\n")
  pheno_covar <- prepare_phenotype_covariates(sample_ids, phenotype_file, covariate_file)
  phenotype <- pheno_covar$phenotype
  covariates <- pheno_covar$covariates
  
  # Create null model
  cat("\n=== Creating Null Model ===\n")
  if (length(unique(phenotype)) == 2) {
    cat("Binary phenotype detected. Using logistic model.\n")
    obj <- SKAT_Null_Model(phenotype ~ ., data = covariates, out_type = "D")
  } else {
    cat("Continuous phenotype detected. Using linear model.\n")
    obj <- SKAT_Null_Model(phenotype ~ ., data = covariates, out_type = "C")
  }
  
  # Get unique genes
  unique_genes <- unique(variant_gene_map$ENSG)
  cat("\n=== Analyzing", length(unique_genes), "genes ===\n")
  
  # Initialize results
  results_list <- list()
  
  # Progress counter
  pb <- txtProgressBar(min = 0, max = length(unique_genes), style = 3)
  
  # Analyze each gene
  for (i in seq_along(unique_genes)) {
    gene <- unique_genes[i]
    setTxtProgressBar(pb, i)
    
    # Get variants for this gene
    gene_variants <- variant_gene_map$variant[variant_gene_map$ENSG == gene]
    
    if (length(gene_variants) >= 2) {  # Need at least 2 variants for SKAT
      # Extract genotype matrix for this gene
      gene_genotype <- genotype_matrix[, gene_variants, drop = FALSE]
      
      # Calculate MAF for each variant
      maf <- apply(gene_genotype, 2, function(x) {
        # Handle missing values
        x <- x[!is.na(x)]
        if (length(x) == 0) return(NA)
        af <- sum(x) / (2 * length(x))
        return(min(af, 1 - af))
      })
      
      # Filter for rare variants
      rare_mask <- !is.na(maf) & maf < maf_threshold & maf > 0
      n_rare <- sum(rare_mask)
      
      if (n_rare >= 2) {
        gene_genotype_rare <- gene_genotype[, rare_mask, drop = FALSE]
        
        tryCatch({
          # Run SKAT-O
          skat_o_result <- SKAT(gene_genotype_rare, obj, method = "optimal.adj")
          
          # Also run individual tests
          skat_result <- SKAT(gene_genotype_rare, obj, method = "SKAT")
          burden_result <- SKAT(gene_genotype_rare, obj, method = "Burden")
          
          # Store results
          results_list[[gene]] <- data.frame(
            gene = gene,
            n_variants_total = length(gene_variants),
            n_variants_tested = n_rare,
            mean_maf = mean(maf[rare_mask], na.rm = TRUE),
            min_maf = min(maf[rare_mask], na.rm = TRUE),
            max_maf = max(maf[rare_mask], na.rm = TRUE),
            p_SKAT = skat_result$p.value,
            p_Burden = burden_result$p.value,
            p_SKATO = skat_o_result$p.value,
            optimal_rho = ifelse(!is.null(skat_o_result$param$rho), 
                               skat_o_result$param$rho, NA),
            n_samples_tested = nrow(gene_genotype_rare),
            stringsAsFactors = FALSE
          )
        }, error = function(e) {
          results_list[[gene]] <- data.frame(
            gene = gene,
            n_variants_total = length(gene_variants),
            n_variants_tested = n_rare,
            mean_maf = mean(maf[rare_mask], na.rm = TRUE),
            min_maf = min(maf[rare_mask], na.rm = TRUE),
            max_maf = max(maf[rare_mask], na.rm = TRUE),
            p_SKAT = NA,
            p_Burden = NA,
            p_SKATO = NA,
            optimal_rho = NA,
            n_samples_tested = nrow(gene_genotype_rare),
            error = as.character(e),
            stringsAsFactors = FALSE
          )
        })
      } else {
        # Not enough rare variants
        results_list[[gene]] <- data.frame(
          gene = gene,
          n_variants_total = length(gene_variants),
          n_variants_tested = n_rare,
          mean_maf = ifelse(n_rare > 0, mean(maf[rare_mask], na.rm = TRUE), NA),
          min_maf = ifelse(n_rare > 0, min(maf[rare_mask], na.rm = TRUE), NA),
          max_maf = ifelse(n_rare > 0, max(maf[rare_mask], na.rm = TRUE), NA),
          p_SKAT = NA,
          p_Burden = NA,
          p_SKATO = NA,
          optimal_rho = NA,
          n_samples_tested = nrow(gene_genotype),
          error = "Less than 2 rare variants",
          stringsAsFactors = FALSE
        )
      }
    } else {
      # Not enough variants total
      results_list[[gene]] <- data.frame(
        gene = gene,
        n_variants_total = length(gene_variants),
        n_variants_tested = 0,
        mean_maf = NA,
        min_maf = NA,
        max_maf = NA,
        p_SKAT = NA,
        p_Burden = NA,
        p_SKATO = NA,
        optimal_rho = NA,
        n_samples_tested = nrow(genotype_matrix),
        error = "Less than 2 variants in gene",
        stringsAsFactors = FALSE
      )
    }
  }
  close(pb)
  
  # Combine results
  cat("\n\n=== Combining Results ===\n")
  results_df <- do.call(rbind, results_list)
  
  # Add multiple testing correction
  results_df$p_SKATO_FDR <- NA
  results_df$p_SKATO_Bonferroni <- NA
  
  # Only adjust p-values that are not NA
  valid_p <- !is.na(results_df$p_SKATO)
  if (sum(valid_p) > 0) {
    results_df$p_SKATO_FDR[valid_p] <- p.adjust(results_df$p_SKATO[valid_p], method = "fdr")
    results_df$p_SKATO_Bonferroni[valid_p] <- p.adjust(results_df$p_SKATO[valid_p], method = "bonferroni")
  }
  
  # Sort by p-value
  results_df <- results_df[order(results_df$p_SKATO), ]
  
  # Save results
  output_file <- paste0(output_prefix, ".txt")
  write.table(results_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("Results saved to:", output_file, "\n")
  
  # Print summary
  cat("\n=== Summary Statistics ===\n")
  cat("Total genes analyzed:", nrow(results_df), "\n")
  cat("Genes with sufficient variants:", sum(valid_p), "\n")
  cat("Genes with p < 0.05:", sum(results_df$p_SKATO < 0.05, na.rm = TRUE), "\n")
  cat("Genes with FDR < 0.05:", sum(results_df$p_SKATO_FDR < 0.05, na.rm = TRUE), "\n")
  cat("Genes with Bonferroni < 0.05:", sum(results_df$p_SKATO_Bonferroni < 0.05, na.rm = TRUE), "\n")
  
  return(results_df)
}

# =====================================================
# STEP 4: Visualization Functions
# =====================================================

plot_SKAT_results <- function(results_df, output_prefix = "SKAT_O_plots") {
  library(ggplot2)
  library(ggrepel)
  
  # Filter to genes with valid p-values
  results_plot <- results_df[!is.na(results_df$p_SKATO), ]
  
  if (nrow(results_plot) == 0) {
    cat("No valid p-values to plot.\n")
    return(NULL)
  }
  
  # Manhattan-like plot
  results_plot$log10_p <- -log10(results_plot$p_SKATO)
  results_plot$rank <- 1:nrow(results_plot)
  
  # Significance thresholds
  sig_nominal <- -log10(0.05)
  sig_bonferroni <- -log10(0.05 / nrow(results_plot))
  
  p1 <- ggplot(results_plot, aes(x = rank, y = log10_p)) +
    geom_point(aes(size = n_variants_tested, 
                   color = log10_p > sig_bonferroni),
               alpha = 0.6) +
    geom_hline(yintercept = sig_nominal, linetype = "dashed", 
               color = "blue", alpha = 0.5) +
    geom_hline(yintercept = sig_bonferroni, linetype = "dashed", 
               color = "red", alpha = 0.5) +
    scale_color_manual(values = c("gray30", "red"), guide = FALSE) +
    scale_size_continuous(range = c(2, 6)) +
    labs(x = "Gene Rank", 
         y = expression(-log[10](p-value)),
         title = "SKAT-O Gene-based Association Results",
         size = "# Variants") +
    theme_minimal() +
    theme(panel.grid.minor = element_blank())
  
  # Add labels for top genes
  top_genes <- head(results_plot, 10)
  if (nrow(top_genes) > 0 && max(top_genes$log10_p) > 2) {
    p1 <- p1 + geom_text_repel(data = top_genes,
                               aes(label = gene),
                               size = 3,
                               box.padding = 0.5)
  }
  
  # QQ plot
  n <- nrow(results_plot)
  expected <- -log10(seq(1/n, 1 - 1/n, length.out = n))
  observed <- sort(results_plot$log10_p, decreasing = TRUE)
  
  qq_data <- data.frame(expected = expected, observed = observed)
  
  p2 <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_point(alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, color = "red", 
                linetype = "dashed") +
    labs(x = expression(Expected ~ -log[10](p)),
         y = expression(Observed ~ -log[10](p)),
         title = "QQ Plot for SKAT-O Results") +
    theme_minimal() +
    coord_fixed()
  
  # P-value distribution
  p3 <- ggplot(results_plot, aes(x = p_SKATO)) +
    geom_histogram(bins = 20, fill = "skyblue", color = "black", 
                   alpha = 0.7) +
    labs(x = "P-value", y = "Count",
         title = "P-value Distribution") +
    theme_minimal()
  
  # Save plots
  library(gridExtra)
  combined_plot <- grid.arrange(p1, p2, p3, 
                               layout_matrix = rbind(c(1, 1),
                                                    c(2, 3)))
  
  ggsave(paste0(output_prefix, ".pdf"), combined_plot, 
         width = 12, height = 10)
  
  cat("Plots saved to:", paste0(output_prefix, ".pdf"), "\n")
  
  return(list(manhattan = p1, qq = p2, hist = p3))
}



# Run SKAT-O analysis with your files
results <- run_SKAT_O_analysis(
  genotype_file = "input_variant_GT_code_transverse.txt",
  variant_tag_file = "variant_tag.txt",
  covariate_file = "ocd_DATA_PCA_matched_result_v8.txt",
  output_prefix = "SKAT_O_results_OCD"
)

# Generate visualization plots
plots <- plot_SKAT_results(results, "SKAT_O_plots_OCD")


