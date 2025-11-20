

library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)

# Get all file names
files <- list.files(pattern = "gene_tpm_v10_brain_.*\\.gct\\.gz$")

# Function to extract tissue name from filename
get_tissue_name <- function(filename) {
  tissue <- gsub("gene_tpm_v10_brain_", "", filename)
  tissue <- gsub(".gct.gz", "", tissue)
  tissue <- gsub("_", " ", tissue)
  # Clean up tissue names
  #tissue <- gsub("ba\\d+", "", tissue)  # Remove Brodmann area numbers
  #tissue <- gsub("\\.", " ", tissue)
  tissue <- tools::toTitleCase(tissue)
  return(tissue)
}


# Alternative approach: Put the panel label as part of the title
create_correlation_plot_v2 <- function(filename, gene1 = "KNCN", gene2 = "MKNK1-AS1", panel_label = NULL) {
  
  # Read the data
  data <- fread(filename, skip = 2, showProgress = FALSE)
  gene_ids <- data[[1]]
  gene_names <- data[[2]]
  expr_matrix <- as.matrix(data[, 3:ncol(data)])
  
  # Filter to genes with median TPM >= 0.1
  #median_tpm <- apply(expr_matrix, 1, median)
  #keep_idx <- which(median_tpm >= 0.1)
  #expr_matrix <- expr_matrix[keep_idx, ]
  #gene_names <- gene_names[keep_idx]
  #gene_ids <- gene_ids[keep_idx]
  
  # Find gene indices
  gene1_idx <- which(gene_names == gene1 | grepl(gene1, gene_ids))[1]
  gene2_idx <- which(gene_names == gene2 | grepl(gene2, gene_ids))[1]
  
  if (is.na(gene1_idx) || is.na(gene2_idx)) {
    return(NULL)
  }
  
  # Extract and log-transform expression
  gene1_expr <- log2(as.numeric(expr_matrix[gene1_idx, ]) + 1)
  gene2_expr <- log2(as.numeric(expr_matrix[gene2_idx, ]) + 1)
  
  # Create data frame
  plot_data <- data.frame(gene1 = gene1_expr, gene2 = gene2_expr)
  
  # Pearson correlation
  cor_test <- cor.test(gene1_expr, gene2_expr, method = "pearson")
  
  # Get tissue name
  tissue_name <- get_tissue_name(filename)
  plot_title <- if (!is.null(panel_label)) paste0("(", panel_label, ") ", tissue_name) else tissue_name
  
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
             label = sprintf("r = %.3f\np = %.3e", cor_test$estimate, cor_test$p.value),
             hjust = 1.1, vjust = -0.5, size = 3)
  
  return(p)
}


# Use version 2 (panel label in title) for cleaner appearance
plot_list <- list()
panel_labels <- letters[1:length(files)]

for(i in 1:length(files)) {
  cat("Processing", files[i], "\n")
  p <- create_correlation_plot_v2(files[i], panel_label = panel_labels[i])
  if(!is.null(p)) {
    plot_list[[i]] <- p
  }
}

# Continue with the rest of the code as before...

# Remove NULL plots
plot_list <- plot_list[!sapply(plot_list, is.null)]

# Arrange plots in a grid
n_plots <- length(plot_list)
n_cols <- 4
n_rows <- ceiling(n_plots / n_cols)

# Create the combined figure
combined_plot <- plot_grid(
  plotlist = plot_list,
  ncol = n_cols,
  nrow = n_rows,
  align = "hv",
  rel_heights = rep(1, n_rows),
  rel_widths = rep(1, n_cols)
)

# Add overall title
#title <- ggdraw() + 
#  draw_label(
#    "KNCN vs MKNK1-AS1 Expression Correlation Across Brain Regions",
#    fontface = 'bold',
#    size = 14
#  )

# Combine title and plots
final_plot <- plot_grid(
  #title, 
  combined_plot,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)


final_plot

# Save the figure
ggsave("KNCN_MKNK1AS1_correlations_all_tissues.pdf", 
       final_plot, 
       width = 12, 
       height = 10,
       dpi = 600)

png("KNCN_MKNK1AS1_correlations_all_tissues.png", 
    width = 12, 
    height = 10, 
    units = "in", 
    res = 600)
print(final_plot)
dev.off()



#########################################


correlation_results <- list()

for (i in 1:length(files)) {
  data <- fread(files[i], skip = 2, showProgress = FALSE)
  gene_ids <- data[[1]]
  gene_names <- data[[2]]
  expr_matrix <- as.matrix(data[, 3:ncol(data)])
  
  # Apply median TPM filter
  #median_tpm <- apply(expr_matrix, 1, median)
  #keep_idx <- which(median_tpm >= 0.1)
  #expr_matrix <- expr_matrix[keep_idx, ]
  #gene_names <- gene_names[keep_idx]
  #gene_ids <- gene_ids[keep_idx]
  
  gene1_idx <- which(gene_names == "KNCN" | grepl("KNCN", gene_ids))[1]
  gene2_idx <- which(gene_names == "MKNK1-AS1" | grepl("MKNK1-AS1", gene_ids))[1]
  
  if (!is.na(gene1_idx) && !is.na(gene2_idx)) {
    gene1_expr <- log2(as.numeric(expr_matrix[gene1_idx, ]) + 1)
    gene2_expr <- log2(as.numeric(expr_matrix[gene2_idx, ]) + 1)
    cor_test <- cor.test(gene1_expr, gene2_expr, method = "pearson")
    
    correlation_results[[i]] <- list(
      file = files[i],
      correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      tissue = get_tissue_name(files[i])
    )
  }
}

correlation_results <- correlation_results[!sapply(correlation_results, is.null)]


# Optional: Create a bar plot showing all correlations
# First, make sure cor_df includes p-values
cor_df <- data.frame(
  Tissue = sapply(correlation_results, function(x) x$tissue),
  Correlation = sapply(correlation_results, function(x) x$correlation),
  P_value = sapply(correlation_results, function(x) x$p_value),
  Significant = sapply(correlation_results, function(x) x$p_value < 0.05)
)

# Add significance stars
cor_df$Significance <- cut(cor_df$P_value, 
                           breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                           labels = c("***", "**", "*", ""))

# Create the bar plot with significance indicators
cor_bar_plot <- ggplot(cor_df, aes(x = reorder(Tissue, Correlation), 
                                   y = Correlation, 
                                   fill = Correlation > 0.7)) +
  geom_bar(stat = "identity") +
  # Add significance stars
  geom_text(aes(label = Significance, 
                y = Correlation + ifelse(Correlation > 0, 0.05, -0.05)),
            size = 5) +
  coord_flip() +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkred"),
                    labels = c("r ≤ 0.7", "r > 0.7")) +
  geom_hline(yintercept = 0.7, linetype = "dashed", color = "red") +
  labs(x = "Brain Region",
       y = "Pearson Correlation Coefficient", size = 13,
       fill = "Correlation Strength") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 14, face = "bold"),
        axis.text.x = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(size = 18, face = "bold"), 
        axis.title.y = element_text(size = 18, face = "bold"),
        # Legend
        legend.title = element_text(size = 13, face = "bold"),  # "Correlation Strength"
        legend.text = element_text(size = 13, face= "bold"),  # "r ≤ 0.7", "r > 0.7"
        legend.position = "bottom")
  # Add caption explaining significance
  #labs(caption = "Significance: * p < 0.05, ** p < 0.01, *** p < 0.001")


cor_bar_plot



png("KNCN_MKNK1AS1_correlation_barplot.png", 
    width = 8, 
    height = 6, 
    units = "in", 
    res = 600)
print(cor_bar_plot)
dev.off()

################################################################################


# Initialize results list
results_list <- list()

# Main loop over tissue files
for (file in files) {
  cat("Processing", file, "\n")
  
  # Read the data
  data <- fread(file, skip = 2, showProgress = FALSE)
  gene_ids <- data[[1]]
  gene_names <- data[[2]]
  expr_matrix <- as.matrix(data[, 3:ncol(data)])
  tissue <- get_tissue_name(file)
  
  # Find full indices BEFORE filtering
  kncn_idx_full <- which(gene_names == "KNCN" | grepl("KNCN", gene_ids))[1]
  mknk_idx_full <- which(gene_names == "MKNK1-AS1" | grepl("MKNK1-AS1", gene_ids))[1]
  
  if (is.na(kncn_idx_full) || is.na(mknk_idx_full)) {
    cat("Skipping", tissue, "- KNCN or MKNK1-AS1 not found\n")
    next
  }
  
  # Apply median TPM filter, but retain KNCN and MKNK1-AS1
  median_tpm <- apply(expr_matrix, 1, median)
  keep_idx <- which(median_tpm >= 0.1)
  required_idx <- c(kncn_idx_full, mknk_idx_full)
  final_idx <- union(keep_idx, required_idx)
  
  expr_matrix <- expr_matrix[final_idx, ]
  gene_names <- gene_names[final_idx]
  gene_ids <- gene_ids[final_idx]
  
  # Get new indices after filtering
  kncn_idx <- which(gene_names == "KNCN" | grepl("KNCN", gene_ids))[1]
  mknk_idx <- which(gene_names == "MKNK1-AS1" | grepl("MKNK1-AS1", gene_ids))[1]
  
  if (is.na(kncn_idx) || is.na(mknk_idx)) {
    cat("Skipping", tissue, "- KNCN or MKNK1-AS1 not retained after filtering\n")
    next
  }
  
  # Extract and log-transform expression vectors
  kncn_expr <- log2(as.numeric(expr_matrix[kncn_idx, ]) + 1)
  mknk_expr <- log2(as.numeric(expr_matrix[mknk_idx, ]) + 1)
  
  # Compute correlations with all genes
  gene_cor_table <- data.frame(Gene = gene_names, Cor_KNCN = NA, Cor_MKNK1AS1 = NA)
  
  for (i in 1:nrow(expr_matrix)) {
    gene_expr <- log2(as.numeric(expr_matrix[i, ]) + 1)
    gene_cor_table$Cor_KNCN[i] <- suppressWarnings(cor(kncn_expr, gene_expr, method = "pearson"))
    gene_cor_table$Cor_MKNK1AS1[i] <- suppressWarnings(cor(mknk_expr, gene_expr, method = "pearson"))
  }
  
  # Filter by correlation threshold
  filtered_genes <- gene_cor_table %>%
    filter(!is.na(Cor_KNCN) & !is.na(Cor_MKNK1AS1)) %>%
    filter(abs(Cor_KNCN) >= 0.7 & abs(Cor_MKNK1AS1) >= 0.7)
  
  # Save results
  output_filename <- paste0("high_cor_genes_", gsub(" ", "_", tissue), ".csv")
  
  if (nrow(filtered_genes) > 0) {
    filtered_genes$Tissue <- tissue
    write.csv(filtered_genes, file = output_filename, row.names = FALSE)
    results_list[[tissue]] <- filtered_genes
  } else {
    cat("No genes passed the threshold for:", tissue, "\n")
    write.csv(data.frame(Gene = character(0),
                         Cor_KNCN = numeric(0),
                         Cor_MKNK1AS1 = numeric(0),
                         Tissue = character(0)),
              file = paste0("high_cor_genes_", gsub(" ", "_", tissue), "_EMPTY.csv"),
              row.names = FALSE)
  }
}




caudate <- read.csv("high_cor_genes_Caudate_Basal_Ganglia.csv", header=T)
nucleus <- read.csv("high_cor_genes_Nucleus_Accumbens_Basal_Ganglia.csv", header=T)
putamen <- read.csv("high_cor_genes_Putamen_Basal_Ganglia.csv", header=T)
hypothal <- read.csv("high_cor_genes_Hypothalamus.csv", header=T)

cn <- caudate[which(caudate$Gene %in% nucleus$Gene),]
cnp <- cn[which(cn$Gene %in% putamen$Gene),]
write.table(cnp, "high_cor_genes_Caudate_Nucleus_Putamen_shared_1823genes.csv", col.names=T, row.names=F, quote=F, sep="\t")

cnph <- cnp[which(cnp$Gene %in% hypothal$Gene),]
write.table(cnph, "high_cor_genes_Caudate_Nucleus_Putamen_Hypothal_shared_868genes.csv", col.names=T, row.names=F, quote=F, sep="\t")



#################################################################################

library(data.table)
library(dplyr)
# OCD genes as a vector
ocd_genes <- read.table("OCD_gene_list.txt", stringsAsFactors = FALSE)[[1]]
# High correlation gene files
high_cor_files <- list.files(pattern = "^high_cor_genes_.*\\.csv$")
# Initialize results dataframe
enrichment_results <- data.frame()

# Initialize list to store OCD genes in high-correlated sets with tissue names
ocd_genes_in_high_cor_list <- list()

for (file in high_cor_files) {
  
  # Skip empty results
  if (grepl("_EMPTY.csv$", file)) next
  
  # Read high-correlation genes
  high_cor_genes <- read.csv(file)$Gene
  
  # Corresponding tissue-specific raw GTEx expression data
  raw_gct_file <- gsub("^high_cor_genes_", "gene_tpm_v10_brain_", gsub(".csv$", ".gct.gz", file))
  if (!file.exists(raw_gct_file)) next
  
  raw_data <- fread(raw_gct_file, skip = 2)
  
  # Tissue-specific expressed genes (e.g., mean TPM ≥ 1)
  expr_matrix <- as.matrix(raw_data[, 3:ncol(raw_data)])
  median_tpm <- apply(expr_matrix, 1, median)
  
  # Use a standard cutoff for tissue-specific background genes
  background_genes <- raw_data[[2]][median_tpm >= 0.1] # Adjust threshold as needed
  
  # Make contingency table values
  a <- sum(high_cor_genes %in% ocd_genes)
  b <- sum(!(high_cor_genes %in% ocd_genes))
  c <- sum((background_genes %in% ocd_genes) & !(background_genes %in% high_cor_genes))
  d <- sum((!(background_genes %in% ocd_genes)) & !(background_genes %in% high_cor_genes))
  
  # Get tissue name
  tissue_name <- gsub("high_cor_genes_|\\.csv", "", file)
  
  # Store OCD genes that are in the high-correlation set for this tissue
  ocd_genes_in_tissue <- high_cor_genes[high_cor_genes %in% ocd_genes]
  if (length(ocd_genes_in_tissue) > 0) {
    for (gene in ocd_genes_in_tissue) {
      ocd_genes_in_high_cor_list[[length(ocd_genes_in_high_cor_list) + 1]] <- 
        data.frame(Gene = gene, Tissue = tissue_name, stringsAsFactors = FALSE)
    }
  }
  
  # Fisher's exact test
  contingency_matrix <- matrix(c(a, b, c, d), nrow = 2)
  fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
  
  # Save results
  enrichment_results <- rbind(enrichment_results, data.frame(
    Tissue = tissue_name,
    OCD_genes_in_high_cor = a,
    Non_OCD_genes_in_high_cor = b,
    OCD_genes_non_in_high_cor = c,
    Non_OCD_genes_non_in_high_cor = d,
    Total_high_cor_genes = length(high_cor_genes),
    Total_OCD_genes = sum(background_genes %in% ocd_genes),
    Total_genes_in_tissue = length(background_genes),
    P_value = fisher_result$p.value,
    Odds_ratio = fisher_result$estimate,
    Enrichment_ratio = (a / length(high_cor_genes)) / (sum(background_genes %in% ocd_genes) / length(background_genes))
  ))
}

# Combine all OCD genes in high-correlation sets
if (length(ocd_genes_in_high_cor_list) > 0) {
  ocd_genes_in_high_cor_df <- do.call(rbind, ocd_genes_in_high_cor_list)
  
  # Save to text file
  write.table(ocd_genes_in_high_cor_df, 
              "OCD_genes_in_high_correlated_sets_by_tissue.txt", 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
  
  # Also save a summary version (unique genes per tissue)
  ocd_genes_summary <- ocd_genes_in_high_cor_df %>%
    group_by(Tissue) %>%
    summarise(OCD_genes = paste(sort(unique(Gene)), collapse = ", "),
              Count = n_distinct(Gene)) %>%
    arrange(desc(Count))
  
  write.table(ocd_genes_summary, 
              "OCD_genes_in_high_correlated_sets_summary.txt", 
              sep = "\t", 
              row.names = FALSE, 
              quote = FALSE)
}

# Adjust p-values for multiple tests
enrichment_results$FDR <- p.adjust(enrichment_results$P_value, method = "BH")
enrichment_results$Bonferroni <- p.adjust(enrichment_results$P_value, method = "bonferroni")
enrichment_results$Significant <- enrichment_results$FDR < 0.05
# Sort by FDR
enrichment_results <- enrichment_results[order(enrichment_results$FDR), ]
enrichment_results
# Save results
write.csv(enrichment_results, "OCD_gene_enrichment_tissue_specific.csv", row.names = FALSE)



# Extract OCD genes in high-correlated set (a)
genes_a <- high_cor_genes[high_cor_genes %in% ocd_genes]

# Create a small data frame with tissue label
if (length(genes_a) > 0) {
  gene_output <- data.frame(
    Tissue = tissue_name,
    OCD_Gene = genes_a
  )
  
  # Append to file (create if doesn't exist)
  write.table(gene_output,
              file = "OCD_genes_in_high_cor_by_tissue.txt",
              row.names = FALSE, col.names = !file.exists("OCD_genes_in_high_cor_by_tissue.txt"),
              sep = "\t", quote = FALSE, append = TRUE)
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

enrichment_plot

ggsave("OCD_gene_enrichment_barplot.png",
       enrichment_plot,
       width = 8,
       height = 6,
       dpi = 600)




# Apply Bonferroni correction with custom threshold
bonferroni_threshold <- 0.05/8  # 0.00625

# Add Bonferroni significance column
enrichment_results$Bonferroni_Significant <- enrichment_results$P_value < bonferroni_threshold

# Create enrichment visualization with Bonferroni correction
enrichment_plot <- ggplot(enrichment_results, 
                          aes(x = reorder(Tissue, -log10(P_value)), 
                              y = -log10(P_value),
                              fill = Bonferroni_Significant)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = -log10(bonferroni_threshold), 
             linetype = "dashed", 
             color = "red") +
  scale_fill_manual(values = c("FALSE" = "gray70", "TRUE" = "darkgreen"),
                    labels = c("NS", "p < 0.00625")) +
  coord_flip() +
  labs(x = "Brain Region",
       y = "-log10(p-value)",
       title = "OCD Gene Enrichment Across Brain Regions",
       subtitle = paste0("Bonferroni correction: p < 0.05/8 = ", 
                         format(bonferroni_threshold, scientific = FALSE)),
       fill = "Significance") +
  theme_classic() +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11, face = "italic"))

# Display the plot
enrichment_plot

# Save the plot
ggsave("OCD_gene_enrichment_barplot_bonferroni_8.png",
       enrichment_plot,
       width = 8,
       height = 6,
       dpi = 600)
