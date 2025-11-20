# Check Antisense lncRNA Overlaps with Canonical Transcripts Only

# =====================================================
# Required packages
# =====================================================
load_required_packages <- function() {
  required_packages <- c("GenomicRanges", "rtracklayer", "biomaRt", 
                         "AnnotationDbi", "dplyr", "stringr")
  
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      if (pkg %in% c("GenomicRanges", "rtracklayer", "AnnotationDbi")) {
        BiocManager::install(pkg)
      } else {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
  }
}

# =====================================================
# Method 1: Using GTF with Canonical Transcripts
# =====================================================

gtf_file = "gencode.v47.annotation.gtf.gz"
output_file = "genecode_v47_overlapping_result.txt"

find_antisense_overlaps_canonical_gtf <- function(gtf_file, 
                                                  output_file = "antisense_canonical_overlaps.txt",
                                                  species = "human") {
  
  load_required_packages()
  
  cat("Loading GTF file...\n")
  gtf <- rtracklayer::import(gtf_file)
  
  # Check if attributes were parsed correctly
  if (!("gene_type" %in% colnames(mcols(gtf)))) {
    cat("Note: GTF attributes need to be parsed from the attribute column\n")
  }
  
  gtf_df <- as.data.frame(gtf)
  
  # Check available columns
  cat("Available metadata columns:", paste(colnames(mcols(gtf)), collapse = ", "), "\n\n")
  
  # Get gene entries directly (full gene length)
  cat("Extracting gene boundaries (full length including introns)...\n")
  
  genes_df <- gtf_df[gtf_df$type == "gene", ]
  
  cat("Found", nrow(genes_df), "genes\n")
  
  # Separate protein-coding and antisense/lncRNA
  protein_coding <- genes_df[genes_df$gene_type == "protein_coding", ]
  
  # For antisense: include lncRNA type OR genes with -AS in the name
  antisense_lnc <- genes_df[
    genes_df$gene_type %in% c("lncRNA", "lincRNA") | 
    grepl("-AS[0-9]*$|-AS-", genes_df$gene_name),
  ]
  
  # Further filter to ensure we're not including protein-coding genes with -AS in name
  antisense_lnc <- antisense_lnc[antisense_lnc$gene_type != "protein_coding", ]
  
  cat("\nGene counts:\n")
  cat("  Protein-coding genes:", nrow(protein_coding), "\n")
  cat("  Antisense/lncRNA genes:", nrow(antisense_lnc), "\n")
  
  # Show how many have -AS in name
  n_with_as <- sum(grepl("-AS[0-9]*$|-AS-", antisense_lnc$gene_name))
  cat("  - With '-AS' in name:", n_with_as, "\n")
  cat("  - lncRNA type only:", nrow(antisense_lnc) - n_with_as, "\n")
  
  # Show gene type distribution
  cat("\nGene types in antisense/lncRNA category:\n")
  print(table(antisense_lnc$gene_type))
  
  # Convert to GRanges
  pc_granges <- makeGRangesFromDataFrame(
    protein_coding,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
  
  antisense_granges <- makeGRangesFromDataFrame(
    antisense_lnc,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand",
    keep.extra.columns = TRUE
  )
  
  # Find overlaps
  cat("\nFinding overlaps with opposite strand genes...\n")
  
  overlaps_list <- list()
  
  # Create progress bar
  pb <- txtProgressBar(min = 0, max = length(antisense_granges), style = 3)
  
  for (i in seq_len(length(antisense_granges))) {
    setTxtProgressBar(pb, i)
    antisense_gene <- antisense_granges[i]
    
    # Find ALL protein-coding genes on the same chromosome, regardless of strand
    same_chr_pc <- pc_granges[seqnames(pc_granges) == seqnames(antisense_gene)]
    
    # Find overlaps (strand-agnostic)
    overlaps <- findOverlaps(antisense_gene, same_chr_pc, type = "any", ignore.strand = TRUE)
    
    if (length(overlaps) > 0) {
      overlap_pc <- same_chr_pc[subjectHits(overlaps)]
      
      for (j in seq_along(overlap_pc)) {
        # Determine if it's same strand or opposite strand
        strand_relationship <- ifelse(
          as.character(strand(antisense_gene)) == as.character(strand(overlap_pc[j])),
          "same_strand",
          "opposite_strand"
        )
        
        overlap_info <- data.frame(
          antisense_name = antisense_gene$gene_name,
          antisense_id = antisense_gene$gene_id,
          antisense_type = antisense_gene$gene_type,
          antisense_chr = as.character(seqnames(antisense_gene)),
          antisense_start = start(antisense_gene),
          antisense_end = end(antisense_gene),
          antisense_strand = as.character(strand(antisense_gene)),
          antisense_length = width(antisense_gene),
          pc_name = overlap_pc[j]$gene_name,
          pc_id = overlap_pc[j]$gene_id,
          pc_chr = as.character(seqnames(overlap_pc[j])),
          pc_start = start(overlap_pc[j]),
          pc_end = end(overlap_pc[j]),
          pc_strand = as.character(strand(overlap_pc[j])),
          pc_length = width(overlap_pc[j]),
          strand_relationship = strand_relationship,
          overlap_start = max(start(antisense_gene), start(overlap_pc[j])),
          overlap_end = min(end(antisense_gene), end(overlap_pc[j])),
          overlap_length = min(end(antisense_gene), end(overlap_pc[j])) - 
            max(start(antisense_gene), start(overlap_pc[j])) + 1,
          stringsAsFactors = FALSE
        )
        
        overlaps_list[[length(overlaps_list) + 1]] <- overlap_info
      }
    }
  }
  close(pb)
  # Combine results
  all_overlaps <- do.call(rbind, overlaps_list)
  
  if (!is.null(all_overlaps) && nrow(all_overlaps) > 0) {
    cat("\nFound", nrow(all_overlaps), "antisense-protein coding overlaps\n")
    
    # Calculate overlap percentages
    all_overlaps$antisense_overlap_pct <- round(
      100 * all_overlaps$overlap_length / 
        (all_overlaps$antisense_end - all_overlaps$antisense_start + 1), 1
    )
    all_overlaps$pc_overlap_pct <- round(
      100 * all_overlaps$overlap_length / 
        (all_overlaps$pc_end - all_overlaps$pc_start + 1), 1
    )
    
    # Add combined name with overlap percentage
    all_overlaps$combined_name <- paste0(
      ifelse(is.na(all_overlaps$antisense_name) | all_overlaps$antisense_name == "", 
             all_overlaps$antisense_id, 
             all_overlaps$antisense_name),
      "-AS-",
      ifelse(is.na(all_overlaps$pc_name) | all_overlaps$pc_name == "", 
             all_overlaps$pc_id, 
             all_overlaps$pc_name)
    )
    
    # Add a version with percentage for display
    all_overlaps$combined_name_with_pct <- paste0(
      all_overlaps$combined_name,
      " (",
      pmax(all_overlaps$antisense_overlap_pct, all_overlaps$pc_overlap_pct),
      "%)"
    )
    
    # Print summary of overlaps
    cat("\nOverlap summary:\n")
    cat("  Mean overlap of antisense genes:", round(mean(all_overlaps$antisense_overlap_pct), 1), "%\n")
    cat("  Mean overlap of protein-coding genes:", round(mean(all_overlaps$pc_overlap_pct), 1), "%\n")
    cat("  Number with >50% overlap:", sum(pmax(all_overlaps$antisense_overlap_pct, 
                                                all_overlaps$pc_overlap_pct) > 50), "\n")
    cat("  Number with >90% overlap:", sum(pmax(all_overlaps$antisense_overlap_pct, 
                                                all_overlaps$pc_overlap_pct) > 90), "\n")
    
    # Sort by overlap percentage
    all_overlaps <- all_overlaps[order(-pmax(all_overlaps$antisense_overlap_pct, 
                                             all_overlaps$pc_overlap_pct)), ]
    
    # Show top overlaps
    cat("\nTop 10 overlaps by percentage:\n")
    top_overlaps <- head(all_overlaps, 10)
    print(data.frame(
      Antisense = top_overlaps$antisense_name,
      PC_Gene = top_overlaps$pc_name,
      AS_Overlap = paste0(top_overlaps$antisense_overlap_pct, "%"),
      PC_Overlap = paste0(top_overlaps$pc_overlap_pct, "%"),
      Length = top_overlaps$overlap_length
    ), row.names = FALSE)
    
    # Save results
    write.table(all_overlaps, output_file, sep = "\t", 
                row.names = FALSE, quote = FALSE)
    cat("\nResults saved to:", output_file, "\n")
  } else {
    cat("No overlaps found\n")
    all_overlaps <- data.frame()
  }
  
  return(all_overlaps)
}

write.table(all_overlaps, "genecode_v47_overlapping_result.txt", col.names=T, row.names=F, quote=F, sep="\t")

overlap_name <- fread("genecode_v47_overlapping_result.txt")



# Filter to keep best antisense lncRNA per protein-coding gene
# Priority: 1) antisense_name contains "-AS", 2) highest overlap percentage
library(dplyr)
library(data.table)

filter_best_antisense_per_pc <- function(overlaps_df = "genecode_v47_overlapping_result.txt", 
                                         output_file = "genecode_v47_overlapping_result_one_pair.txt") {
  
  # Read file if path is provided
  if (is.character(overlaps_df) && length(overlaps_df) == 1) {
    overlaps_df <- read.table(overlaps_df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  
  # Convert to data.frame to avoid conflicts
  overlaps_df <- as.data.frame(overlaps_df)
  
  cat("Filtering to keep best antisense lncRNA per protein-coding gene...\n")
  cat("Total overlaps before filtering:", nrow(overlaps_df), "\n")
  
  # Add a column to identify if antisense_name contains "-AS"
  overlaps_df$has_AS <- grepl("-AS", overlaps_df$antisense_name)
  
  # For each PC gene, apply filtering logic
  best_overlaps <- overlaps_df %>%
    group_by(pc_name) %>%
    mutate(
      # Count how many -AS entries exist for this PC gene
      n_AS = sum(has_AS)
    ) %>%
    filter(
      # If any -AS exists for this PC gene, keep only -AS entries
      if (any(has_AS)) has_AS else TRUE
    ) %>%
    # Among the remaining entries, keep the one with highest PC overlap percentage
    slice_max(pc_overlap_pct, n = 1, with_ties = FALSE) %>%
    # If there are still ties, keep the one with longer overlap
    slice_max(overlap_length, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    # Remove helper columns
    select(-has_AS, -n_AS)
  
  cat("\nResults:\n")
  cat("Unique PC genes:", n_distinct(best_overlaps$pc_name), "\n")
  cat("Selected antisense-PC pairs:", nrow(best_overlaps), "\n")
  cat("Pairs with -AS:", sum(grepl("-AS", best_overlaps$antisense_name)), "\n")
  cat("Removed overlaps:", nrow(overlaps_df) - nrow(best_overlaps), "\n")
  
  # Show which ones were removed
  removed_overlaps <- overlaps_df %>%
    anti_join(best_overlaps, by = c("antisense_name", "pc_name"))
  
  if (nrow(removed_overlaps) > 0) {
    cat("\nRemoved pairs:\n")
    removed_summary <- removed_overlaps[, c("pc_name", "antisense_name", "pc_overlap_pct")]
    removed_summary$has_AS <- grepl("-AS", removed_summary$antisense_name)
    removed_summary <- removed_summary[order(removed_summary$pc_name, -removed_summary$pc_overlap_pct), ]
    print(head(removed_summary, 20), row.names = FALSE)
  }
  
  # Sort by PC overlap percentage
  best_overlaps <- best_overlaps %>%
    arrange(-pc_overlap_pct)
  
  # Save results
  write.table(best_overlaps, output_file, sep = "\t", 
              row.names = FALSE, quote = FALSE)
  cat("\nFiltered results saved to:", output_file, "\n")
  
  return(best_overlaps)
}

# Alternative version without dplyr if you prefer
filter_best_antisense_per_pc_base <- function(overlaps_df, 
                                              output_file = "best_antisense_overlaps.txt") {
  
  # Read file if path is provided
  if (is.character(overlaps_df) && length(overlaps_df) == 1) {
    overlaps_df <- read.table(overlaps_df, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  }
  
  overlaps_df <- as.data.frame(overlaps_df)
  
  cat("Filtering to keep best antisense lncRNA per protein-coding gene...\n")
  cat("Total overlaps before filtering:", nrow(overlaps_df), "\n")
  
  # Add a column to identify if antisense_name contains "-AS"
  overlaps_df$has_AS <- grepl("-AS", overlaps_df$antisense_name)
  
  # Split by PC gene
  split_data <- split(overlaps_df, overlaps_df$pc_name)
  
  # For each PC gene, apply the filtering logic
  best_overlaps <- lapply(split_data, function(gene_data) {
    # Check if any antisense has "-AS" suffix
    has_AS_entries <- any(gene_data$has_AS)
    
    if (has_AS_entries) {
      # If there are -AS entries, filter to keep only those
      gene_data <- gene_data[gene_data$has_AS, ]
    }
    
    # Order by pc_overlap_pct (descending) and overlap_length (descending)
    gene_data <- gene_data[order(-gene_data$pc_overlap_pct, -gene_data$overlap_length), ]
    
    # Keep first row and remove helper column
    result <- gene_data[1, ]
    result$has_AS <- NULL
    return(result)
  })
  
  # Combine back
  best_overlaps <- do.call(rbind, best_overlaps)
  rownames(best_overlaps) <- NULL
  
  cat("\nResults:\n")
  cat("Unique PC genes:", length(unique(best_overlaps$pc_name)), "\n")
  cat("Selected antisense-PC pairs:", nrow(best_overlaps), "\n")
  cat("Pairs with -AS:", sum(grepl("-AS", best_overlaps$antisense_name)), "\n")
  cat("Removed overlaps:", nrow(overlaps_df) - nrow(best_overlaps), "\n")
  
  # Sort by PC overlap percentage
  best_overlaps <- best_overlaps[order(-best_overlaps$pc_overlap_pct), ]
  
  # Save results
  write.table(best_overlaps, output_file, sep = "\t", 
              row.names = FALSE, quote = FALSE)
  cat("\nFiltered results saved to:", output_file, "\n")
  
  return(best_overlaps)
}

# Example usage:
# result <- filter_best_antisense_per_pc("genecode_v47_overlapping_result.txt")
# or with a data frame:
# result <- filter_best_antisense_per_pc(my_overlaps_df)

best_overlaps <- filter_best_antisense_per_pc(overlap_name)



overlap_name <- fread("genecode_v47_overlapping_result_one_pair.txt")
overlap_name$merged_name <- paste0(overlap_name$pc_name, "/", overlap_name$antisense_name)
names(overlap_name)[9] <- c("gene_symbol")
head(overlap_name)
