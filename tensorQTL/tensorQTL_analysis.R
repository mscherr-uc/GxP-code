#!/usr/bin/env Rscript
# TensorQTL Results Analysis and Visualization
# Adapted for 8 treatment conditions: BPA_100nM_24, BPA_100nM_6, EtOH_24, EtOH_6, H2O_24, H2O_6, MBP_500nM_24, MBP_500nM_6

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(viridis)
library(corrplot)
library(VennDiagram)

# =============================================================================
# CONFIGURATION
# =============================================================================
base_dir <- "/rs/rs_grp_gxp/RNAseq_analysis/GxP_20250730/eQTL_mapping/tensor"
input_dir <- file.path(base_dir, "tensorqtl_output_cis")
output_dir <- file.path(base_dir, "analysis_results")
figures_dir <- file.path(output_dir, "figures")

# Create output directories
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Define treatment conditions
conditions <- c("BPA_100nM_24", "BPA_100nM_6", "EtOH_24", "EtOH_6", 
                "H2O_24", "H2O_6", "MBP_500nM_24", "MBP_500nM_6")

# FDR thresholds
fdr_strict <- 0.05
fdr_lenient <- 0.1

cat("Starting tensorQTL analysis for", length(conditions), "conditions\n")

# =============================================================================
# 1. LOAD AND PROCESS DATA
# =============================================================================
cat("Loading tensorQTL results...\n")

qtl_data <- list()
qtl_summary <- data.frame()

for (condition in conditions) {
  cat("Processing:", condition, "\n")
  
  # Load QTL results
  qtl_file <- file.path(input_dir, paste0(condition, "_cis.cis_qtl.txt.gz"))
  
  if (file.exists(qtl_file)) {
    qtl_results <- fread(qtl_file)
    qtl_results$condition <- condition
    qtl_data[[condition]] <- qtl_results
    
    # Calculate summary statistics
    n_tests <- nrow(qtl_results)
    n_genes <- length(unique(qtl_results$phenotype_id))
    n_variants <- length(unique(qtl_results$variant_id))
    
    # Count significant eQTLs (using pval_nominal since qvalue might not be available)
    # Note: Without multiple testing correction, we'll use a stricter p-value threshold
    sig_threshold <- 0.05 / n_tests  # Bonferroni correction as approximation
    n_sig_bonferroni <- sum(qtl_results$pval_nominal < sig_threshold, na.rm = TRUE)
    n_sig_pval001 <- sum(qtl_results$pval_nominal < 0.001, na.rm = TRUE)
    n_sig_pval0001 <- sum(qtl_results$pval_nominal < 0.0001, na.rm = TRUE)
    
    # Count eGenes (genes with at least one significant eQTL)
    egenes_bonferroni <- length(unique(qtl_results$phenotype_id[qtl_results$pval_nominal < sig_threshold]))
    egenes_001 <- length(unique(qtl_results$phenotype_id[qtl_results$pval_nominal < 0.001]))
    egenes_0001 <- length(unique(qtl_results$phenotype_id[qtl_results$pval_nominal < 0.0001]))
    
    summary_row <- data.frame(
      condition = condition,
      n_tests = n_tests,
      n_genes = n_genes,
      n_variants = n_variants,
      egenes_bonferroni = egenes_bonferroni,
      egenes_p001 = egenes_001,
      egenes_p0001 = egenes_0001,
      eqtls_bonferroni = n_sig_bonferroni,
      eqtls_p001 = n_sig_pval001,
      eqtls_p0001 = n_sig_pval0001,
      bonferroni_threshold = sig_threshold
    )
    
    qtl_summary <- rbind(qtl_summary, summary_row)
  } else {
    cat("WARNING: File not found:", qtl_file, "\n")
  }
}

# Combine all data
all_qtl_data <- do.call(rbind, qtl_data)

cat("Data loading complete. Total tests:", nrow(all_qtl_data), "\n")

# =============================================================================
# 2. SUMMARY STATISTICS TABLE
# =============================================================================
cat("Generating summary statistics...\n")

# Add treatment and timepoint columns for easier grouping
qtl_summary$treatment <- sapply(strsplit(qtl_summary$condition, "_"), function(x) x[1])
qtl_summary$timepoint <- sapply(strsplit(qtl_summary$condition, "_"), function(x) paste(x[2:length(x)], collapse = "_"))

# Write summary table
fwrite(qtl_summary, file.path(output_dir, "eQTL_summary_by_condition.txt"), 
       sep = '\t', quote = FALSE, row.names = FALSE)

# =============================================================================
# 3. VISUALIZATION: SUMMARY BARPLOTS
# =============================================================================
cat("Creating summary visualizations...\n")

# eGenes by condition
p1 <- ggplot(qtl_summary, aes(x = condition, y = egenes_p001)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = egenes_p001), vjust = -0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of eGenes by Condition (p < 0.001)",
       x = "Condition", y = "Number of eGenes") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figures_dir, "eGenes_by_condition.png"), p1, 
       width = 10, height = 6, dpi = 300)

# eQTLs by condition
p2 <- ggplot(qtl_summary, aes(x = condition, y = eqtls_p001)) +
  geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
  geom_text(aes(label = eqtls_p001), vjust = -0.3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Significant eQTLs by Condition (p < 0.001)",
       x = "Condition", y = "Number of eQTLs") +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(figures_dir, "eQTLs_by_condition.png"), p2, 
       width = 10, height = 6, dpi = 300)

# =============================================================================
# 4. QQ PLOTS FOR EACH CONDITION
# =============================================================================
cat("Creating QQ plots...\n")

create_qq_plot <- function(condition_data, condition_name, output_path) {
  # Prepare data for QQ plot
  pvals <- condition_data$pval_nominal[!is.na(condition_data$pval_nominal)]
  pvals <- pvals[pvals > 0]  # Remove zero p-values
  pvals <- sort(pvals)
  
  n <- length(pvals)
  expected <- (1:n) / n
  
  # Clip very small p-values for visualization
  pvals[pvals < 1e-20] <- 1e-20
  
  qq_data <- data.frame(
    expected = -log10(expected),
    observed = -log10(pvals)
  )
  
  # Calculate confidence intervals (optional)
  ci <- 0.95
  qq_data$clower <- -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n - (1:n) + 1))
  qq_data$cupper <- -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n - (1:n) + 1))
  
  p <- ggplot(qq_data, aes(x = expected, y = observed)) +
    geom_ribbon(aes(ymin = clower, ymax = cupper), alpha = 0.2, fill = "gray") +
    geom_point(alpha = 0.6, size = 0.8) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    theme_classic() +
    labs(title = paste("QQ Plot -", condition_name),
         x = expression(Expected -log[10](p)),
         y = expression(Observed -log[10](p))) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  ggsave(output_path, p, width = 8, height = 8, dpi = 300)
  
  return(p)
}

for (condition in conditions) {
  if (condition %in% names(qtl_data)) {
    output_file <- file.path(figures_dir, paste0("QQ_plot_", condition, ".png"))
    create_qq_plot(qtl_data[[condition]], condition, output_file)
  }
}

# =============================================================================
# 5. TREATMENT COMPARISON ANALYSIS
# =============================================================================
cat("Performing treatment comparisons...\n")

# Extract significant eQTLs for each condition (using p < 0.001)
sig_eqtls <- list()
for (condition in conditions) {
  if (condition %in% names(qtl_data)) {
    sig_data <- qtl_data[[condition]][qtl_data[[condition]]$pval_nominal < 0.001, ]
    sig_eqtls[[condition]] <- unique(sig_data$phenotype_id)
  }
}

# Create overlap matrix for eGenes
overlap_matrix <- matrix(0, nrow = length(conditions), ncol = length(conditions))
rownames(overlap_matrix) <- conditions
colnames(overlap_matrix) <- conditions

for (i in 1:length(conditions)) {
  for (j in 1:length(conditions)) {
    if (conditions[i] %in% names(sig_eqtls) && conditions[j] %in% names(sig_eqtls)) {
      overlap_matrix[i, j] <- length(intersect(sig_eqtls[[conditions[i]]], 
                                               sig_eqtls[[conditions[j]]]))
    }
  }
}

# Correlation heatmap of eGene overlaps
png(file.path(figures_dir, "eGene_overlap_heatmap.png"), 
    width = 12, height = 10, units = "in", res = 300)
corrplot(overlap_matrix, method = "color", type = "full", 
         order = "hclust", tl.cex = 0.8, tl.col = "black",
         addCoef.col = "black", number.cex = 0.7,
         title = "eGene Overlap Between Conditions", mar = c(0,0,1,0))
dev.off()

# =============================================================================
# 6. SAVE SIGNIFICANT RESULTS
# =============================================================================
cat("Saving significant results...\n")

for (condition in conditions) {
  if (condition %in% names(qtl_data)) {
    # Save top significant eQTLs
    sig_results <- qtl_data[[condition]][qtl_data[[condition]]$pval_nominal < 0.001, ]
    sig_results <- sig_results[order(sig_results$pval_nominal), ]
    
    output_file <- file.path(output_dir, paste0(condition, "_significant_eQTLs_p001.txt"))
    fwrite(sig_results, output_file, sep = '\t', quote = FALSE, row.names = FALSE)
    
    # Save just the gene-variant pairs for downstream analysis
    pairs <- sig_results[, .(phenotype_id, variant_id)]
    pairs_file <- file.path(output_dir, paste0(condition, "_significant_pairs_p001.txt"))
    fwrite(pairs, pairs_file, sep = '\t', quote = FALSE, row.names = FALSE)
  }
}

# =============================================================================
# 7. FINAL SUMMARY REPORT
# =============================================================================
cat("Generating final summary report...\n")

total_egenes <- sum(qtl_summary$egenes_p001)
total_eqtls <- sum(qtl_summary$eqtls_p001)
total_tests <- sum(qtl_summary$n_tests)

summary_report <- paste0(
  "TensorQTL Analysis Summary\n",
  "=" %&% paste(rep("=", 50), collapse = "") %&% "\n",
  "Total conditions analyzed: ", length(conditions), "\n",
  "Total statistical tests: ", format(total_tests, big.mark = ","), "\n",
  "Total eGenes (p < 0.001): ", total_egenes, "\n",
  "Total significant eQTLs (p < 0.001): ", format(total_eqtls, big.mark = ","), "\n\n",
  "Results saved to: ", output_dir, "\n",
  "Figures saved to: ", figures_dir, "\n\n",
  "Files generated:\n",
  "- eQTL_summary_by_condition.txt: Summary statistics table\n",
  "- *_significant_eQTLs_p001.txt: Significant eQTL results per condition\n",
  "- *_significant_pairs_p001.txt: Gene-variant pairs per condition\n",
  "- figures/: QQ plots and summary visualizations\n"
)

# Fix the string concatenation
summary_report <- paste0(
  "TensorQTL Analysis Summary\n",
  paste(rep("=", 50), collapse = ""), "\n",
  "Total conditions analyzed: ", length(conditions), "\n",
  "Total statistical tests: ", format(total_tests, big.mark = ","), "\n",
  "Total eGenes (p < 0.001): ", total_egenes, "\n",
  "Total significant eQTLs (p < 0.001): ", format(total_eqtls, big.mark = ","), "\n\n",
  "Results saved to: ", output_dir, "\n",
  "Figures saved to: ", figures_dir, "\n\n",
  "Files generated:\n",
  "- eQTL_summary_by_condition.txt: Summary statistics table\n",
  "- *_significant_eQTLs_p001.txt: Significant eQTL results per condition\n",
  "- *_significant_pairs_p001.txt: Gene-variant pairs per condition\n",
  "- figures/: QQ plots and summary visualizations\n"
)

writeLines(summary_report, file.path(output_dir, "ANALYSIS_SUMMARY.txt"))
cat(summary_report)

cat("Analysis complete! Check the output directory for results.\n")
