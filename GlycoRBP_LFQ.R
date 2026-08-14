## ============================================================================
## GlycoRBP MS Pipeline: rGAL vs UL pulldown analysis as the example
## Sections:
##   1. Setup & data import
##   2. Keratin contaminant removal 
##   3. Protein ID assignment
##   4. LFQ extraction, log2 transform, detection filtering
##   5. Group-specific / condition-exclusive protein definitions
##   6. MNAR / MAR classification
##   7. Imputation (KNN for MAR, QRILC for MNAR)
##   8. Imputation diagnostic plot
##   9. PCA (full imputed matrix)
##   10. Differential expression (limma, excluding strict condition-specific)
##   11. Volcano plot
##   12. GSEA (limma on full imputed matrix, rank by moderated t)
## ============================================================================

## ---- 1. Setup & data import ------------------------------------------------

library(openxlsx)
library(tidyverse)
library(readr)
library(limma)
library(dplyr)
library(VennDiagram)
library(grid)
library(ggplot2)
library(ggrepel)
library(msigdbr)
library(fgsea)
library(impute)
library(imputeLCMD)
library(clusterProfiler)
library(org.Hs.eg.db)

protein_df <- read_tsv(
  "LFQ_rGAL/combined_protein.tsv",
  show_col_types = FALSE
)

treatment_cols <- c(
  "rGAL_1 MaxLFQ Intensity",
  "rGAL_2 MaxLFQ Intensity",
  "rGAL_3 MaxLFQ Intensity"
)

control_cols <- c(
  "UL_1 MaxLFQ Intensity",
  "UL_2 MaxLFQ Intensity",
  "UL_3 MaxLFQ Intensity"
)

treatment_name <- "rGAL"
control_name   <- "UL"

sample_cols <- c(treatment_cols, control_cols)

group <- factor(
  c(rep(treatment_name, length(treatment_cols)),
    rep(control_name, length(control_cols))),
  levels = c(control_name, treatment_name)
)

## ---- 2. Keratin contaminant removal ---------------------------------------
## Keratins (KRT#) and keratin-associated proteins (KRTAP) are common
## environmental/skin/hair contaminants in MS sample prep and are removed
## before any downstream processing.

keratin_gene_pattern <- "^KRT[0-9]|^KRTAP"

is_keratin <- grepl(
  keratin_gene_pattern,
  protein_df$Gene,
  ignore.case = TRUE
)

# Backup check against protein description field, if present, in case
# a keratin entry lacks a standard KRT/KRTAP gene symbol
if ("Protein Description" %in% colnames(protein_df)) {
  is_keratin <- is_keratin | grepl(
    "keratin",
    protein_df$`Protein Description`,
    ignore.case = TRUE
  )
}

cat("Proteins before keratin removal:", nrow(protein_df), "\n")
cat("Keratin / keratin-associated proteins removed:", sum(is_keratin, na.rm = TRUE), "\n")

removed_keratins <- protein_df[which(is_keratin), , drop = FALSE]

protein_df <- protein_df[!is_keratin, , drop = FALSE]

cat("Proteins after keratin removal:", nrow(protein_df), "\n")

## ---- 3. Protein ID assignment ---------------------------------------------

if ("Gene" %in% colnames(protein_df)) {
  protein_id <- protein_df$Gene
} else if ("Gene Name" %in% colnames(protein_df)) {
  protein_id <- protein_df$`Gene Name`
} else if ("Protein" %in% colnames(protein_df)) {
  protein_id <- protein_df$Protein
} else {
  stop("Specify the correct protein or gene identifier column.")
}

missing_id <- is.na(protein_id) | trimws(protein_id) == ""

protein_id[missing_id] <- paste0("Unannotated_", which(missing_id))

protein_id <- make.unique(protein_id)

## ---- 4. LFQ extraction, log2 transform, detection filtering ---------------

lfq_raw <- protein_df[, sample_cols, drop = FALSE]

lfq_raw[] <- lapply(lfq_raw, function(x) as.numeric(as.character(x)))

lfq_raw <- as.matrix(lfq_raw)

rownames(lfq_raw) <- protein_id
colnames(lfq_raw) <- sample_cols

# Replace zero / non-finite with NA
lfq_raw[!is.finite(lfq_raw)] <- NA
lfq_raw[lfq_raw <= 0] <- NA

# Log2 transform
lfq_log2 <- log2(lfq_raw)

# Detection counts per group
treatment_idx <- seq_along(treatment_cols)
control_idx   <- length(treatment_cols) + seq_along(control_cols)

n_treatment <- rowSums(!is.na(lfq_log2[, treatment_idx, drop = FALSE]))
n_control   <- rowSums(!is.na(lfq_log2[, control_idx, drop = FALSE]))

# 2-of-3 filter: keep protein if detected in >=2/3 replicates in at least one group
keep <- n_treatment >= 2 | n_control >= 2

lfq_filtered <- lfq_log2[keep, , drop = FALSE]

n_treatment_filtered <- n_treatment[keep]
n_control_filtered   <- n_control[keep]

cat("Proteins before detection filtering:", nrow(lfq_log2), "\n")
cat("Proteins after detection filtering:", nrow(lfq_filtered), "\n")

detection_table <- data.frame(
  Protein = rownames(lfq_filtered),
  Treatment_detected = n_treatment_filtered,
  Control_detected = n_control_filtered,
  stringsAsFactors = FALSE
)

colnames(detection_table)[2:3] <- c(
  paste0(treatment_name, "_detected"),
  paste0(control_name, "_detected")
)

write.csv(
  lfq_filtered,
  paste0(treatment_name,"_filtered.csv"),
  row.names = TRUE
)

write.csv(
  detection_table,
  paste0(treatment_name, "_protein_detection_counts_before_imputation.csv"),
  row.names = FALSE
)

## ---- 5. Group-specific / condition-exclusive protein definitions ---------

# Venn: proteins detected (>=2/3) in each group
treatment_detected <- n_treatment >= 2
control_detected   <- n_control >= 2

treatment_proteins <- unique(rownames(lfq_log2)[treatment_detected])
control_proteins   <- unique(rownames(lfq_log2)[control_detected])

treatment_proteins <- treatment_proteins[!is.na(treatment_proteins) & treatment_proteins != ""]
control_proteins   <- control_proteins[!is.na(control_proteins) & control_proteins != ""]

venn.plot <- venn.diagram(
  x = list(UL = control_proteins, rGAL = treatment_proteins),
  filename = NULL,
  fill = c("skyblue", "salmon"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)

ggsave(paste0(treatment_name,"_venn.png"), plot = venn.plot, width = 6, height = 6, dpi = 300)

write.csv(
  data.frame(Protein = treatment_proteins),
  paste0(treatment_name, "_detected_2of3.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Protein = control_proteins),
  paste0(control_name, "_detected_2of3.csv"),
  row.names = FALSE
)

# Strict group-specific lists (>=2/3 in one group & 0/3 in the other)
treatment_specific_idx <- n_treatment >= 2 & n_control == 0
control_specific_idx   <- n_control >= 2 & n_treatment == 0

treatment_specific <- rownames(lfq_log2)[treatment_specific_idx]
control_specific   <- rownames(lfq_log2)[control_specific_idx]

write.csv(
  data.frame(Protein = treatment_specific),
  paste0(treatment_name, "_specific_2of3_vs_0of3.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(Protein = control_specific),
  paste0(control_name, "_specific_2of3_vs_0of3.csv"),
  row.names = FALSE
)

## ---- 6. MNAR / MAR classification -----------------------------------------
## MNAR: protein detected >=2/3 in one group & <=1/3 in the other group.
## MAR: all remaining missing values.

mnar_mask <- matrix(
  FALSE,
  nrow = nrow(lfq_filtered),
  ncol = ncol(lfq_filtered),
  dimnames = dimnames(lfq_filtered)
)

# Detected in control, missing in treatment
treatment_mnar_rows <- n_treatment_filtered <= 1 & n_control_filtered >= 2
mnar_mask[treatment_mnar_rows, treatment_idx] <- is.na(
  lfq_filtered[treatment_mnar_rows, treatment_idx, drop = FALSE]
)

# Detected in treatment, missing in control
control_mnar_rows <- n_control_filtered <= 1 & n_treatment_filtered >= 2
mnar_mask[control_mnar_rows, control_idx] <- is.na(
  lfq_filtered[control_mnar_rows, control_idx, drop = FALSE]
)

mar_mask <- is.na(lfq_filtered) & !mnar_mask

cat("Likely MAR entries:", sum(mar_mask), "\n")
cat("Likely MNAR entries:", sum(mnar_mask), "\n")

## ---- 7. Imputation (KNN for MAR, QRILC for MNAR) --------------------------

set.seed(12345)
knn_candidate <- impute::impute.knn(
  lfq_filtered,
  k = 3,
  rowmax = 0.9,
  colmax = 0.9,
  maxp = 1500,
  rng.seed = 12345
)$data

set.seed(12345)
qrilc_result <- imputeLCMD::impute.QRILC(lfq_filtered, tune.sigma = 1)
qrilc_candidate <- qrilc_result[[1]]
rownames(qrilc_candidate) <- rownames(lfq_filtered)
colnames(qrilc_candidate) <- colnames(lfq_filtered)

cat(
  "Missing values remaining in QRILC candidate matrix:",
  sum(is.na(qrilc_candidate)), "\n"
)

lfq_imputed <- lfq_filtered
lfq_imputed[mar_mask]  <- knn_candidate[mar_mask]
lfq_imputed[mnar_mask] <- qrilc_candidate[mnar_mask]

cat("Missing values after mixed imputation:", sum(is.na(lfq_imputed)), "\n")

imputation_summary <- data.frame(
  Protein = rownames(lfq_filtered),
  Treatment_detected = n_treatment_filtered,
  Control_detected = n_control_filtered,
  MAR_values_imputed = rowSums(mar_mask),
  MNAR_values_imputed = rowSums(mnar_mask),
  stringsAsFactors = FALSE
)

colnames(imputation_summary)[2:3] <- c(
  paste0(treatment_name, "_detected"),
  paste0(control_name, "_detected")
)

write.csv(
  imputation_summary,
  paste0(treatment_name, "_vs_", control_name, "_missingness_summary.csv"),
  row.names = FALSE
)

write.csv(
  lfq_imputed,
  paste0(treatment_name, "_vs_", control_name, "_ALL_imputed_log2_LFQ.csv")
)

## ---- 8. Imputation diagnostic plot ----------------------------------------

observed_values <- lfq_filtered[!is.na(lfq_filtered)]
mar_values  <- lfq_imputed[mar_mask]
mnar_values <- lfq_imputed[mnar_mask]

imputation_plot_df <- bind_rows(
  data.frame(Intensity = as.numeric(observed_values), Value_type = "Observed LFQ values"),
  data.frame(Intensity = as.numeric(mar_values),  Value_type = "KNN-imputed MAR values"),
  data.frame(Intensity = as.numeric(mnar_values), Value_type = "QRILC-imputed MNAR values")
) %>%
  filter(is.finite(Intensity)) %>%
  mutate(
    Value_type = factor(
      Value_type,
      levels = c("Observed LFQ values", "KNN-imputed MAR values", "QRILC-imputed MNAR values")
    )
  )

distribution_summary <- imputation_plot_df %>%
  group_by(Value_type) %>%
  summarise(
    n = n(),
    median_intensity = median(Intensity),
    mean_intensity = mean(Intensity),
    .groups = "drop"
  )

distribution_summary

p<-ggplot(imputation_plot_df, aes(x = Intensity)) +
  geom_density(fill = "grey70", alpha = 0.6, linewidth = 0.8) +
  geom_vline(
    data = distribution_summary,
    aes(xintercept = median_intensity),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  facet_wrap(~Value_type, ncol = 1, scales = "fixed") +
  labs(
    x = "Log2 LFQ intensity",
    y = "Density",
    title = "Distribution of observed and imputed LFQ values"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(paste0(treatment_name,"_imputation.png"), plot = p, width = 7, height = 4, dpi = 300)

## ---- 9. PCA (full imputed matrix) -----------------------------------------

pca_input <- lfq_imputed

protein_variance <- apply(pca_input, 1, var, na.rm = TRUE)

pca_input <- pca_input[
  is.finite(protein_variance) & protein_variance > 0,
  ,
  drop = FALSE
]

pca <- prcomp(t(pca_input), center = TRUE, scale. = FALSE)

variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  Sample = rownames(pca$x),
  Group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  stringsAsFactors = FALSE
)

pca_df$Sample <- gsub("\\s*MaxLFQ\\s*Intensity", "", pca_df$Sample)

ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = Sample), size = 4, show.legend = FALSE) +
  labs(
    x = paste0("PC1 (", round(variance_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(variance_explained[2], 1), "%)")
  ) +
  theme_classic() +
  theme(legend.title = element_blank(), legend.position = "right")

# (optional) PC1 & PC2 group tests
pc1_test <- summary(aov(PC1 ~ Group, data = pca_df))
pc2_test <- summary(aov(PC2 ~ Group, data = pca_df))
pc1_test
pc2_test

# (optional) within-group PCA dispersion test
pca_df <- pca_df %>%
  group_by(Group) %>%
  mutate(
    centroid_PC1 = mean(PC1),
    centroid_PC2 = mean(PC2),
    distance_to_centroid = sqrt((PC1 - centroid_PC1)^2 + (PC2 - centroid_PC2)^2)
  ) %>%
  ungroup()

pca_df

var.test(distance_to_centroid ~ Group, data = pca_df)
wilcox.test(distance_to_centroid ~ Group, data = pca_df)

pca_df %>%
  group_by(Group) %>%
  summarise(
    mean_distance = mean(distance_to_centroid),
    sd_distance = sd(distance_to_centroid),
    variance_distance = var(distance_to_centroid),
    n = n()
  )

## ---- 10. Differential expression (limma, excluding strict condition-specific) ----

condition_specific <- union(treatment_specific, control_specific)
condition_specific_found <- intersect(condition_specific, rownames(lfq_imputed))

lfq_limma <- lfq_imputed[
  !rownames(lfq_imputed) %in% condition_specific_found,
  ,
  drop = FALSE
]

write.csv(
  lfq_limma,
  paste0(treatment_name, "_vs_", control_name, "_KNN_QRILC_excluding_strict_specific_for_limma.csv")
)

group <- factor(
  c(rep(treatment_name, length(treatment_cols)),
    rep(control_name, length(control_cols))),
  levels = c(control_name, treatment_name)
)

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design

contrast_formula <- paste0(treatment_name, " - ", control_name)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
contrast_matrix

fit <- lmFit(lfq_limma, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)

results <- topTable(fit, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
results$Protein <- rownames(results)
rownames(results) <- NULL

write.csv(
  results,
  paste0(treatment_name, "_vs_", control_name, "_limma_results.csv"),
  row.names = FALSE
)

## ---- 11. Volcano plot ------------------------------------------------------

results$logFC <- as.numeric(results$logFC)
results$P.Value <- as.numeric(results$P.Value)

fc_cutoff <- 0.59
p_cutoff <- 0.05
adjp_cutoff <- 0.05

results <- results %>%
  mutate(
    negLog10P = -log10(P.Value),
    Significant = case_when(
      logFC > fc_cutoff & adj.P.Val < adjp_cutoff ~ "BH-adjusted P < 0.05",
      logFC > fc_cutoff & P.Value < p_cutoff       ~ "P < 0.05",
      logFC < -fc_cutoff & P.Value < p_cutoff       ~ "Down",
      TRUE ~ "Not significant"
    )
  )

p <- ggplot(results, aes(x = logFC, y = negLog10P, color = Significant)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(
    values = c(
      "Down" = "#2C7BB6",
      "Not significant" = "grey70",
      "P < 0.05" = "orange",
      "BH-adjusted P < 0.05" = "#8B0000"
    ),
    breaks = c("BH-adjusted P < 0.05", "P < 0.05", "Down", "Not significant")
  ) +
  geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
  theme_classic(base_size = 14) +
  labs(
    title = "Volcano Plot",
    x = expression(log[2] ~ FC ~ "(rGAL vs UL)"),
    y = expression(-log[10] ~ "P-value"),
    color = NULL
  )

p

ggsave(paste0(treatment_name, "_Volcano_plot.png"), plot = p, width = 10, height = 6, dpi = 300)

## ---- 12. GSEA (limma on full imputed matrix, rank by moderated t) --------

lfq_gsea <- lfq_imputed

fit_gsea <- lmFit(lfq_gsea, design)
fit_gsea <- contrasts.fit(fit_gsea, contrast_matrix)
fit_gsea <- eBayes(fit_gsea, trend = TRUE, robust = TRUE)

results_gsea <- topTable(fit_gsea, coef = 1, number = Inf, adjust.method = "BH", sort.by = "P")
results_gsea$Protein <- rownames(results_gsea)
rownames(results_gsea) <- NULL

ranked_list <- results_gsea$t
names(ranked_list) <- results_gsea$Protein

ranked_list <- ranked_list[
  is.finite(ranked_list) & !is.na(names(ranked_list)) & names(ranked_list) != ""
]

ranked_list <- tapply(ranked_list, names(ranked_list), function(x) x[which.max(abs(x))])
ranked_list <- sort(unlist(ranked_list), decreasing = TRUE)

run_gsea_collection <- function(collection_name, msig_table, ranked_stats, output_prefix,
                                min_size = 5, max_size = 500) {
  
  pathways <- split(msig_table$gene_symbol, msig_table$gs_name)
  pathways <- lapply(pathways, unique)
  
  set.seed(12345)
  gsea_result <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranked_stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0
  )
  
  gsea_result <- as.data.frame(gsea_result) %>%
    mutate(
      Collection = collection_name,
      Enriched_in = case_when(
        NES > 0 ~ treatment_name,
        NES < 0 ~ control_name,
        TRUE ~ NA_character_
      ),
      leadingEdge = vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
    ) %>%
    arrange(pval, desc(abs(NES)))
  
  significant_result <- gsea_result %>% filter(!is.na(pval), pval < 0.05)
  
  write.csv(
    significant_result,
    paste0(output_prefix, "_", collection_name, "_GSEA_nominal_P005.csv"),
    row.names = FALSE
  )
  
  cat(collection_name, ":", nrow(significant_result), "terms with nominal P < 0.05\n")
  
  return(list(all_results = gsea_result, significant_results = significant_result))
}

gobp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gsea_gobp <- run_gsea_collection("GOBP", gobp, ranked_list, paste0(treatment_name, "_vs_", control_name))

gocc <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:CC")
gsea_gocc <- run_gsea_collection("GOCC", gocc, ranked_list, paste0(treatment_name, "_vs_", control_name))

gomf <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:MF")
gsea_gomf <- run_gsea_collection("GOMF", gomf, ranked_list, paste0(treatment_name, "_vs_", control_name))

c2_all <- msigdbr(species = "Homo sapiens", collection = "C2")
kegg <- c2_all %>% filter(grepl("KEGG", gs_subcollection, ignore.case = TRUE))
unique(kegg$gs_subcollection)
length(unique(kegg$gs_name))
gsea_kegg <- run_gsea_collection("KEGG", kegg, ranked_list, paste0(treatment_name, "_vs_", control_name))

gsea_global_p005 <- bind_rows(
  gsea_gobp$significant_results,
  gsea_gocc$significant_results,
  gsea_gomf$significant_results,
  gsea_kegg$significant_results
) %>%
  arrange(Collection, pval, desc(abs(NES)))

write.csv(
  gsea_global_p005,
  paste0(treatment_name, "_vs_", control_name, "_global_GSEA.csv"),
  row.names = FALSE
)
