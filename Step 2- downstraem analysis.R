setwd("D:/proteomics/proteomics")

############################################################
# 0. Setup
############################################################

library(tidyverse)
library(readr)
library(limma)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(readxl)
library(imputeLCMD)
library(ggrepel)


protein_file <- "D:/proteomics/proteomics/proteinGroups.txt"
out_dir      <- "D:/proteomics/proteomics"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

############################################################
# 1. Read proteinGroups and identify LFQ columns
############################################################

prot <- read_tsv(protein_file, show_col_types = FALSE)

# LFQ columns
intensity_prefix <- "LFQ intensity "
lfq_cols <- colnames(prot)[startsWith(colnames(prot), intensity_prefix)]
stopifnot(length(lfq_cols) > 0)

sample_names <- str_remove(lfq_cols, fixed(intensity_prefix)) %>%
  str_replace_all(" ", "")

print(sample_names)


############################################################
# 2. Create sample metadata 
############################################################

meta <- tibble(
  sample = sample_names,
  time = factor(str_extract(sample_names, "d0|d5|d10"),
                levels = c("d0", "d5", "d10")),
  replicate = as.numeric(str_extract(sample_names, "(?<=rep)[0-9]+"))
)

print(meta)


###############################################################
#                 5. Create LFQ matrix
###############################################################
lfq_mat <- prot %>% select(all_of(lfq_cols)) %>% as.matrix()
mode(lfq_mat) <- "numeric"

rownames(lfq_mat) <- prot$`Protein IDs`
colnames(lfq_mat) <- sample_names

###############################################################
#                 6. Filter bad proteins
###############################################################
# Remove reverse hits & contaminants
prot_filt <- prot %>%
  filter(is.na(`Potential contaminant`)) %>%
  filter(is.na(`Reverse`)) %>%
  filter(is.na(`Only identified by site`))

lfq_mat <- prot_filt %>% select(all_of(lfq_cols)) %>% as.matrix()
mode(lfq_mat) <- "numeric"
rownames(lfq_mat) <- prot_filt$`Protein IDs`

###############################################################
#                 7. Missing value filtering
###############################################################
colnames(lfq_mat)
meta$sample

colnames(lfq_mat) <- colnames(lfq_mat) %>%
  str_remove("LFQ intensity ") %>%
  str_replace_all(" ", "")
colnames(lfq_mat)


# keep proteins quantified in >= 2 replicates in any group
keep <- rowSums(!is.na(lfq_mat[, meta$sample[meta$time=="d0"]])) >= 2 |
  rowSums(!is.na(lfq_mat[, meta$sample[meta$time=="d5"]])) >= 2 |
  rowSums(!is.na(lfq_mat[, meta$sample[meta$time=="d10"]])) >= 2

lfq_mat <- lfq_mat[keep, ]


dim(lfq_mat)        # number of proteins × number of samples
rownames(lfq_mat)[1:10]  # first 10 protein IDs
colnames(lfq_mat)        # sample names
summary(lfq_mat)         # min, max, median, etc.

sum(lfq_mat == 0)   # how many zero values
sum(is.na(lfq_mat)) # how many NA values

head(lfq_mat)       # first 6 rows
lfq_mat[1:10, 1:5]  # first 10 proteins × first 5 samples

boxplot(lfq_mat, main = "LFQ intensities per sample", las = 2)

###############################################################
#                 8. Log2 transform & & QRILC imputation
###############################################################
# Replace zeros with NA
lfq_mat[lfq_mat == 0] <- NA

# Keep proteins with at least 1 observed value
lfq_mat <- lfq_mat[rowSums(!is.na(lfq_mat)) > 0, ]

# Log2 transform
lfq_log <- log2(lfq_mat)

# Ensure numeric matrix
lfq_log <- as.matrix(lfq_log)
mode(lfq_log) <- "numeric"

# Impute missing values (QRILC)
imp_res <- impute.QRILC(lfq_log)   # returns a list
# Get the imputed matrix
lfq_imputed <- imp_res[[1]]   # first element of the list

class(lfq_imputed)            # should be "matrix"
dim(lfq_imputed)              # (proteins × samples)
sum(is.na(lfq_imputed))       # should be 0

boxplot(lfq_imputed, main = "LFQ intensities per sample", las = 2)

###############################################################
#                 9. PCA plot
###############################################################


pca <- prcomp(t(lfq_imputed), scale.=TRUE)

pca_df <- data.frame(
  PC1 = pca$x[,1],
  PC2 = pca$x[,2],
  time = meta$time,
  sample = meta$sample
)

ggplot(pca_df, aes(PC1, PC2, color = time, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel() +
  theme_bw()
ggsave(file.path(out_dir, "PCA_plot.png"))


#OR

# Define custom colors for your time points (example)
custom_colors <- c("blue", "green", "orange")  # adjust based on your time points

ggplot(pca_df, aes(PC1, PC2, color = time, label = sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 5, box.padding = 0.5, point.padding = 0.3) +
  scale_color_manual(values = custom_colors) +
  theme_minimal(base_size = 16) +
  theme(
    panel.grid = element_blank(),           # remove grid
    axis.line = element_line(color = "black", size = 0.8),  # show axis lines
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 14)
  ) +
  labs(title = "PCA Plot", x = "PC1", y = "PC2", color = "Time")

# Save the plot
ggsave(file.path(out_dir, "PCA_plot.png"), width = 8, height = 6)


###############################################################
#        10. Sample correlation heatmap
###############################################################
cor_mat <- cor(lfq_imputed)

pheatmap(cor_mat, clustering_distance_cols = "correlation")
ggsave(file.path(out_dir, "correlation_heatmap.png"))



###############################################################
# 10. limma design and contrasts
###############################################################

meta$time <- factor(meta$time, levels = c("d0", "d5", "d10"))
stopifnot(all(colnames(lfq_imputed) == meta$sample))

design <- model.matrix(~ 0 + time, data = meta)
colnames(design) <- levels(meta$time)
design

contrast.matrix <- makeContrasts(
  d5_vs_d0  = d5 - d0,
  d10_vs_d0 = d10 - d0,
  d10_vs_d5 = d10 - d5,
  levels = design
)
contrast.matrix

###############################################################
# 11. Fit limma model
###############################################################

fit  <- lmFit(lfq_imputed, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

###############################################################
# 12. Extract differential proteins
###############################################################

res_d5_vs_d0  <- topTable(fit2, coef = "d5_vs_d0",  number = Inf) %>%
  rownames_to_column("Protein_ID")
res_d10_vs_d0 <- topTable(fit2, coef = "d10_vs_d0", number = Inf) %>%
  rownames_to_column("Protein_ID")

write_tsv(res_d5_vs_d0,  file.path(out_dir, "DE_d5_vs_d0.tsv"))
write_tsv(res_d10_vs_d0, file.path(out_dir, "DE_d10_vs_d0.tsv"))


# Read the same proteinGroups.txt you used earlier
prot_annot <- read_tsv(protein_file, show_col_types = FALSE) %>%
  select(`Protein IDs`, `Gene names`) %>%
  distinct()

# Add gene names to results by Protein_ID
res_d5_vs_d0  <- res_d5_vs_d0  %>%
  left_join(prot_annot, by = c("Protein_ID" = "Protein IDs"))

res_d10_vs_d0 <- res_d10_vs_d0 %>%
  left_join(prot_annot, by = c("Protein_ID" = "Protein IDs"))

# Now res_* tables have a "Gene names" column
write_tsv(res_d5_vs_d0,  file.path(out_dir, "DE_d5_vs_d0_withGenes.tsv"))
write_tsv(res_d10_vs_d0, file.path(out_dir, "DE_d10_vs_d0_withGenes.tsv"))




###############################################################
# 13. Volcano plots
###############################################################

plot_volcano <- function(res, title, fc_cut = 1, p_cut = 0.05) {
  res2 <- res %>%
    mutate(
      sig = case_when(
        adj.P.Val < p_cut & logFC >  fc_cut ~ "Up",
        adj.P.Val < p_cut & logFC < -fc_cut ~ "Down",
        TRUE                                ~ "NS"
      )
    )
  
  ggplot(res2, aes(logFC, -log10(adj.P.Val), color = sig)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c(Up = "red", Down = "blue", NS = "grey70")) +
    geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_cut), linetype = "dashed") +
    theme_bw() +
    labs(title = title, x = "log2 fold-change", y = "-log10 adj. P-value")
}

pdf(file.path(out_dir, "volcano_plots_LFQ.pdf"), width = 10, height = 4)
print(plot_volcano(res_d5_vs_d0,  "d5 vs d0"))
print(plot_volcano(res_d10_vs_d0, "d10 vs d0"))
dev.off()

###############################################################
# 14. Save imputed matrix
###############################################################

lfq_imputed_df <- as.data.frame(lfq_imputed) %>%
  rownames_to_column("Protein_ID")

write_tsv(lfq_imputed_df, file.path(out_dir, "LFQ_log2_imputed_matrix.tsv"))

