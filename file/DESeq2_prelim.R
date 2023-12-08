# Stella Lin 
# MICB405 DESeq2 Analysis Final Project
# Nov 24, 2023

# Load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# Loading count data
# Control replicate 1
ctrl_R1 <- read_tsv("aln_STAR/control_n1/control_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Control replicate 2
ctrl_R2 <- read_tsv("aln_STAR/control_n2/control_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Control replicate 3
ctrl_R3 <- read_tsv("aln_STAR/control_n3/control_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Gradient replicate 1
gradient_R1 <- read_tsv("aln_STAR/gradient_n1/gradient_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Gradient replicate 2
gradient_R2 <- read_tsv("aln_STAR/gradient_n2/gradient_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Gradient replicate 3
gradient_R3 <- read_tsv("aln_STAR/gradient_n3/gradient_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 1
low_R1 <- read_tsv("aln_STAR/low_n1/low_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 2
low_R2 <- read_tsv("aln_STAR/low_n2/low_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 3
low_R3 <- read_tsv("aln_STAR/low_n3/low_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)

# Making metadata
dat <- data.frame(
  row.names = ctrl_R1$gene_id,
  control_rep1 = ctrl_R1$sense,
  control_rep2 = ctrl_R2$sense,
  control_rep3 = ctrl_R3$sense,
  low_high_1h_rep1 = gradient_R1$sense,
  low_high_1h_rep2 = gradient_R2$sense,
  low_high_1h_rep3 = gradient_R3$sense,
  low_rep1 = low_R1$sense,
  low_rep2 = low_R2$sense,
  low_rep3 = low_R3$sense
)
dat_matrix <- as.matrix(dat) # Transform into matrix

# Set up metadata file to contain column information
metadata <- data.frame(row.names = colnames(dat_matrix),
                       condition = c("control", "control", "control", "low_high_1h", "low_high_1h", "low_high_1h", "low", "low", "low"))

#### Running DESeq2 ####

# Create DESeq2 object
dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix,
                                     colData = metadata,
                                     design = ~condition)

# Setting control condition as reference
dds_matrix$condition <- relevel(dds_matrix$condition, ref = "control")

# Running DESeq2 on dataset
dds <- DESeq(dds_matrix)
saveRDS(dds, "dds.rds") #save


#### QA: Sample Clustering ####

# Generate PCA plot
rld <- rlog(dds) # Log transformation on count data
pca_plot <- plotPCA(rld, intgroup = "condition") +
  labs(title = "PCA Plot")

# Generate Heatmap
sample_dists <- dist(t(assay(rld))) # Calculates distances between log-transformed data samples
sample_dist_matrix <- as.matrix(sample_dists) # Convert to matrix
colnames(sample_dist_matrix) <- NULL # Removes column names
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) # Colours!
heatmap <- pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colours)


#### Extracting data ####
res_low <- results(dds, name = "condition_low_vs_control") %>% as.data.frame() 
res_grad <- results(dds, name = "condition_low_high_1h_vs_control") %>% as.data.frame() 
res_comp <- results(dds, contrast = c("condition", "low_high_1h", "low"))
res_comp_df <- as.data.frame(res_comp)


#### Data wrangling ####
res_no_NA_low <- res_low %>%
  drop_na() # Filter out NA
res_no_NA_grad <- res_grad %>%
  drop_na()
res_no_NA <- res_comp_df %>%
  drop_na()

res_filtered <- res_no_NA %>%
  filter(padj <= 0.05) 

res_filtered_final <- res_filtered %>%
  filter(log2FoldChange <- -1 | log2FoldChange >= 1) %>%
  rownames_to_column("gene_id") 

top15_genes <- res_filtered_final %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 15)

bot15_genes <- res_filtered_final %>%
  arrange(log2FoldChange) %>%
  head(n = 15)

greatest30_genes <- bind_rows(top15_genes, bot15_genes) %>%
  mutate(up_down = if_else(log2FoldChange > 0, "Up", "Down")) %>%
  arrange(log2FoldChange)

order_30genes <- greatest30_genes %>%
  select(gene_id) %>%
  pull()

CSV <- write_csv(res_filtered_final, "nitrogen_results.csv")
csv_top <- write.csv(top10_genes, "nitrogen_top_results.csv")
csv_bot <- write.csv(bot10_genes, "nitrogen_bot_results.csv")

#### Plots ####
# Volcano plots
volcano_plot <- ggplot(res_filtered_final, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  labs(title = "Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(p-value)")


# Mean-Dispersion Plots
Mean_Dispersion <- plotDispEsts(dds)


# DESeq plot
DESeq_plot <- ggplot(greatest30_genes, aes(x = gene_id, y = log2FoldChange, fill = up_down)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.4) +  # Adjust width as needed
  scale_fill_manual(values = c("#9FE2BF", "#F08080"), name = "Regulation") +
  labs(title = "DESeq Bar Plot",
       x = "Gene ID",
       y = "log2(Fold Change)",
       fill = "Regulation") +
  coord_flip() +
  geom_text(aes(label = sprintf("%.2f", log2FoldChange)), vjust = -0.5, size = 3) 


DESeq_plot

# Selecting gene_IDs
selected_gene_ids <- c("#LOC9270637", "#LOC9270637", "#LOC4339749", "#LOC4348505")
selected_data <- greatest30_genes %>%
  filter(gene_id %in% selected_gene_ids)
DESeq_plot_selected <- ggplot(selected_data, aes(x = gene_id, y = log2FoldChange, fill = up_down)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.4) +
  scale_fill_manual(values = c("#F08080", "#9FE2BF"), name = "Regulation") +
  labs(title = "Up/Down Regulation of Gene ID's",
       x = "Gene ID",
       y = "log2(Fold Change)",
       fill = "Regulation") +
  coord_flip() +
  geom_text(aes(label = sprintf("%.2f", log2FoldChange)), vjust = -0.5, size = 3)

DESeq_plot_selected
