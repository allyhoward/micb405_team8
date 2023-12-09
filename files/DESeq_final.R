# Stella Lin 
# MICB405 DESeq2 Analysis Final Project - No Gradient version
# Dec 7, 2023

## Reran without including gradient data - realized to not be relevant to our experiment
## Control = high; changed low to be control as it represented crops with no added fertilizer and decided to study high instead

# Load packages
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(dplyr))

# Loading count data
# High replicate 1
high_R1 <- read_tsv("aln_STAR/control_n1/control_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# High replicate 2
high_R2 <- read_tsv("aln_STAR/control_n2/control_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# High replicate 3
high_R3 <- read_tsv("aln_STAR/control_n3/control_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 1
low_R1 <- read_tsv("aln_STAR/low_n1/low_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 2
low_R2 <- read_tsv("aln_STAR/low_n2/low_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 3
low_R3 <- read_tsv("aln_STAR/low_n3/low_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)

# Making metadata
dat_new <- data.frame(
  row.names = high_R1$gene_id,
  control_rep1 = high_R1$sense,
  control_rep2 = high_R2$sense,
  control_rep3 = high_R3$sense,
  low_rep1 = low_R1$sense,
  low_rep2 = low_R2$sense,
  low_rep3 = low_R3$sense
)
dat_matrix_new <- as.matrix(dat_new) # Transform into matrix

# Set up metadata file to contain column information
metadata_new <- data.frame(row.names = colnames(dat_matrix_new),
                       condition = c("high", "high", "high", "low", "low", "low"))

#### Running DESeq2 ####

# Create DESeq2 object
dds_matrix_new <- DESeqDataSetFromMatrix(countData = dat_matrix_new,
                                     colData = metadata_new,
                                     design = ~condition)

# Setting control condition as reference
dds_matrix_new$condition <- relevel(dds_matrix_new$condition, ref = "low")

# Running DESeq2 on dataset
dds_new <- DESeq(dds_matrix_new)
saveRDS(dds_new, "dds_new.rds") #save


#### QA: Sample Clustering ####

# Generate PCA plot
rld_new <- rlog(dds_new) # Log transformation on count data
pca_plot_new <- plotPCA(rld_new, intgroup = "condition") +
  labs(title = "PCA Plot")

# Generate Heatmap
sample_dists_new <- dist(t(assay(rld_new))) # Calculates distances between log-transformed data samples
sample_dist_matrix_new <- as.matrix(sample_dists_new) # Convert to matrix
colnames(sample_dist_matrix_new) <- NULL # Removes column names
colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255) # Colours!
heatmap_new <- pheatmap(sample_dist_matrix_new,
                    clustering_distance_rows = sample_dists_new,
                    clustering_distance_cols = sample_dists_new,
                    col = colours)

#### Extracting data ####
res_new <- results(dds_new, name = "condition_low_vs_high") %>% as.data.frame() 


#### Data wrangling ####
res_no_NA_new <- res_new %>%
  drop_na() # Filter out NA


res_filtered_new <- res_no_NA_new %>%
  filter(padj <= 0.05) 

res_filtered_final_new <- res_filtered_new %>%
  filter(log2FoldChange <- -1 | log2FoldChange >= 1) %>%
  rownames_to_column("gene_id") 

top15_genes_new <- res_filtered_final_new %>%
  arrange(desc(log2FoldChange)) %>%
  head(n = 15)

bot15_genes_new <- res_filtered_final_new %>%
  arrange(log2FoldChange) %>%
  head(n = 15)

greatest30_genes_new <- bind_rows(top15_genes_new, bot15_genes_new) %>%
  mutate(up_down = if_else(log2FoldChange > 0, "Up", "Down")) %>%
  arrange(log2FoldChange)

order_30genes_new <- greatest30_genes_new %>%
  select(gene_id) %>%
  pull()

CSV_new <- write_csv(res_filtered_final_new, "nitrogen_results_new.csv")
csv_top_new <- write.csv(top10_genes_new, "nitrogen_top_results_new.csv")
csv_bot_new <- write.csv(bot10_genes_new, "nitrogen_bot_results_new.csv")

#### Plots ####
# Volcano plots
volcano_plot_new <- ggplot(res_filtered_final_new, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point() + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  labs(title = "Volcano Plot",
       x = "log2(Fold Change)",
       y = "-log10(p-value)")


# Mean-Dispersion Plots
Mean_Dispersion_new <- plotDispEsts(dds_new)


# DESeq plot
sorted_30_genes_new <- greatest30_genes_new[order(-greatest30_genes_new$log2FoldChange), ] # sort by log2 Fold Change

DESeq_plot_new <- ggplot(greatest30_genes_new, aes(x = gene_id, y = log2FoldChange, fill = up_down)) +
  geom_bar(stat="identity") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), width = 0.4) +  # Adjust width as needed
  scale_fill_manual(values = c("#9FE2BF", "#F08080"), name = "Regulation") +
  labs(title = "DESeq Bar Plot",
       x = "Gene ID",
       y = "log2(Fold Change)" +
       fill = "Regulation") +
  coord_flip() +
  geom_text(aes(label = sprintf("%.2f", log2FoldChange)), vjust = -0.5, size = 3) 


DESeq_plot
