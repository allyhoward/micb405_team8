### Vera Pu ###
### Upset Plot ###
### GOAL: filter out all the upregulated/downregulated genes in each condition and compare their intersections ###

## Load Packages
library(tidyverse)
library(DESeq2)
library(UpSetR)

## Load Data
# Control replicate 1
ctrl_R1 <- read_tsv("aln_STAR/control_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Control replicate 2
ctrl_R2 <- read_tsv("aln_STAR/control_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Control replicate 3
ctrl_R3 <- read_tsv("aln_STAR/control_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Gradient replicate 1
gradient_R1 <- read_tsv("aln_STAR/gradient_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Gradient replicate 2
gradient_R2 <- read_tsv("aln_STAR/gradient_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Gradient replicate 3
gradient_R3 <- read_tsv("aln_STAR/gradient_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 1
low_R1 <- read_tsv("aln_STAR/low_R1ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 2
low_R2 <- read_tsv("aln_STAR/low_R2ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)
# Low replicate 3
low_R3 <- read_tsv("aln_STAR/low_R3ReadsPerGene.out.tab", col_names = c("gene_id", "total", "antisense", "sense"), skip = 4)

## Making metadata
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

## Running DESeq2
# Create DESeq2 object
dds_matrix <- DESeqDataSetFromMatrix(countData = dat_matrix,
                                     colData = metadata,
                                     design = ~condition)
# Setting control condition as reference
dds_matrix$condition <- relevel(dds_matrix$condition, ref = "control")
# Running DESeq2 on dataset
dds <- DESeq(dds_matrix)

## Extracting data
res_low <- results(dds, name = "condition_low_vs_control") %>% as.data.frame() # for low
res_grad <- results(dds, name = "condition_low_high_1h_vs_control") %>% as.data.frame() # for gradient
res_ctrl <- results(dds, contrast = c("condition", "low_high_1h", "low")) %>% as.data.frame() # for control

## Data wrangling
# low_upregulated
res_no_NA_low <- res_low %>%
  drop_na() # Filter out NA
res_filtered_low <- res_no_NA_low %>%
  filter(padj <= 0.05)
res_filtered_final_low_up <- res_filtered_low %>%
  filter(log2FoldChange >= 1) %>%
  rownames_to_column("gene_id")
up_low_genes <- res_filtered_final_low_up %>% 
  pull("gene_id") 
# low_downregulated
res_no_NA_low <- res_low %>%
  drop_na() # Filter out NA
res_filtered_low <- res_no_NA_low %>%
  filter(padj <= 0.05)
res_filtered_final_low_down <- res_filtered_low %>%
  filter(log2FoldChange < 1) %>%
  rownames_to_column("gene_id")
down_low_genes <- res_filtered_final_low_down %>% 
  pull("gene_id") 
# gradient_upregulated
res_no_NA_grad <- res_grad %>%
  drop_na()
res_filtered_grad <- res_no_NA_grad %>%
  filter(padj <= 0.05)
res_filtered_final_grad_up <- res_filtered_grad %>%
  filter(log2FoldChange >= 1) %>%
  rownames_to_column("gene_id")
up_grad_genes <- res_filtered_final_grad_up %>% 
  pull("gene_id") 
# gradient_downregulated
res_no_NA_grad <- res_grad %>%
  drop_na()
res_filtered_grad <- res_no_NA_grad %>%
  filter(padj <= 0.05)
res_filtered_final_grad_down <- res_filtered_grad %>%
  filter(log2FoldChange < 1) %>%
  rownames_to_column("gene_id")
down_grad_genes <- res_filtered_final_grad_down %>% 
  pull("gene_id") 
# control_upregulated
res_no_NA_ctrl <- res_ctrl %>%
  drop_na() # Filter out NA
res_filtered_ctrl <- res_no_NA_ctrl %>%
  filter(padj <= 0.05)
res_filtered_final_ctrl_up <- res_filtered_ctrl %>%
  filter(log2FoldChange >= 1) %>%
  rownames_to_column("gene_id")
up_ctrl_genes <- res_filtered_final_ctrl_up %>% 
  pull("gene_id") 
# ctrl_downregulated
res_no_NA_ctrl <- res_ctrl %>%
  drop_na() # Filter out NA
res_filtered_ctrl <- res_no_NA_ctrl %>%
  filter(padj <= 0.05)
res_filtered_final_ctrl_down <- res_filtered_ctrl %>%
  filter(log2FoldChange < 1) %>%
  rownames_to_column("gene_id")
down_ctrl_genes <- res_filtered_final_ctrl_down %>% 
  pull("gene_id") 

## Make a list of all features then turn into dataframe
total_features_up <- list(gradient = up_grad_genes, low = up_low_genes, control = up_ctrl_genes)
feature_df_up <- stack(total_features_up)
total_features_down <- list(gradient = down_grad_genes, low = down_low_genes, control = down_ctrl_genes)
feature_df_down <- stack(total_features_down)

## Turn into a binary dataframe
# upregulated
binary_df_up <- feature_df_up %>% 
  mutate(value = 1) %>% 
  complete(ind, values = unique(unlist(total_features_up)), fill = list(value = 0)) %>%
  pivot_wider(names_from = ind, values_from = value) %>% 
  as.data.frame() 
# downregulated
binary_df_down <- feature_df_down %>% 
  mutate(value = 1) %>% 
  complete(ind, values = unique(unlist(total_features_down)), fill = list(value = 0)) %>%
  pivot_wider(names_from = ind, values_from = value) %>% 
  as.data.frame()

## Generate upset plot
upset(binary_df_up, sets = c("gradient", "low", "control")) 
upset(binary_df_down, sets = c("gradient", "low", "control"))
