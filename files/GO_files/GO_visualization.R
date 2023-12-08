#### Ally Howard
# MICB405 GO Term enrichment and visualization
# Nov 24, 2023

##### All requisite files are in micb405_team8/files/GO_files ####
# please download all files into the same directory before running the R script
# Or download the zipped file in micb405_team8

## Load packages
library(tidyverse)
library(topGO)

## Load data
nitrogen_DESeq_results <- read_csv("nitrogen_results.csv")
DESeq_by_log2fold <- nitrogen_DESeq_results[order(nitrogen_DESeq_results$log2FoldChange), ]
write.csv(DESeq_by_log2fold, "DESeq_nitrogen")

#### Upregulated data ####
pos_log2fold_GO <- read_csv("pos_log2fold_GO.csv")

# Visualize
up_GO_filtered <- pos_log2fold_GO %>%
  mutate(GeneRatio = intersection_size/query_size, adjusted_p_value = as.numeric(adjusted_p_value)) %>%
  filter(adjusted_p_value < 0.05) %>%
  head(15)

sorted_up_GO <- up_GO_filtered[order(up_GO_filtered$GeneRatio, decreasing = TRUE),]

sorted_up_GO %>%
  ggplot(aes(x = fct_reorder(term_name, GeneRatio), y = GeneRatio, colour = adjusted_p_value)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = intersection_size)) +
  coord_flip() +
  scale_colour_gradient(low = "red", high = "blue") +
  ggtitle("Upregulated Expression in O. sativa at High Nitrogen Concentrations") + 
  labs(y = "Enrichment Ratio", x = "GO Term Description") + 
  labs(colour = "P-value") +
  labs(size = "Number of Significant Genes")

#### Downregulated data ####
neg_log2fold_GO <- read_csv("neg_log2fold_GO.csv")

# Visualize
down_GO_filtered <- neg_log2fold_GO %>%
  mutate(GeneRatio = intersection_size/query_size, adjusted_p_value = as.numeric(adjusted_p_value)) %>%
  filter(adjusted_p_value < 0.05) %>%
  head(15)

sorted_down_GO <- down_GO_filtered[order(down_GO_filtered$GeneRatio, decreasing = TRUE),]

sorted_down_GO %>%
  ggplot(aes(x = fct_reorder(term_name, GeneRatio), y = GeneRatio, colour = adjusted_p_value)) +
  geom_col(width = 0.05) +
  geom_point(aes(size = intersection_size)) +
  coord_flip() +
  scale_colour_gradient(low = "red", high = "blue") +
  ggtitle("Downregulated Expression in O. sativa at High Nitrogen Concentrations") + 
  labs(y = "Enrichment Ratio", x = "GO Term Description") + 
  labs(colour = "P-value") +
  labs(size = "Number of Significant Genes")

#### make a plot with both up and downregulated genes ####
# Add labels to upregulated and downregulated dataframes
up_GO <- sorted_up_GO %>% 
  mutate(up_down = "UP")

down_GO <- sorted_down_GO %>% 
  mutate(up_down = "DOWN")

# Make a joined dataframe
joined_GO_filtered_arranged <- bind_rows(up_GO, down_GO) %>%
  filter(adjusted_p_value <= 0.05) %>%
  mutate(GeneRatio = intersection_size/query_size, adjusted_p_value = as.numeric(adjusted_p_value)) %>%
  arrange(GeneRatio) %>%
  mutate(Term = factor(term_name)) %>%
  head(n = 40)

# Extract the column order
order_term_joined <- joined_GO_filtered_arranged %>% 
  pull(Term)

joined_GO_filtered_arranged %>% 
  ggplot(aes(x= Term, y = GeneRatio, color = adjusted_p_value)) +
  geom_point(aes(size= intersection_size)) +
  coord_flip() +
  scale_x_discrete(limits = order_term_joined) +
  scale_color_gradient(low = "red", high = "blue") +
  theme_light() +
  labs(x = "GO Term Description", y = "Enrichment Ratio", color = "P-value", size = "Number of Significant Genes") +
  theme(panel.border = element_rect(color = "black"), 
        panel.grid = element_line(colour = "grey96"), 
        strip.background = element_rect(colour = "black")) +
  scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  facet_grid(.~ up_down) +
  ggtitle("Up and Downregulated Pathways in Rice Exposed to High Nitrogen") +
  theme(plot.title = element_text(hjust = 2.0))

