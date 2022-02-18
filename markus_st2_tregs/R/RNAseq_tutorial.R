# dataset
library(airway)

# packages
library(tidyverse)
library(plotly)
library(ggrepel)
library(GGally)
library(tidyHeatmap)
library(tidybulk)

# Use colourblind-friendly colours
friendly_cols <- dittoSeq::dittoColors()

# Set theme
custom_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        # legend.position = "bottom",
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
      )
  )

# load data
data("airway")
airway

airway@colData$SampleName %>% n_distinct()
airway@colData$cell %>% n_distinct()

# SummarizedExperiment data object does not work with tidyverse vocabulary directly
# A more convenient usage comes through tidySummarizedExperiment
library(tidySummarizedExperiment)
airway %>% 
  select(.sample) %>% 
  n_distinct()

# SummarizedExperiment functions can still work
assay(airway)

# tidyverse functions
airway %>% 
  filter(dex == "trt") %>% 
  select(.sample, cell, dex)

# check the total number of reads for each sample:
# the variation of total reads across samples necessitate normalisation to ensure gene expressions are comparable
airway %>% 
  group_by(.sample) %>% 
  summarise(total_counts = sum(counts))

airway %>% 
  mutate(sample_name = str_remove(.sample, "SRR1039")) %>% 
  select(.sample, sample_name)

# convert ensembl codes to symbols
counts <-
  airway %>%
  mutate(sample_name = str_remove(.sample, "SRR1039")) %>%
  mutate(symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                        keys = .feature,
                                        keytype = "ENSEMBL",
                                        column = "SYMBOL",
                                        multiVals = "first"
                                        )
         )

x <-
  airway %>%
  mutate(sample_name = str_remove(.sample, "SRR1039")) %>%
  mutate(symbol = AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                        keys = .feature,
                                        keytype = "ENSEMBL",
                                        column = "SYMBOL",
                                        multiVals = "first"
                                        )
  ) %>% 
  keep_abundant(factor_of_interest = dex) %>% 
  scale_abundance()
  

# filter lowly expressed transcripts
counts_filtered <- counts %>% 
  keep_abundant(factor_of_interest = dex)

# normalisation
counts_scaled <- counts_filtered %>% 
  scale_abundance()

counts_scaled %>% 
  pivot_longer(cols = c("counts", "counts_scaled"), 
               names_to = "source", values_to = "abundance") %>% 
  ggplot(aes(x = abundance + 1, color = sample_name)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10() +
  custom_theme

# dimension reduction
counts_scaled_PCA <- counts_scaled %>% 
  reduce_dimensions(method = "PCA")

# PCA plot
counts_scaled_PCA %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = dex, shape = cell)) +
  geom_point() +
  geom_text_repel(aes(label = sample_name), show.legend = FALSE) +
  custom_theme

# hierarchical clustering with heatmap
counts_scaled %>% 
  # extract 500 most variable genes 
  # why would counts still have transcripts less than 10 after keep_abundant??
  keep_variable(.abundance = counts_scaled, top = 500) %>%
  as_tibble() %>%
  
  # create heatmap
  heatmap(
    .column = sample_name,
    .row = .feature,
    .value = counts_scaled,
    transform = log1p
  ) %>%
  add_tile(dex) %>%
  add_tile(cell)

# Differential gene expression
de_all <- counts_scaled %>%
  
  # edgeR QLT
  test_differential_abundance(
    ~ dex + cell,
    method = "edgeR_quasi_likelihood",
    prefix = "edgerQLT_"
  ) %>%
  
  # edgeR LRT
  test_differential_abundance(
    ~ dex + cell,
    method = "edgeR_likelihood_ratio",
    prefix = "edgerLR_"
  ) %>%
  
  # limma-voom
  test_differential_abundance(
    ~ dex + cell,
    method = "limma_voom",
    prefix = "voom_"
  ) %>%
  
  # DESeq2
  test_differential_abundance(
    ~ dex + cell,
    method = "deseq2",
    prefix = "deseq2_"
  )

de_all %>%
  pivot_transcript() %>%
  select(edgerQLT_PValue, edgerLR_PValue, voom_P.Value, deseq2_pvalue, feature) %>%
  ggpairs(1:4)

de_all %>%
  pivot_transcript() %>%
  select(contains(c("FDR", "adj")), feature) %>% 
  ggpairs(1:4)

counts_de <- counts_scaled %>% 
  test_differential_abundance(~ dex + cell)

# volcano plot, minimal
counts_de %>%
  ggplot(aes(x = logFC, y = PValue, colour = FDR < 0.05)) +
  geom_point() +
  scale_y_continuous(trans = "log10_reverse") +
  custom_theme

top6symbols <- counts_de %>% 
  pivot_transcript() %>% 
  slice_min(PValue, n = 6) %>% 
  pull(symbol)

counts_de %>% 
  pivot_transcript() %>% 
  select(1:8) %>% 
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>% 
  mutate(symbol = ifelse(symbol %in% top6symbols, as.character(symbol), "")) %>% 
  ggplot(aes(logFC, PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() + 
  scale_y_continuous(trans = "log10_reverse") +
  custom_theme +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2))

get_bibliography(counts_de)  

# cell type composition analysis
