library(tidyverse)
library(tidyHeatmap)
library(tidybulk)
library(readxl)
library(ggupset)
library(ggrepel)

# untar("GSE130879_RAW.tar")
GSM3755558_4_Fat_Treg_genes_tsv <- read_delim("GSM3755558_4_Fat_Treg_genes.tsv.gz", 
                                              delim = "\t", escape_double = FALSE, trim_ws = TRUE)

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

counts <- read_xlsx("data/GSE130842_Count_table_Delacher_et_al_2019.xlsx")
rpkm <- read_xlsx("data/GSE130842_Rpkm_table_Delacher_et_al_2019.xlsx")

counts_preprocessed <- 
  # annotate genes and add tissue and mouse columns
  counts %>% 
  left_join(rpkm %>% select(ID, symbol), ., by = "ID") %>% 
  # mutate(symbol = AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,
  #                                       keys = ID,
  #                                       keytype = "ENSEMBL",
  #                                       column = "SYMBOL",
  #                                       multiVals = "first"
  #                                       ), .after = ID)
  select(ID, symbol, contains("tisTregST2")) %>% 
  pivot_longer(-c(ID, symbol), names_to = "sample", values_to = "count") %>% 
  mutate(tissue = str_extract(sample, "(?<=_).+(?=_)")) %>% 
  mutate(mouse = str_extract(sample, "R\\d+")) %>% 
  tidybulk(sample, symbol, count) %>% 
  
  # create contrasts for different analysis
  mutate(analysis = tissue) %>% 
  nest(data=-analysis) %>% 
  mutate(analysis2 = analysis) %>%
  mutate(data2 = data) %>% 
  tidyr::expand(nesting(analysis, data), nesting(analysis2, data2)) %>% 
  filter(analysis != analysis2) %>% 
  mutate(contrast = map2_chr(
    analysis, analysis2, 
    ~ sort(c(.x, .y)) %>% 
      paste0("tissue", .) %>% 
      paste(collapse = " - ")), 
    .before = 1) %>% 
  distinct(contrast, .keep_all = TRUE) %>% 
  mutate(data = map2(data, data2, ~ bind_rows(.x, .y))) %>% 
  select(contrast, data) %>% 
  
  # preprocessing
  mutate(data = map(
    data,
    ~.x %>% 
      aggregate_duplicates() %>% 
      keep_abundant(factor_of_interest = tissue) %>% 
      scale_abundance()
  ))

# Differential expression analysis
de_all <- counts_preprocessed %>% 
  
  mutate(data = map2(
    contrast, data,
    ~ .y %>% 
      test_differential_abundance(
        ~ 0 + tissue + mouse,
        .abundance = count_scaled,
        # # When you use pipes the object on the left is the first input to the function by default. To stop that use curly braces 
        # .contrasts = .x %>% distinct(tissue) %>% pull %>% {sprintf("tissue%s - tissue%s", .[1], .[2])},
        .contrast = .x,
        omit_contrast_in_colnames = TRUE,
        action = "get"
      )
  ))

job::job({de_all %>% saveRDS("data/de_all.rds", compress = "xz")})

# filter significantly DE genes
de_all_filtered <- de_all %>% 
  mutate(data = map(
    data,
    ~.x %>% 
      filter(FDR < 0.05 & abs(logFC) > 2) %>%
      filter(logCPM > mean(logCPM)) %>% 
      arrange(desc(logFC))
  ))

saveRDS(de_all_filtered, "intermediate_data/de_all_filtered.rds", compress = "xz")

# Volcano plot
# Fat vs Spleen
top_genes_vatspl <- de_all_filtered %>% 
  pluck("data", 9) %>% 
  head(30) %>% 
  pull(symbol)

de_all %>% 
  pluck("data", 9) %>% 
  # Subset data
  mutate(significant = (FDR < 0.05) & (abs(logFC) >= 2) & (logCPM > mean(logCPM))) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_vatspl, symbol, "")) %>%
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 20) +
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "VAT vs Spleen") +
  theme(plot.title = element_text(hjust = 0.5))

# Skin vs Spleen
top_genes_skinspl <- de_all_filtered %>% 
  pluck("data", 15) %>% 
  head(30) %>% 
  pull(symbol)

de_all %>% 
  pluck("data", 15) %>% 
  # Subset data
  mutate(significant = (FDR < 0.05) & (abs(logFC) >= 2) & (logCPM > mean(logCPM))) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_skinspl, symbol, "")) %>%
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 20) +
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "Skin vs Spleen") +
  theme(plot.title = element_text(hjust = 0.5))

# Fat vs Skin
top_genes_vatskin <- de_all_filtered %>% 
  pluck("data", 8) %>% 
  head(30) %>% 
  pull(symbol)

de_all %>% 
  pluck("data", 8) %>% 
  # Subset data
  mutate(significant = (FDR < 0.05) & (abs(logFC) >= 2) & (logCPM > mean(logCPM))) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_vatskin, symbol, "")) %>%
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 20) +
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "Fat vs Skin") +
  theme(plot.title = element_text(hjust = 0.5))

# ggupset plot
de_all_filtered %>% 
  unnest(data) %>% 
  filter(logFC>0) %>% 
  group_by(symbol) %>% 
  summarise(enriched_tissue = list(contrast)) %>% 
  ggplot(aes(enriched_tissue)) +
  geom_bar() +
  scale_x_upset() +
  labs(title="genes enriched in target tissue") +
  theme(plot.title = element_text(hjust = 0.5))

de_all_filtered %>% 
  unnest(data) %>% 
  filter(logFC>0) %>% 
  group_by(symbol) %>% 
  summarise(enriched_tissue = list(contrast)) %>% 
  ggplot(aes(enriched_tissue)) +
  geom_bar() +
  scale_x_mergelist(sep = "_") +
  axis_combmatrix(sep = "_")

de_all_filtered %>% 
  unnest(data) %>% 
  filter(logFC>0) %>% 
  group_by(symbol) %>% 
  summarise(enriched_tissue = list(contrast)) %>% 
  mutate(n_tissue = lengths(enriched_tissue)) %>% 
  filter(n_tissue==1) %>% 
  unnest(enriched_tissue) %>% 
  separate(enriched_tissue, c("target", "background"), sep = " - ") %>% 
  filter(target == "tissueFat" & background == "tissueSpleen") %>% 
  pull(symbol) %>% 
  cat(sep = "\t")

# PCA analysis
pca_all <- counts_preprocessed %>% 
  mutate(data = map(
    data,
    ~.x %>% 
      reduce_dimensions(
        method = "PCA",
        .abundance = count_scaled,
        action = "get"
      )
  ))

# liver vs lung
pca_all %>% 
  pluck("data", 10) %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Fat vs skin
pca_all %>% 
  pluck("data", 8) %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Liver vs lung
pca_all %>% 
  pluck("data", 10) %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# liver vs skin
pca_all %>% 
  pluck("data", 11) %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# liver vs spleen
pca_all %>% 
  pluck("data", 12) %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Heatmap
# Fat vs skin
counts_preprocessed %>% 
  pluck("data", 8) %>% 
  keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>% 
  tidyHeatmap::heatmap(symbol, sample, count_scaled, transform = log1p) %>% 
  add_tile(tissue) %>% 
  add_tile(mouse)

  
  

# Old preliminary analysis =============================

# preprocessing
counts_preprocessed <- counts_annotated %>% 
  aggregate_duplicates() %>% 
  keep_abundant(factor_of_interest = tissue) %>% 
  scale_abundance()

# generate contrast
tissue_contrast <- counts_preprocessed %>% 
  distinct(tissue) %>% 
  mutate(tissue = paste0("tissue", tissue)) %>%
  mutate(background = "tissueSpleen") %>% 
  mutate(contrast = sprintf("%s - %s", tissue, background)) %>% 
  filter(tissue != background) %>% 
  pull(contrast)

# DE analysis
de_all <- counts_preprocessed %>% 
  test_differential_abundance(
    ~ 0 + tissue + mouse,
    .abundance = count_scaled,
    .contrasts = tissue_contrast
  )

x <- de_all %>% 
  pivot_longer(contains("___"), 
               names_to = c("stats", "contrast"), 
               names_sep = "___",
               values_to = "value") %>% 
  nest(stat_df = - c(tissue, contrast)) %>% 
  mutate(stat_df = map(
    stat_df, 
    ~ .x %>% 
      pivot_wider(names_from = stats, values_from = value)
  ))

# Spleen vs VAT
counts_splvat <- counts_annotated %>% 
  select(ID, symbol, contains(c("tisTregST2_Spleen", "tisTregST2_Fat"))) %>%
  pivot_longer(-c(ID, symbol), names_to = "sample", values_to = "count") %>% 
  mutate(tissue = str_extract(sample, "(?<=_).+(?=_)")) %>% 
  mutate(mouse = str_extract(sample, "R\\d+")) %>% 
  tidybulk(sample, symbol, count)

# compare the total reads between samples
counts_splvat %>% 
  group_by(sample) %>% 
  summarise(total_counts = sum(count))

# check that each gene is only present in a sample once so that there are no duplicates
counts_splvat %>% 
  count(symbol, sample) %>% 
  arrange(desc(n))

counts_splvat %>% 
  filter(symbol == "A230057D06Rik" & sample == "tisTregST2_Fat_R1")

# preprocess: aggregate_duplicates with the sample symbol, filter low transcripts, normalisation
counts_scaled_splvat <- counts_splvat %>% 
  aggregate_duplicates() %>% 
  keep_abundant(factor_of_interest = tissue) %>% 
  scale_abundance()

# check data quality after scaling
counts_scaled_splvat %>% 
  pivot_longer(contains("count"), names_to = "source", values_to = "abundance") %>% 
  ggplot(aes(x = abundance+1, colour = sample)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~source) +
  custom_theme

# PCA
counts_scaled_splvat %>% 
  reduce_dimensions(sample, symbol, count_scaled, method = "PCA") %>% 
  pivot_sample() %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Heatmap
counts_scaled_splvat %>% 
  keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>% 
  tidyHeatmap::heatmap(symbol, sample, count_scaled, transform = log1p) %>% 
  add_tile(tissue) %>% 
  add_tile(mouse)

# Differential Gene Expression
de_splvat <- counts_scaled_splvat %>% 
  test_differential_abundance(
    ~ tissue + mouse,
    .abundance = count_scaled
  )

top_genes_splvat <- de_splvat %>% 
  pivot_transcript() %>% 
  slice_min(PValue, n = 30) %>% 
  pull(symbol)


de_splvat %>%
  pivot_transcript() %>%
  
  # Subset data
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_splvat, symbol, "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 20) +
  
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "Spleen vs VAT") +
  theme(plot.title = element_text(hjust = 0.5))

# Spleen vs Lung
counts_spllung <- counts_annotated %>% 
  select(ID, symbol, contains(c("tisTregST2_Spleen", "tisTregST2_Lung"))) %>%
  pivot_longer(-c(ID, symbol), names_to = "sample", values_to = "count") %>% 
  mutate(tissue = str_extract(sample, "(?<=_).+(?=_)")) %>% 
  mutate(mouse = str_extract(sample, "R\\d+")) %>% 
  tidybulk(sample, symbol, count)

# compare the total reads between samples
counts_spllung %>% 
  group_by(sample) %>% 
  summarise(total_counts = sum(count))

# check that each gene is only present in a sample once so that there are no duplicates
counts_spllung %>% 
  count(symbol, sample) %>% 
  arrange(desc(n))

counts_spllung %>% 
  filter(symbol == "A230057D06Rik" & sample == "tisTregST2_Lung_R1")

# preprocess: aggregate_duplicates with the sample symbol, filter low transcripts, normalisation
counts_scaled_spllung <- counts_spllung %>% 
  aggregate_duplicates() %>% 
  keep_abundant(factor_of_interest = tissue) %>% 
  scale_abundance()

# check data quality after scaling
counts_scaled_spllung %>% 
  pivot_longer(contains("count"), names_to = "source", values_to = "abundance") %>% 
  ggplot(aes(x = abundance+1, colour = sample)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~source) +
  custom_theme

# PCA
counts_scaled_spllung %>% 
  reduce_dimensions(sample, symbol, count_scaled, method = "PCA") %>% 
  pivot_sample() %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Heatmap
counts_scaled_spllung %>% 
  keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>% 
  tidyHeatmap::heatmap(symbol, sample, count_scaled, transform = log1p) %>% 
  add_tile(tissue) %>% 
  add_tile(mouse)

# Differential Gene Expression
de_spllung <- counts_scaled_spllung %>% 
  test_differential_abundance(
    ~ tissue + mouse,
    .abundance = count_scaled
  )

top_genes_spllung <- de_spllung %>% 
  pivot_transcript() %>% 
  slice_min(PValue, n = 30) %>% 
  pull(symbol)


de_spllung %>%
  pivot_transcript() %>%
  # Subset data
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_spllung, symbol, "")) %>%
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 30) +
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "Spleen vs Lung") +
  theme(plot.title = element_text(hjust = 0.5))

# Spleen vs Skin
counts_splskin <- counts_annotated %>% 
  select(ID, symbol, contains(c("tisTregST2_Spleen", "tisTregST2_Skin"))) %>%
  pivot_longer(-c(ID, symbol), names_to = "sample", values_to = "count") %>% 
  mutate(tissue = str_extract(sample, "(?<=_).+(?=_)")) %>% 
  mutate(mouse = str_extract(sample, "R\\d+")) %>% 
  tidybulk(sample, symbol, count)

# compare the total reads between samples
counts_splskin %>% 
  group_by(sample) %>% 
  summarise(total_counts = sum(count))

# check that each gene is only present in a sample once so that there are no duplicates
counts_splskin %>% 
  count(symbol, sample) %>% 
  arrange(desc(n))

counts_splskin %>% 
  filter(symbol == "A230057D06Rik" & sample == "tisTregST2_Skin_R1")

# preprocess: aggregate_duplicates with the sample symbol, filter low transcripts, normalisation
counts_scaled_splskin <- counts_splskin %>% 
  aggregate_duplicates() %>% 
  keep_abundant(factor_of_interest = tissue) %>% 
  scale_abundance()

# check data quality after scaling
counts_scaled_splskin %>% 
  pivot_longer(contains("count"), names_to = "source", values_to = "abundance") %>% 
  ggplot(aes(x = abundance+1, colour = sample)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~source) +
  custom_theme

# PCA
counts_scaled_splskin %>% 
  reduce_dimensions(sample, symbol, count_scaled, method = "PCA") %>% 
  pivot_sample() %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Heatmap
counts_scaled_splskin %>% 
  keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>% 
  tidyHeatmap::heatmap(symbol, sample, count_scaled, transform = log1p) %>% 
  add_tile(tissue) %>% 
  add_tile(mouse)

# Differential Gene Expression
de_splskin <- counts_scaled_splskin %>% 
  test_differential_abundance(
    ~ tissue + mouse,
    .abundance = count_scaled
  )

top_genes_splskin <- de_splskin %>% 
  pivot_transcript() %>% 
  slice_min(PValue, n = 30) %>% 
  pull(symbol)


de_splskin %>%
  pivot_transcript() %>%
  # Subset data
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_splskin, symbol, "")) %>%
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 20) +
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "Spleen vs Skin") +
  theme(plot.title = element_text(hjust = 0.5))

# # Spleen vs Liver
counts_splliver <- counts_annotated %>% 
  select(ID, symbol, contains(c("tisTregST2_Spleen", "tisTregST2_Liver"))) %>%
  pivot_longer(-c(ID, symbol), names_to = "sample", values_to = "count") %>% 
  mutate(tissue = str_extract(sample, "(?<=_).+(?=_)")) %>% 
  mutate(mouse = str_extract(sample, "R\\d+")) %>% 
  tidybulk(sample, symbol, count)

# compare the total reads between samples
counts_splliver %>% 
  group_by(sample) %>% 
  summarise(total_counts = sum(count))

# check that each gene is only present in a sample once so that there are no duplicates
counts_splliver %>% 
  count(symbol, sample) %>% 
  arrange(desc(n))

counts_splliver %>% 
  filter(symbol == "A230057D06Rik" & sample == "tisTregST2_Fat_R1")

# preprocess: aggregate_duplicates with the sample symbol, filter low transcripts, normalisation
counts_scaled_splliver <- counts_splliver %>% 
  aggregate_duplicates() %>% 
  keep_abundant(factor_of_interest = tissue) %>% 
  scale_abundance()

# check data quality after scaling
counts_scaled_splliver %>% 
  pivot_longer(contains("count"), names_to = "source", values_to = "abundance") %>% 
  ggplot(aes(x = abundance+1, colour = sample)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~source) +
  custom_theme

# PCA
counts_scaled_splliver %>% 
  reduce_dimensions(sample, symbol, count_scaled, method = "PCA") %>% 
  pivot_sample() %>% 
  ggplot(aes(PC1, PC2, color = tissue, shape = mouse)) +
  geom_point() +
  custom_theme

# Heatmap
counts_scaled_splliver %>% 
  keep_variable(.abundance = count_scaled, top = 100) %>%
  as_tibble() %>% 
  tidyHeatmap::heatmap(symbol, sample, count_scaled, transform = log1p) %>% 
  add_tile(tissue) %>% 
  add_tile(mouse)

# Differential Gene Expression
de_splliver <- counts_scaled_splliver %>% 
  test_differential_abundance(
    ~ tissue + mouse,
    .abundance = count_scaled
  )

top_genes_splliver <- de_splliver %>% 
  pivot_transcript() %>% 
  slice_min(PValue, n = 30) %>% 
  pull(symbol)


de_splliver %>%
  pivot_transcript() %>%
  # Subset data
  mutate(significant = FDR < 0.05 & abs(logFC) >= 2) %>%
  mutate(symbol = ifelse(symbol %in% top_genes_splliver, symbol, "")) %>%
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = symbol)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel(max.overlaps = 20) +
  # Custom scales
  custom_theme +
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  labs(title = "Spleen vs Liver") +
  theme(plot.title = element_text(hjust = 0.5))

