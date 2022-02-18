library(msigdbr)
library(fgsea)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(pathview)
library(org.Mm.eg.db)
library(ReactomePA)
library(tidyverse)
library(GO.db)
library(biomartr)


# Functions
ORA_GO <- function(.de_genes, .bg_genes, .row_num){
  .de_genes %>% 
    pluck("data", .row_num) %>% 
    filter(logFC>0) %>% 
    pull(ID) %>% 
    enrichGO(
      universe = .bg_genes %>% pluck("data", .row_num) %>% pull(ID),
      OrgDb = org.Mm.eg.db,
      keyType = "ENSEMBL",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05
    ) %>% 
    setReadable(OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
}

GSEA_GO <- function(.bg_genes, .row_num, ...){
  .bg_genes %>% 
    pluck("data", .row_num) %>% 
    arrange(desc(logFC)) %>% 
    select(ID, logFC) %>% 
    deframe() %>% 
    gseGO(
      ont = "BP",
      OrgDb = org.Mm.eg.db,
      keyType = "ENSEMBL",
      pvalueCutoff = 0.01,
      eps = 0,
      pAdjustMethod = "BH",
      ...
    ) %>% 
    setReadable(OrgDb = org.Mm.eg.db, keyType = "ENSEMBL")
}

ORA_KEGG <- function(.de_genes, .bg_genes, .row_num){
  .de_genes %>% 
    pluck("data", .row_num) %>% 
    filter(logFC>0) %>% 
    pull(ID) %>% 
    bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>% 
    drop_na() %>% 
    pull(ENTREZID) %>% 
    enrichKEGG(
      organism = "mmu",
      universe = .bg_genes %>% 
        pluck("data", .row_num) %>% 
        pull(ID) %>% 
        bitr(fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db) %>% 
        drop_na() %>% 
        pull(ENTREZID),
      keyType = "kegg",
      pAdjustMethod = "BH",
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05
    ) %>% 
    setReadable(OrgDb = org.Mm.eg.db,
                keyType = "ENTREZID")
}

GSEA_KEGG <- function(.bg_genes, .row_num, ...){
  .bg_genes %>% 
    pluck("data", .row_num) %>% 
    arrange(desc(logFC)) %>% 
    mutate(ENTREZID = AnnotationDbi::mapIds(
      org.Mm.eg.db,
      keys = ID,
      keytype = "ENSEMBL",
      column = "ENTREZID"
    )) %>% 
    drop_na() %>% 
    select(ENTREZID, logFC) %>% 
    deframe() %>% 
    gseKEGG(
      organism = "mmu",
      keyType = "kegg",
      pvalueCutoff = 0.01,
      pAdjustMethod = "BH",
      eps = 0,
      ...
    ) %>% 
    setReadable(OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
}

ORA_MSigDb <- function(.de_genes, .row_num, .gene_set){
  .de_genes %>% 
    pluck("data", .row_num) %>% 
    filter(logFC>0) %>% 
    pull(ID) %>% 
    enricher(
      pvalueCutoff = 0.01,
      qvalueCutoff = 0.05,
      TERM2GENE = .gene_set %>% select(gs_name, ensembl_gene),
      TERM2NAME = .gene_set %>% select(gs_name, gs_description)
    ) %>% 
    setReadable("org.Mm.eg.db", "ENSEMBL")
}

GSEA_MSigDb <- function(.bg_genes, .row_num, .gene_set, ...){
  .bg_genes %>% 
    pluck("data", .row_num) %>% 
    arrange(desc(logFC)) %>% 
    select(ID, logFC) %>% 
    deframe() %>% 
    GSEA(
      pvalueCutoff = 0.01,
      pAdjustMethod = "BH",
      eps = 0,
      TERM2GENE = .gene_set %>% select(gs_name, ensembl_gene),
      TERM2NAME = .gene_set %>% select(gs_name, gs_description),
      ...
    ) %>% 
    setReadable(OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL")
}

# VAT vs Spleen
# GO enrichment analysis
## GO annotation
go_vatspl <- de_all_filtered %>% 
  pluck("data", 9) %>% 
  filter(logFC > 0) %>% 
  pull(symbol) %>% 
  groupGO(
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",
    level = 4
  )

View(go_vatspl@result)

# de_all_filtered %>% 
#   pluck("data", 9) %>% 
#   arrange(desc(logFC)) %>% 
#   pull(ID) %>% 
#   getGO(
#     organism = "Mus musculus",
#     genes = .,
#     filters = "ensembl_gene_id"
#         )
# x %>% 
#   as_tibble() %>% 
#   mutate(terms = goslim_goa_accession %>% 
#            mget(GOTERM, ifnotfound = NA)) %>% 
#   unnest(terms)
# a <- x$goslim_goa_accession %>% 
#   mget(GOTERM, ifnotfound = NA)
  
sapply(sample(GOTERM, 10), GOID) 
sapply(sample(GOTERM, 10), Term) %>% as.data.frame()
sapply(sample(GOTERM, 10), Synonym)
sapply(sample(GOTERM, 10), Secondary)
sapply(sample(GOTERM, 10), Definition)
sapply(sample(GOTERM, 10), Ontology)

GOTerm = tibble(
  GOID = sapply(GOTERM, GOID),
  Term = sapply(GOTERM, Term),
  Synonym = sapply(GOTERM, Synonym),
  Secondary = sapply(GOTERM, Secondary),
  Definition = sapply(GOTERM, Definition),
  Ontology = sapply(GOTERM, Ontology)
  ) %>% 
  mutate(Synonym = map_chr(Synonym, ~ paste(.x, collapse = "/"))) %>% 
  mutate(Secondary = map_chr(Secondary, ~ifelse(length(.x)==0, NA, .x)))

job::job({saveRDS(GOTerm, "intermediate_data/GOTerm.rds", compress = "xz")})

## Over-representation analysis
ORA_GO_vatspl <- ORA_GO(de_all_filtered, de_all, 9)
ORA_GO_skinspl <- ORA_GO(de_all_filtered, de_all, 15)
ORA_GO_vatskin <- ORA_GO(de_all_filtered, de_all, 8)


ORA_GO@result %>% 
  filter(qvalue<0.05) %>% 
  View()

## GSEA analysis
GSEA_GO_vatspl <- GSEA_GO(de_all, 9)
GSEA_GO_skinspl <- GSEA_GO(de_all, 15, nPermSimple = 1000000)
GSEA_GO_vatskin <- GSEA_GO(de_all, 8, nPermSimple = 100000)

GSEA_GO_vatspl@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues < 0.05) %>% 
  View()

GSEA_GO_vatskin@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues < 0.05) %>% 
  View()

GSEA_GO_skinspl@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues < 0.05) %>% 
  View()

goplot(GSEA_GO_vatskin)
gseaplot2(GSEA_GO)

# KEGG enrichment analysis
search_kegg_organism('Mus musculus', by='scientific_name')

ORA_KEGG_vatspl <- ORA_KEGG(de_all_filtered, de_all, 9)
ORA_KEGG_skinspl <- ORA_KEGG(de_all_filtered, de_all, 15)
ORA_KEGG_vatskin <- ORA_KEGG(de_all_filtered, de_all, 8)

ORA_KEGG_vatspl@result %>% 
  filter(qvalue<0.05) %>% 
  View()

ORA_KEGG_skinspl@result %>% 
  filter(qvalue<0.05) %>% 
  View()

ORA_KEGG_vatskin@result %>% 
  filter(qvalue<0.05) %>% 
  View()


# KEGG GSEA analysis
GSEA_KEGG_vatspl <- GSEA_KEGG(de_all, 9)
GSEA_KEGG_skinspl <- GSEA_KEGG(de_all, 15, nPermSimple = 100000)
GSEA_KEGG_vatskin <- GSEA_KEGG(de_all, 8, nPermSimple = 10000)

GSEA_KEGG_vatspl@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues<0.05) %>% 
  View()

GSEA_KEGG_skinspl@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues<0.05) %>% 
  View()

GSEA_KEGG_vatskin@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues<0.05) %>% 
  View()

browseKEGG(GSEA_KEGG, "mmu04915")
pathview(de_all %>% 
           pluck("data", 9) %>% 
           select(symbol, logFC) %>% 
           deframe(),
         pathway.id = "mmu04915", 
         species = "mmu",
         limit = list(gene=(de_all %>% 
                        pluck("data", 9) %>% 
                        select(symbol, logFC) %>% 
                        deframe() %>% 
                        abs() %>% 
                        max()),
                      cpd = 1
                        )
         )
  
# KEGG module gene set enrichment analysis
ORA_KEGGM <- de_all_filtered %>% 
  pluck("data", 9) %>% 
  filter(logFC>0) %>%
  mutate(ENTREZID = AnnotationDbi::mapIds(org.Mm.eg.db, 
                                          keys = symbol, 
                                          keytype = "SYMBOL", 
                                          column = "ENTREZID")
  ) %>% 
  drop_na(ENTREZID) %>% 
  pull(ENTREZID) %>% 
  enrichMKEGG(
    organism = "mmu",
    universe = de_all %>% 
      pluck("data", 9) %>% 
      pull(symbol) %>% 
      bitr(fromType = "SYMBOL", toType = "ENTREZID", org.Mm.eg.db),
    keyType = "kegg",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.01,
    qvalueCutoff = 0.05
  ) %>% 
  setReadable(OrgDb = org.Mm.eg.db,
              keyType = "ENTREZID")

ORA_KEGGM@result %>% 
  filter(qvalue<0.05)

# Reactome

# Universal enrichment analysis
# enricher() for hypergeometric test
# GSEA() for gene set enrichment analysis
# TERM2GENE: .df: term ID | mapped gene
# TERM2NAME (optional) is a .df: term ID | term name 

# Molecular Signature Database
# H: hallmark gene sets
# C1: positional gene sets
# C2: curated gene sets
# C3: motif gene sets
# C4: computational gene sets
# C5: GO gene sets
# C6: oncogenic signatures
# C7: immunologic signatures

# users can download GMT files from Broad Institute and use the read.gmt() function
# to parse the file to be used in enricher() and GSEA()
# alternatively msigdbr already packed the MSigDB gene sets in tidydata format that can be
# used directly with clusterProfiler

msigdbr_show_species()

# retrieve mouse gene sets of interest
m_H <- msigdbr("Mus musculus", "H")
m_C2 <- msigdbr("Mus musculus", "C2")
m_C5 <- msigdbr("Mus musculus", "C5")
m_C7 <- msigdbr("Mus musculus", "C7")

# MSigDb over-representation analysis
ORA_MSigDb_H_vatspl <- ORA_MSigDb(de_all_filtered, 9, m_H)
ORA_MSigDb_H_skinspl <- ORA_MSigDb(de_all_filtered, 15, m_H)
ORA_MSigDb_H_vatskin <- ORA_MSigDb(de_all_filtered, 8, m_H)

ORA_MSigDb_C2_vatspl <- ORA_MSigDb(de_all_filtered, 9, m_C2)
ORA_MSigDb_C2_skinspl <- ORA_MSigDb(de_all_filtered, 15, m_C2)
ORA_MSigDb_C2_vatskin <- ORA_MSigDb(de_all_filtered, 8, m_C2)

ORA_MSigDb_C5_vatspl <- ORA_MSigDb(de_all_filtered, 9, m_C5)
ORA_MSigDb_C5_skinspl <- ORA_MSigDb(de_all_filtered, 15, m_C5)
ORA_MSigDb_C5_vatskin <- ORA_MSigDb(de_all_filtered, 8, m_C5)

ORA_MSigDb_C7_vatspl <- ORA_MSigDb(de_all_filtered, 9, m_C7)
ORA_MSigDb_C7_skinspl <- ORA_MSigDb(de_all_filtered, 15, m_C7)
ORA_MSigDb_C7_vatskin <- ORA_MSigDb(de_all_filtered, 8, m_C7)

ORA_MSigDb_H_vatspl@result %>% 
  filter(qvalue<0.05) %>% 
  View()

ORA_MSigDb_H_skinspl@result %>% 
  filter(qvalue<0.05) %>% 
  View()

ORA_MSigDb_H_vatskin@result %>% 
  filter(qvalue<0.05) %>% 
  View()

# MSigDb GSEA analysis
GSEA_MSigDb_H_vatspl <- GSEA_MSigDb(de_all, 9, m_H)
GSEA_MSigDb_H_skinspl <- GSEA_MSigDb(de_all, 15, m_H, nPermSimple = 10000)
GSEA_MSigDb_H_vatskin <- GSEA_MSigDb(de_all, 8, m_H)

GSEA_MSigDb_C2_vatspl <- GSEA_MSigDb(de_all, 9, m_C2)
GSEA_MSigDb_C2_skinspl <- GSEA_MSigDb(de_all, 15, m_C2, nPermSimple = 100000)
GSEA_MSigDb_C2_vatskin <- GSEA_MSigDb(de_all, 8, m_C2, nPermSimple = 10000)

GSEA_MSigDb_C5_vatspl <- GSEA_MSigDb(de_all, 9, m_C5)
GSEA_MSigDb_C5_skinspl <- GSEA_MSigDb(de_all, 15, m_C5, nPermSimple = 10000)
GSEA_MSigDb_C5_vatskin <- GSEA_MSigDb(de_all, 8, m_C5, nPermSimple = 10000)

GSEA_MSigDb_C7_vatspl <- GSEA_MSigDb(de_all, 9, m_C7)
GSEA_MSigDb_C7_skinspl <- GSEA_MSigDb(de_all, 15, m_C7, nPermSimple = 100000)
GSEA_MSigDb_C7_vatskin <- GSEA_MSigDb(de_all, 8, m_C7, nPermSimple = 10000)

GSEA_MSigDb_H_vatspl@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues<0.05) %>% 
  View()

GSEA_MSigDb_H_skinspl@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues<0.05) %>% 
  View()

GSEA_MSigDb_H_vatskin@result %>% 
  arrange(desc(NES)) %>% 
  filter(qvalues<0.05) %>% 
  View()


# Visualisation

# barplot (only available for enrichment results)
barplot_vatspl <- barplot(ORA_MSigDb_H_vatspl, 
        x = "GeneRatio",
        color = "qvalue", 
        showCategory = 20,
        title = "Over-representated pathways for genes enriched in VAT vs Spleen"
        ) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/Over-representated pathways for genes enriched in VAT vs Spleen.png", barplot_vatspl)

barplot_skinspl <- barplot(ORA_MSigDb_H_skinspl, 
        x = "GeneRatio",
        color = "qvalue", 
        showCategory = 10,
        title = "Over-representated pathways for genes enriched in Skin vs Spleen"
        ) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/Over-representated pathways for genes enriched in Skin vs Spleen.png", barplot_skinspl)

# dot plot analysis (works for both ORA and GSEA and compareCluster)
dotplot_vatspl <- GSEA_MSigDb_H_vatspl %>% 
  rename(qvalue = qvalues) %>% 
  dotplot(x = "GeneRatio",
          color = "qvalue",
          showCategory = 15,
          font.size = 10,
          title = "GSEA for ST2+ Tregs in VAT vs Spleen"
          ) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/GSEA for ST2+ Tregs in VAT vs Spleen.png", dotplot_vatspl)

dotplot_skinspl <- GSEA_MSigDb_H_skinspl %>% 
  rename(qvalue = qvalues) %>% 
  dotplot(x = "GeneRatio",
          color = "qvalue",
          showCategory = 15,
          font.size = 10,
          title = "GSEA for ST2+ Tregs in Skin vs Spleen"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("results/GSEA for ST2+ Tregs in VAT vs Spleen.png", dotplot_vatspl)

# Gene Concept Network
GSEA_MSigDb_H_vatspl %>% 
  rename(qvalue = qvalues) %>% 
  cnetplot(
    showCategory = 5,
    foldChange = de_all %>% pluck("data", 9) %>% select(symbol, logFC) %>% arrange(desc(logFC)) %>% deframe() %>% head(),
    categorySize = "qvalue"
  ) +
  labs(title = "Core enriched genes in top5 enriched pathways (ST2+Treg in VAT vs Spleen)",
    size = "size by qvalue") +
  theme(plot.title = element_text(hjust = 0.5))

GSEA_MSigDb_H_skinspl %>% 
  rename(qvalue = qvalues) %>% 
  cnetplot(
    showCategory = 5,
    foldChange = de_all %>% pluck("data", 15) %>% select(symbol, logFC) %>% arrange(desc(logFC)) %>% deframe(),
    categorySize = "qvalue"
  ) +
  labs(title = "Core enriched genes in top5 enriched pathways (ST2+Treg in Skin vs Spleen)",
       size = "size by qvalue") +
  theme(plot.title = element_text(hjust = 0.5))

# Tree plot
GSEA_MSigDb_H_vatspl %>% 
  rename(qvalue = qvalues) %>% 
  pairwise_termsim() %>% 
  treeplot(
    color = "qvalue"
  ) +
  labs(
    title = "Hierarchical clustering of enriched pathways (ST2+Treg in VAT vs Spleen)"
       ) +
  theme(plot.title = element_text(hjust = 0.5))

GSEA_MSigDb_H_skinspl %>% 
  rename(qvalue = qvalues) %>% 
  pairwise_termsim() %>% 
  treeplot(
    color = "qvalue"
  ) +
  labs(
    title = "Hierarchical clustering of enriched pathways (ST2+Treg in Skin vs Spleen)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# GSEA plot
gseaplot2(GSEA_MSigDb_H_vatspl, geneSetID = 1, title = GSEA_MSigDb_H_vatspl$Description[1]) +
  theme(plot.title = element_text(hjust = 0.5))

gseaplot2(GSEA_MSigDb_H_vatspl, geneSetID = 1:3, title = "GSEA_MSigDb_H_vatspl") +
  theme(plot.title = element_text(hjust = 0.5))


# Biological theme comparison



# creating dataframe
de_filtered_vat_skin_spl <- 
  de_all_filtered %>% 
  slice(9, 8, 15) %>% 
  unnest(data)


GOID = biomartr::getGO(organism = "Mus musculus", genes = de_filtered_vat_skin_spl$ID, filters = "ensembl_gene_id")


output <- GOID %>% 
  as_tibble() %>% 
  left_join(GOTerm, by = c("goslim_goa_accession" = "GOID")) %>% 
  select(-goslim_goa_description) %>% 
  rename(GOID = goslim_goa_accession) %>% 
  nest(GO = - ensembl_gene_id) %>% 
  left_join(de_filtered_vat_skin_spl, ., by = c("ID" = "ensembl_gene_id")) %>% 
  
  mutate(description = map_chr(
    GO,
    ~ if(!is.null(.x)){
      .x %>% 
        pull(Term) %>% 
        unique() %>% 
        paste(collapse = "/")
    }else{NA}
  )) %>% 
  select(contrast, ENSEMBL=ID, symbol, description, merged_transcripts, .abundant, logFC, logCPM, `F`, PValue, FDR) %>% 
  nest(data = -contrast)

saveRDS(output, "results/annotated_de_filtered_spl_vat_skin.rds", compress = "xz")

write.xlsx2(
  output %>% pluck("data", 1),
  file = "results/ST2_Treg_DE.xlsx",
  sheetName = "Fat_vs_Spleen",
  append = FALSE
)

write.xlsx2(
  output %>% pluck("data", 2),
  file = "results/ST2_Treg_DE.xlsx",
  sheetName = "Fat_vs_Skin",
  append = TRUE
)

write.xlsx2(
  output %>% pluck("data", 3),
  file = "results/ST2_Treg_DE.xlsx",
  sheetName = "Skin_vs_Spleen",
  append = TRUE
)

#+++++++++++++++++++++++++++
# xlsx.writeMultipleData
#+++++++++++++++++++++++++++++
# file : the path to the output file
# ... : a list of data to write to the workbook
xlsx.writeMultipleData <- function (file, ...)
{
  require(xlsx, quietly = TRUE)
  objects <- list(...)
  fargs <- as.list(match.call(expand.dots = TRUE))
  objnames <- as.character(fargs)[-c(1, 2)]
  nobjects <- length(objects)
  for (i in 1:nobjects) {
    if (i == 1)
      write.xlsx(objects[[i]], file, sheetName = objnames[i])
    else write.xlsx(objects[[i]], file, sheetName = objnames[i],
                    append = TRUE)
  }
}

