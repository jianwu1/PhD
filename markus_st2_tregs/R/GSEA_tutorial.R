library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(ReactomePA)
library(fgsea)

# install_github("ctlab/fgsea")

packageVersion("msigdbr")

# GO analyses (groupGO(), enrichGO() and gseGO()) enricher() buildGOmap()
data(geneList, package="DOSE")
# notice that gene used for groupGO and enrichGO are filtered genes with logFC > 2
gene <- names(geneList)[abs(geneList)>2]

ggo <- groupGO(
  gene = gene, # a vector of gene IDs (can be any ID type that is supported by the corresponding OrgDb)
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "CC", # one of "MF", "BP", "CC"
  level = 3, # specific GO level
  readable = TRUE # if readable is TRUE, gene IDs will map to gene symbols
)

# GO Over-representation analysis
ego <- enrichGO(
  gene = gene,
  universe = names(geneList),
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE
)

gene.df <- bitr(
  gene,
  fromType = "ENTREZID",
  toType = c("ENSEMBL", "SYMBOL"),
  OrgDb = org.Hs.eg.db
  # drop = TRUE by default (drop NAs)
)

ego2 <- enrichGO(
  gene.df$ENSEMBL,
  universe = names(geneList) %>% bitr(fromType = "ENTREZID", toType = "ENSEMBL", org.Hs.eg.db) %>% pull(ENSEMBL),
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

# Gene set enrichment analysis: GSE is the general option, there are specific versions: gseGO(), gseKEGG
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "CC",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)



