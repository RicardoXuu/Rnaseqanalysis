rm(list = ls())


library("pathview")
library("gage")
library("gageData")
library("dplyr")
library("clusterProfiler")
library("DOSE")
library("stringr")
library("org.Hs.eg.db")


data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs =  kegg.sets.hs[sigmet.idx.hs]
head(kegg.sets.hs,3)

sig.gene <- read.csv(file="ensemble_ID.csv")
head(sig.gene)
gene <- sig.gene[,2]
head(gene)
gene.df<-bitr(gene, fromType = "ENSEMBL", 
                    toType = c("SYMBOL","ENTREZID"),
                    OrgDb = org.Hs.eg.db)
head(gene.df)

gkd <- read.csv(file = "GKD.csv")
gkd <- gkd[,-1]
head(gkd)
foldchanges <- gkd$log2FoldChange
names(foldchanges)= gene.df$ENTREZID
head(foldchanges)

keggres = gage(foldchanges, gsets = kegg.sets.hs, 
                            same.dir = TRUE)
# Look at both up (greater), down (less), and statatistics.
lapply(keggres, head)

keggrespathways = data.frame(id=rownames(keggres$greater), keggres$greater) %>% 
  tbl_df() %>% 
  filter(row_number()<=10) %>% 
  .$id %>% 
  as.character()
keggrespathways

# Get the IDs.
keggresids = substr(keggrespathways, start=1, stop=8)
keggresids

# 下面开始画KEGG通路，先定义画图函数
plot_pathway = function(pid) pathview(gene.data=foldchanges, pathway.id=pid, 
                                      species="hsa", new.signature=FALSE)
# 然后同时画多个pathways，这些plots自动存到工作目录
tmp = sapply(keggresids, function(pid) pathview(gene.data=foldchanges, pathway.id=pid, 
                                                species="hsa"))














