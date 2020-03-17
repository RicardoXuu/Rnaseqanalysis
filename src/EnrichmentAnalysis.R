
library(clusterProfiler)
library(DOSE)
library(org.Hs.eg.db)
library(stringr)

#首先，得到差异基因的ENTREZID，在之前存储的ensemble_ID.csv文件中有
sig.gene <- read.csv(file="ensemble_ID.csv")
head(sig.gene)
gene <- sig.gene[,2]
head(gene)
gene.df <- bitr(gene, fromType = "ENSEMBL", 
                      toType = c("SYMBOL","ENTREZID"),
                      OrgDb = org.Hs.eg.db)
head(gene.df)

#GO分析，GO:CC是细胞组分富集分析，GO:BP是生物过程富集分析
ego_cc <- enrichGO(gene       = gene.df$ENSEMBL,
                   OrgDb      = org.Hs.eg.db,
                   keyType    = 'ENSEMBL',
                   ont        = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
ego_bp <- enrichGO(gene       = gene.df$ENSEMBL,
                   OrgDb      = org.Hs.eg.db,
                   keyType    = 'ENSEMBL',
                   ont        = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.01,
                   qvalueCutoff = 0.05)
barplot(ego_bp,showCategory = 18,title="The GO_BP enrichment analysis of all DEGs ")+ 
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(ego_bp) str_wrap(ego_bp,width = 50))


#KEGG富集分析
kk <- enrichKEGG(gene = gene.df$ENTREZID,
                 organism = 'hsa',
                 pvalueCutoff = 0.05)
kk[1:30]
barplot(kk,showCategory = 30, title="The KEGG enrichment analysis of all DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(kk) str_wrap(kk,width = 30))
dotplot(kk,showCategory = 30, title="The KEGG enrichment analysis of all DEGs")+
  scale_size(range=c(2, 12))+
  scale_x_discrete(labels=function(kk) str_wrap(kk,width = 30))


# Gene Set Enrichment Analysis（GSEA）
# 从之前存储的差异基因表GKD.csv，获取按照log2FC大小来排序的基因列表
gkd <- read.csv(file = "GKD.csv")
gkd <- gkd[,-1]
head(gkd)
genelist <- gkd$log2FoldChange
names(genelist) <- gkd[,1]
head(genelist)
genelist <- sort(genelist, decreasing = TRUE)
head(genelist)
# 根据前面得到的基因列表，进行GSEA分析
gsemf <- gseGO(genelist,
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               ont="BP"
)
head(gsemf)
str(gsemf)
# 从GSEA的结果中，找出一个，画出GSEA图，目前没有看懂这个图的含义，后期学习
gseaplot(gsemf, geneSetID="GO:0006810")
gseaplot(gsemf, geneSetID="GO:0000302")

#以上，完成差异基因的GO，KEGG，GSEA富集分析


















