rm(list = ls())

library(tidyverse)
library(DESeq2)

mycounts<-read.csv("gene_count.csv")
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
rownames(mycounts)<-mycounts[,1]
mycounts<-mycounts[,-1]
head(mycounts)
condition <- factor(c(rep("control",2),rep("treat",2)), levels = c("control","treat"))
condition
colData <- data.frame(row.names = colnames(mycounts), condition)
colData
#以上，得到两张表，mycounts,以及cloData,用于进行DEseq2分析

dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds

res = results(dds)
res = res[order(res$pvalue),]
head(res)
summary(res)
write.csv(res,file="All_results.csv")
#以上，得到全部基因按照p值大小排列的表达差异表

table(res$padj<0.05)
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
dim(diff_gene_deseq2)
head(diff_gene_deseq2)
write.csv(diff_gene_deseq2,file= "DEG_treat_vs_control.csv")
#最终得到显著差异表达基因矩阵，存储为DEG_treat_vs_control.csv。

RefSeq_ID <-row.names(diff_gene_deseq2)
RefSeq_ID
NM_ID <- tibble(RefSeq_ID)
write.csv(NM_ID,file = 'NM_ID.csv')
#暂时未找到biomaRt直接转换RefSeqID的函数，故将RefSeqID导出，在网页处理，然后以ensemble_ID.csv导入
#转换ID的时候发现存在对应不上和重复对应的情况，将这些数据删除后，进行了下面的biomaRt注释

ensemble_ID <- read.csv("ensemble_ID.csv")
head(ensemble_ID)
rownames(ensemble_ID)<-ensemble_ID[,1]
ensemble_ID<-ensemble_ID[,-1]
head(ensemble_ID)
rownames(ensemble_ID)<-ensemble_ID[,1]
head(ensemble_ID)
my_ensembl_gene_id <- row.names(ensemble_ID)
my_ensembl_gene_id

library(biomaRt)
library(curl)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
options(timeout = 4000000)
homo_symbols<- getBM(attributes=c('ensembl_gene_id','external_gene_name',"description"),
                     filters = 'ensembl_gene_id', values = my_ensembl_gene_id, mart = mart)
head(homo_symbols)
write.csv(homo_symbols,file="homo_symbols.csv")

#目前得到了了两张表，一张是根据原始数据的Refseq_ID作的差异分析，一张是对得到的差异
#分析结果，将Refseq_ID先转换为ensemble_ID,再利用biomaRt进行了注释，现在问题是，在Refseq_ID换
#到ensemble_ID时，删除了十一行信息，这些行要么找不到对应的ensemble_ID,要么存在多个RefseqID对
#应同一个ensemble_ID的现象。所以下一步如何合并呢














