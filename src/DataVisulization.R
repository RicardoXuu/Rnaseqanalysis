
library(DESeq2)
library(ggplot2)
library(dplyr)

plotMA(res,ylim=c(-2,2))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=6, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

res_order<-res[order(row.names(res)),]
res = res_order
res.shrink <- lfcShrink(dds, contrast = c("condition","treat","control"), res=res)
plotMA(res.shrink, ylim = c(-5,5))
topGene <- rownames(res)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
})

# 不画图，只显示数据
plotCounts(dds, gene=which.min(res$padj), intgroup="condition", returnData=TRUE)
#只画图，不显示数据
plotCounts(dds, gene="NM_006536", intgroup="condition", returnData=FALSE)

plotCounts(dds, gene="NM_006536", intgroup="condition", returnData=TRUE) %>%
    ggplot(aes(condition, count)) + geom_boxplot(aes(fill=condition)) + scale_y_log10() + ggtitle("CLCA2")

d <- plotCounts(dds, gene="NM_006536", intgroup="condition", returnData=TRUE)
ggplot(d, aes(x=condition, y=count)) + 
  geom_point(aes(color= condition),size= 4, position=position_jitter(w=0.5,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))+ ggtitle("CLCA2")

vsdata <- vst(dds, blind=FALSE)
plotPCA(vsdata, intgroup="condition")

library("pheatmap")
select<-order(rowMeans(counts(dds, normalized = TRUE)),
              decreasing = TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])
ntd <- normTransform(dds)
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsdata)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsdata$condition, vsdata$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


