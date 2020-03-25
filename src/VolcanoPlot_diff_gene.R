rm(list = ls())

#导入差异分析全表，并筛选存在Padj值的数据行
data_Volcano <- read.csv("All_results.csv")
table(data_Volcano$padj<0.05)
data_Volcano <-subset(data_Volcano,padj > 0)

library(ggplot2)
library(ggrepel)

#对数据集进行处理，标记出显著下调和显著上调的基因
data_Volcano$threshold = factor(
  ifelse(data_Volcano$padj < 0.05 & abs(data_Volcano$log2FoldChange) >= 1, 
         ifelse(data_Volcano$log2FoldChange>= 1 ,'Up','Down'),'NoSignifi'),
  levels=c('Up','Down','NoSignifi')
  )

#画出火山图
ggplot(data_Volcano,aes(x=log2FoldChange,y=-log2(padj),color=threshold))+
  geom_point()+
  scale_color_manual(values=c("#DC143C","#00008B","#808080"))+      #确定点的颜色
  geom_text_repel(
    data = data_Volcano[data_Volcano$padj<0.000000000000000000000000000000000000000005,],
    aes(label = X),
    size = 3,
    segment.color = "black", show.legend = FALSE )+                 #添加关注的点的基因名
  theme_bw()+                                                       #修改图片背景
  theme(
    legend.title = element_blank()                                  #不显示图例标题
  )+
  ylab('-log2 (p-adj)')+                                            #修改y轴名称
  xlab('log2 (FoldChange)')+                                        #修改x轴名称
  geom_vline(xintercept=c(-1,1),lty=3,col="black",lwd=0.5) +        #添加横线|FoldChange|>2
  geom_hline(yintercept = -log2(0.05),lty=3,col="black",lwd=0.5)    #添加竖线padj<0.05

