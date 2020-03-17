
head(diff_gene_deseq2)
str(diff_gene_deseq2)
dim(diff_gene_deseq2)
diff_gene_deseq2_design <- diff_gene_deseq2
head(diff_gene_deseq2_design)
ID_tran <- read.csv("ID_tran.csv")
head(ID_tran)
ID_tran <- ID_tran[,-1]
head(ID_tran)
#这里复制了diff_gene_deseq2,并创建表ID_tran，用这两个表完成refseqID向ensembleID的转换

RefSeq_ID <- rownames(diff_gene_deseq2_design)
diff_gene_deseq2_design <- cbind(RefSeq_ID,diff_gene_deseq2_design)
head(diff_gene_deseq2_design)
colnames(diff_gene_deseq2_design)[1]<-c("RefSeq_ID")
head(diff_gene_deseq2_design)
diff_name <- merge(ID_tran,diff_gene_deseq2_design,by="RefSeq_ID")
table(diff_name$padj<0.05)
diff_name = diff_name[order(diff_name$pvalue),]
head(diff_name)
dim(diff_name)
diff_name <- diff_name[-353,] 
diff_name <- diff_name[-334,]
diff_name <- diff_name[-323,] 
diff_name <- diff_name[-227,] 
diff_name <- diff_name[-225,] 
diff_name <- diff_name[-221,] 
diff_name <- diff_name[-177,] 
diff_name <- diff_name[-176,] 
diff_name <- diff_name[-151,] 
diff_name <- diff_name[-93,] 
diff_name <- diff_name[-67,] 
dim(diff_name)
head(diff_name)
diff_name <- diff_name[,-1]
head(diff_name)
write.csv(diff_name,file = 'diff_name.csv')
#现在解决了两个表数据行不对等的问题，得到新表diff_name，已经完成了ID转换

GKD <- merge(diff_name,homo_symbols,by="ensembl_gene_id")
head(GKD)
GKD = GKD[order(GKD$pvalue),]
head(GKD)

CXCL8 <- GKD[GKD$external_gene_name == "CXCL8",]
head(CXCL8)
write.csv(GKD,file = 'GKD.csv')
#以上，对新表以及基因注释表进行合并，得到并保存文件GKD，作为最终的差异基因表





