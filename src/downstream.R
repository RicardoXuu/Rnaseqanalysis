rm(list = ls())

options(stringsAsFactors = FALSE)
control1 <- read.table("SRR1210370.count",sep = "\t",col.names = c("gene_id","control1"))
control2 <- read.table("SRR1210373.count",sep = "\t",col.names = c("gene_id","control2"))
treat1 <- read.table("SRR1210368.count",sep = "\t",col.names = c("gene_id","treat1"))
treat2 <- read.table("SRR1210372.count",sep = "\t",col.names = c("gene_id","treat2"))
head(control1)
tail(control1)

raw_count <- merge(merge(control1,control2,by="gene_id"),merge(treat1,treat2,by="gene_id"))
head(raw_count)
tail(raw_count)
raw_count_filt <- raw_count[-1:-5,]
head(raw_count_filt)

write.csv(raw_count_filt,file = 'read_count.csv')
read_count <- read.csv(file = 'read_count.csv')
gene_count <- read_count[,2:6]
write.csv(gene_count,file = 'gene_count.csv')
#最终得到一个合并后并删除掉无关信息的gene_count.csv文件，用于下一步进行差异分析





