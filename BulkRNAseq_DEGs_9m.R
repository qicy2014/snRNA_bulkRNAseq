library(DESeq2)
data <- read.table("WT_PS19_9_rawcount.txt", header=TRUE, quote="\t", row.names=1)
countData <- as.matrix(data)
rownames(countData) <- rownames(data)
database <- data.frame(name=colnames(data), condition=c("WT","WT","WT","WT","WT","PS19","PS19","PS19","PS19","PS19"))
rownames(database) <-colnames(data)
dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
dds <- dds[rowSums(counts(dds)) > 5,]
dds <- estimateSizeFactors(dds)
dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","PS19","WT"))
write.table(res, "PS19-WT-9m-DESeq2-rerun_output.txt",sep="\t")


