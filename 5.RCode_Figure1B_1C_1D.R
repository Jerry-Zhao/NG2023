setwd("/Users/sohyunmoon/Desktop/Collaboration_Work/NJ/Project4_Abdel_RNAseq_NG_S961/Code_Manuscript")

library("DESeq2")


data_A <- read.table("Merged_Abdel2023_NG_3Ctrl_3Treat_exonic_gene_Jerry.csv", sep=",", header=T, row.names=1)
group1=rep(c("Ctrl","Treat"), c(3,3))
row1 <- colnames(data_A)
label1=c("Ctrl","Treat")

df <- data.frame(data_A)
df1 <- data.frame(condition=group1, row.names=row1, check.names=FALSE)
dds <- DESeqDataSetFromMatrix(countData=df, colData=df1, design=~condition)
keep <- rowSums(counts(dds) >= 10) >= 3 # 3 replicate ## at least 5 samples with >=10 counts; 5 = least replicate number
dds <- dds[keep,]
dds <- DESeq(dds) # Wald significance tests

res <- results(dds,contrast=c("condition","Treat","Ctrl"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
sum(resOrdered$padj < 0.05, na.rm=TRUE)
print(dim(as.data.frame(resOrdered))) # 16129     6

total <- resOrdered$padj
print(sum(total < 0.05, na.rm=TRUE)) ## Total DEG 90   
print(sum(total < 0.05 & resOrdered$log2FoldChange>0, na.rm=TRUE)) # Up 66
print(sum(total < 0.05 & resOrdered$log2FoldChange<0, na.rm=TRUE)) # Down 24
write.table(as.data.frame(resOrdered), file="DESeq2_Wald_NG_Treat_Ctrl_exonic_gene_results.xls", sep="\t", quote=F)


######## Normalized counts and CPMs ######### 
## CPM values
cpm_counts <- round(fpm(dds, robust = TRUE),3) # default, median of ratios method
write.table(cpm_counts, file="DESeq2_CPM_NG_Treat_Ctrl_exonic_gene.xls", sep="\t", quote=F)


####################################
####### Figure 1B PCA plot #########
####################################

### PCA plot
rld <- rlog(dds)
print(plotPCA(rld, intgroup=c("condition"),ntop = 1000))

result2 <- plotPCA(rld, intgroup=c("condition"),ntop = 1000, returnData = TRUE)
CtrlPCA <- result2[1:3,]
TreatPCA <- result2[4:6,]

png("Fig1B_DESeq2_PCA_NG_Treat_Ctrl_exonic_gene_top1000.png", width = 300, height = 300)
plot(CtrlPCA$PC1, CtrlPCA$PC2,  xlim=c(-25,25),ylim=c(-10,10),pch=21, cex=2, lwd=2, bg="cyan3", col="black", xlab="Principal component 1 (67%)", ylab="Principal component 2 (20%)")
points(TreatPCA$PC1, TreatPCA$PC2, pch=21,cex=2,lwd=2, bg="brown1", col="black")
dev.off()

result2
#                   PC1        PC2 group condition       name
# Ctrl_rep1   -1.043116  6.1768352  Ctrl      Ctrl  Ctrl_rep1
# Ctrl_rep2  -10.181005  0.2105473  Ctrl      Ctrl  Ctrl_rep2
# Ctrl_rep3    2.362852  8.7152258  Ctrl      Ctrl  Ctrl_rep3
# Treat_rep1  -7.210748 -3.9911287 Treat     Treat Treat_rep1
# Treat_rep2  -5.743235 -7.3158825 Treat     Treat Treat_rep2
# Treat_rep3  21.815252 -3.7955971 Treat     Treat Treat_rep3


####################################
###### Figure 1C Volcano plot ######
####################################

total <- read.table("DESeq2_Wald_NG_Treat_Ctrl_exonic_gene_results.xls", sep="\t", header=T, row.names=NULL)
total <- na.omit(total) # Remove NA
total <- total[,c(1,3,7)] # 1:ID, 3: log2FC, 7:FDR
dim(total) # 16123     3

names(total)[1] <- "genes"
names(total)[2] <- "logFC"
names(total)[3] <- "FDR"
head(total)

# cutoff 0.05
fdr <- total[total$FDR < 0.05,]
up <- fdr[fdr$logFC > 0,]
down <- fdr[fdr$logFC < 0,]

print(dim(fdr)) # 90   
print(dim(up)) # 66	
print(dim(down)) # 24  

 
pdf("Fig1C_Volcano.pdf")
plot(total$logFC,-log(total$FDR,10),col="grey", pch=21, xlab="log2FC",ylab="-log10FDR",xlim=c(-10, 10),ylim=c(0, 30))
points(up$logFC, -log(up$FDR,10), pch=21,  cex=1, bg="brown1", lwd=0.8,col="black")
points(down$logFC, -log(down$FDR,10),  pch=21, cex=1, bg="cyan3", lwd=0.8,col="black")
dev.off()



####################################
#### Figure 1D Z-scoreHeatmap  ####
####################################

library("pheatmap")
library("RColorBrewer")

cpm_total <- fpm(dds)

# Call CPM data in FDR<0.05 (90 genes)
deg_cpm <- cpm_total[c("ENSMUSG00000058400", "ENSMUSG00000050359", "ENSMUSG00000026247", "ENSMUSG00000099839", "ENSMUSG00000031722", "ENSMUSG00000042816", "ENSMUSG00000076434", "ENSMUSG00000027737", "ENSMUSG00000037578", "ENSMUSG00000031074", "ENSMUSG00000006014", "ENSMUSG00000067235", "ENSMUSG00000026628", "ENSMUSG00000115016", "ENSMUSG00000051985", "ENSMUSG00000037095", "ENSMUSG00000037035", "ENSMUSG00000005373", "ENSMUSG00000069793", "ENSMUSG00000020592", "ENSMUSG00000019890", "ENSMUSG00000063632", "ENSMUSG00000002897", "ENSMUSG00000030898", "ENSMUSG00000001473", "ENSMUSG00000024366", "ENSMUSG00000036251", "ENSMUSG00000025473", "ENSMUSG00000019647", "ENSMUSG00000023903", "ENSMUSG00000073418", "ENSMUSG00000036390", "ENSMUSG00000048490", "ENSMUSG00000024907", "ENSMUSG00000005087", "ENSMUSG00000079466", "ENSMUSG00000021453", "ENSMUSG00000052684", "ENSMUSG00000023067", "ENSMUSG00000059991", "ENSMUSG00000015243", "ENSMUSG00000046186", "ENSMUSG00000000303", "ENSMUSG00000059456", "ENSMUSG00000039037", "ENSMUSG00000008393", "ENSMUSG00000022952", "ENSMUSG00000025348", "ENSMUSG00000051379", "ENSMUSG00000040118", "ENSMUSG00000033350", "ENSMUSG00000029608", "ENSMUSG00000041594", "ENSMUSG00000026773", "ENSMUSG00000031543", "ENSMUSG00000035673", "ENSMUSG00000027858", "ENSMUSG00000046480", "ENSMUSG00000037428", "ENSMUSG00000022212", "ENSMUSG00000034040", "ENSMUSG00000020277", "ENSMUSG00000030095", "ENSMUSG00000021703", "ENSMUSG00000028832", "ENSMUSG00000037907", "ENSMUSG00000036760", "ENSMUSG00000058624", "ENSMUSG00000098650", "ENSMUSG00000040147", "ENSMUSG00000051650", "ENSMUSG00000034755", "ENSMUSG00000056427", "ENSMUSG00000045954", "ENSMUSG00000020182", "ENSMUSG00000022665", "ENSMUSG00000005803", "ENSMUSG00000024232", "ENSMUSG00000020241", "ENSMUSG00000023094", "ENSMUSG00000030638", "ENSMUSG00000024027", "ENSMUSG00000037994", "ENSMUSG00000021032", "ENSMUSG00000067786", "ENSMUSG00000022376", "ENSMUSG00000039488", "ENSMUSG00000061816", "ENSMUSG00000075394", "ENSMUSG00000021886"),]

dim(deg_cpm)
deg_cpm[1:5,]
deg_zscore <- t(scale(t(deg_cpm)))

pdf("Fig1D_heatmap_zscore.pdf")
pheatmap(deg_zscore, color = colorRampPalette(rev(brewer.pal(n =5, name ="RdYlBu")))(500), breaks = seq(-2, 2, length.out = 500), cluster_cols = FALSE, cluster_rows=FALSE)
dev.off()


