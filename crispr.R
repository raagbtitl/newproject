
#R
dt<-read.table("crispr_guide_data.tsv",sep="\t",header=T,stringsAsFactors = FALSE,)
dt1 <- apply( dt[ , colnames(dt)[1:3] ] , 1 , paste , collapse = "-" )
dt<-cbind(dt1,dt)
rownames(dt)<-dt[,1];dt<-dt[,-1]
dt<-dt[,-c(1:3)]
col_name<-colnames(dt)
a<-sapply(strsplit(col_name,"_"),"[", 3)
pdf("variability.pdf")
boxplot(log2(v),notch=TRUE,col=(c("gold","darkgreen")),las=2)
dev.off()
conds=factor(a)
pdf("pca.pdf")
samples.pca <- prcomp(t(log10(dt+1)), scale. = TRUE) # Performs principal component analysis after scaling the data.
plot(samples.pca$x[,1], samples.pca$x[,2],col=conds,xlab="PC1",ylab="PC2", pch=19, cex=2)
text(samples.pca$x[,1], samples.pca$x[,2], 1:dim(v)[2], cex=0.8, pos=2, col="red")
legend("top", levels(conds), pch=19, cex=1, col=1:length(levels(conds)))
grid()
dev.off()
require(pheatmap)
library("RColorBrewer")
sampleDists <- dist( t( log2(v+1) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
x11();pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,
               col=colors)
## compare d14 versus plx14
chose<-c(4,5,8,9)
dt_newy<-log2(dt[, chose])
cond <- as.factor(a[chose])
design <- model.matrix(~cond+0)
colnames(design) <- levels(cond)
require(limma)
fit <- lmFit(dt_newy, design)
cont.matrix <- makeContrasts(PLX14-D14, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit,0.01)
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", n=nrow(fit2$t))
write.table(tT2, file="crispr_plasmid_d14_vs_plx14.csv", sep="\t")
cn<-read.table("CCLE_copynumber_byGene_2013-12-03.txt",sep="\t",header=T,row.names=1)
ind<-grep("A375", colnames(cn))
cn_a375<-cn[, c(1:4,ind)]
gain<- (cn_a375[which(cn_a375[,5]>0.3),])
loss<- cn_a375[which(cn_a375[,5]< -0.3),]
seg.dat.fn<-read.table("CCLE_copynumber_2013-12-03.seg.a375.txt",sep="\t",header=T,stringsAsFactors = FALSE)
sif<-read.table("CCLE_SNP.Arrays.sif_2013-12-03.txt",sep="\t",header=T,row.names=1,stringsAsFactors = FALSE)




