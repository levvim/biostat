# Levim Jul 2017 Wolchok Lab
# Deseq2 analysis for differential expression
# featurecounts import and formatting for deseq2 and agnostic graphs taken from steven turner: https://gist.github.com/stephenturner/f60c1934405c127f09a6
# featurecounts input by default assumes that the first line is trashed (featurecounts command input) and counts are in the 7th column
# This example utilizes 8 samples, treatments A vs F, with 4 replicates in each group (5,6,7,8)

################################################################################
library("DESeq2")
library("ggplot2")
library("reshape2")
library("gplots")
library("gridExtra")
library("genefilter")
library("PoiClaClu")
library("pheatmap")
library("RColorBrewer")
library("AnnotationDbi")
library("org.Hs.eg.db") #for human
library("org.Mm.eg.db") #for mouse
library(data.table)

#################################################################################
# Navigate to folder containing featurecounts output
setwd("~/git/biostat/deseq2/counts/")

# Specify sample names and their corresponding counts files
samples=c(A5 = "A5.counts.txt",
          A6 = "A6.counts.txt",
          A7 = "A7.counts.txt",
          A8 = "A8.counts.txt",
          F5 = "F5.counts.txt",
          F6 = "F6.counts.txt",
          F9 = "F9.counts.txt",
          F8 = "F8.counts.txt")

# Specify experimental design
# This can be adapted do the design, but essentially requires a chr vector converted to factors for sample grouping
# The design is defaulted to `~condition+sample`. to change this edit the DESeqDataSetFromMatrix() fn

(condition=c("A","A","A","A","F","F","F","F"))
(sample=c("5","6","7","8","5","6","7","8"))

# can also be done using rep
(condition=c(rep("A", 4), rep("F",4)))
(sample=c(rep(c("5","6","7","8"),2)))

################################################################################
# Import multiple samples and create counts df

# loop across samples and create a list of dataframes. using `str(dfl)`, confirm that colnames match the content (specifically Gene and Counts)
dfl=list()
for (i in 1:length(samples)) { 
        dfl[[i]]=read.csv(samples[i], sep="\t", stringsAsFactors=FALSE, skip=1, header=TRUE, check.names=FALSE, quote="")
        colnames(dfl[[i]]) = c("Gene", "Chr", "Start", "End", "Strand", "Length", "Counts")
}

str(dfl)

# create counts dataframe with corresponding sample names
counts=data.frame(dfl[[1]][,c("Counts")]); 
rownames(counts)=dfl[[1]]$Gene; colnames(counts)=names(samples[1])

for (i in 2:length(samples)) { 
        counts[,names(samples[i])] = dfl[[i]]$Counts
}

head(counts)

################################################################################
# Set up experimental design
### Whole experiment

head(counts)

### Assign design variables
(condition = factor(as.character(condition)))
(sample= factor(as.character(sample)))

# navigate to an output directory if necessary
#setwd("~/git/biostat/deseq2/")
#################################################################################
#Analysis with DESeq2

# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata = data.frame(row.names=colnames(counts), cbind.data.frame(condition, sample)))
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~condition+sample)

# for single factors
#(coldata = data.frame(row.names=colnames(counts), cbind.data.frame(sample)))
#dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~sample)

#prefilter dataset
nrow(dds)
dds = dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

#rlog transformation
rld = rlog(dds, blind=FALSE)
head(assay(rld), 3)
dds = estimateSizeFactors(dds)

pdf(file = "count_scatterplot.pdf")
par( mfrow = c( 1, 2 ) )
plot(log2(counts(dds, normalized=TRUE)[,1:2] + 1), main = "raw", pch=16, cex=0.3)
plot(assay(rld)[,1:2], main = "rlog transformed", pch=16, cex=0.3)
dev.off()

#sample distances
sampleDists = dist( t( assay(rld) ) )
sampleDists

dds = estimateSizeFactors(dds)
sampleDistMatrix = as.matrix( sampleDists )
rownames(sampleDistMatrix) = paste( rld$condition, rld$sample, sep="-" )
colnames(sampleDistMatrix) = NULL
colors = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#heatmap of sample to sample distance using rlog transformed values
poisd = PoissonDistance(t(counts(dds)))
samplePoisDistMatrix = as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) = paste( rld$condition, rld$sample, sep="-" )
colnames(samplePoisDistMatrix) = NULL

pdf(file = "sample_distance_heatmaps.pdf")
par(mfrow=c(1,2))
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

pheatmap(samplePoisDistMatrix, clustering_distance_rows=poisd$dd, clustering_distance_cols=poisd$dd, col=colors)
dev.off()

#pca
pdf(file = "pca_raw.pdf")
plotPCA(rld, intgroup = c("condition", "sample"))
dev.off()

#pca using rlog transformed values
data = plotPCA(rld, intgroup = c("condition", "sample"), returnData=TRUE)

pdf(file = "pca.pdf")
percentVar = round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=sample)) + geom_point(size=3) +
          xlab(paste0("PC1: ",percentVar[1],"% variance")) +
          ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

#mds plot
pdf("mds_raw.pdf")
mdsData = data.frame(cmdscale(sampleDistMatrix))
mds = cbind(mdsData, as.data.frame(colData(rld)))
ggplot(mds, aes(X1,X2,color=condition,shape=sample)) + geom_point(size=3)
dev.off()

#mds using rlog transformed values
pdf("mds_rlog.pdf")
mdsPoisData = data.frame(cmdscale(samplePoisDistMatrix))
mdsPois = cbind(mdsPoisData, as.data.frame(colData(dds)))
ggplot(mdsPois, aes(X1,X2,color=condition,shape=sample)) + geom_point(size=3)
dev.off()

################################################################################
#differential expression
#lower the false discovery rate threshold (the threshold on padj in the results table)
#change fdr threshold via:
#res.05 <- results(dds, alpha=.05)
#table(res.05$padj < .05)
#raise the log2 fold change threshold from 0 using the lfcThreshold argument of results
#raise log2 fold change threshold via
#resLFC1 <- results(dds, lfcThreshold=1)
#table(resLFC1$padj < 0.1)

dds = DESeq(dds)

(res = results(dds))
summary(res)
mcols(res, use.names=TRUE) #structure of results table

#if ran twice using BioMart: annotations
#rownames(res) = make.names(res$symbol, unique = TRUE)
#rownames(res)=rownames(counts)

#to compare specific levels of a variable:
results(dds, contrast=c("sample", "ctl", "exp"))

#multiple testing of FDR = .1 (using BH correction)
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))


resSig = subset(res, padj < 0.05)
head(resSig[ order(resSig$log2FoldChange), ]) #strongest down regulation
head(resSig[ order(resSig$log2FoldChange, decreasing=TRUE), ])#strongest upregulation

#visualizing results per gene cohort
topGene = rownames(res)[which.min(res$padj)]
head(topGene)

#for results regarding single genes:
##normalized counts for single gene over treatment group
#pdf("rlog_topGene_v_treatment.pdf")
#plotCounts(dds, gene=topGene, intgroup=c("condition"))
#dev.off()
##normalized counts indicating cell line with color
#data = plotCounts(dds, gene=topGene, intgroup=c("condtion","sample"), returnData=TRUE)
#ggplot(data, aes(x=condition, y=count, color=sample)) +
#          scale_y_log10() + 
#            geom_point(position=position_jitter(width=.1,height=0), size=3)
##normalized counts using a more structural arrangement
#ggplot(data, aes(x=condition y=count, fill=condition)) +
#          scale_y_log10() + 
#            geom_dotplot(binaxis="y", stackdir="center")
##normalized counts with lines connecting celll lines
#ggplot(data, aes(x=condition, y=count, color=sample, group=sample)) +
#          scale_y_log10() + geom_point(size=3) + geom_line()

resLFC1 <- results(dds, lfcThreshold=1)

#ma plot with changes induced by treatment
pdf("ma_changesByTreatment.pdf")
plotMA(res, ylim=c(-5,5))
dev.off()

#An MA-plot of a test for large log2 fold changes. 
pdf("ma_largeLog2FoldChange.pdf")
plotMA(resLFC1, ylim=c(-5,5))
topGene = rownames(resLFC1)[which.min(resLFC1$padj)]
with(resLFC1[topGene, ], {
          points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
          text(baseMean, log2FoldChange, topGene, pos=2, col="dodgerblue")
          })
dev.off()

#histogram of p vals for gene with mean normalized countt larger than 1
pdf("hist_genePval_meanNormalizedCount_largerThan1.pdf")
hist(res$pvalue[res$baseMean > 1], breaks=0:20/20, col="grey50", border="white")
dev.off()

#gene clustering
topVarGenes = head(order(rowVars(assay(rld)),decreasing=TRUE),20)
mat = assay(rld)[ topVarGenes, ]
mat = mat - rowMeans(mat)
df = as.data.frame(colData(rld)[,c("sample","condition")])
pdf("heatmap_relRlogVals_samples.pdf")
pheatmap(mat, annotation_col=df)
dev.off()

#independent filtering: the ratio of small p values for genes binned by mean normalized count
qs <- c(0, quantile(resLFC1$baseMean[resLFC1$baseMean > 0], 0:6/6))
bins <- cut(resLFC1$baseMean, qs)
levels(bins) <- paste0("~",round(signif(.5*qs[-1] + .5*qs[-length(qs)],2)))
ratios <- tapply(resLFC1$pvalue, bins, function(p) mean(p < .05, na.rm=TRUE))
pdf("ratio_smallpVal_genesBinnedByMeanRlogCount.pdf")
barplot(ratios, xlab="mean normalized count", ylab="ratio of small p values")
dev.off()

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
# Annotate output if necessary

##This is the organism annotation package (“org”) for Homo sapiens (“Hs”), organized as an AnnotationDbi database package (“db”), using Entrez Gene IDs (“eg”) as primary key. To get a list of all available key types, use columns("org.Hs.eg.db")
#res$symbol = mapIds(org.Mm.eg.db,
#    keys=row.names(res),
#    column="SYMBOL",
#    keytype="ENSEMBL",
#    multiVals="first")
#res$entrez = mapIds(org.Mm.eg.db,
#    keys=row.names(res),
#    column="ENTREZID",
#    keytype="ENSEMBL",
#    multiVals="first")

resOrdered = res[order(res$padj),]
head(resOrdered)

# Export table/results
resOrderedDF <- as.data.frame(resOrdered)
write.csv(resOrderedDF, file="results.csv")

# Get differential expression results and export:
table(res$padj<0.05)

# Order by adjusted p-value
# Merge with normalized count data
resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
head(resdata)

write.csv(resdata, file="diffexpr_results.csv")

#################################################################################

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
