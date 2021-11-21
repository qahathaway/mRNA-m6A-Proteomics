
####Salmon/Tximport####

library(GenomicFeatures)
library(tximport)
library(readr)
library(rjson)
library(dplyr)
library(AnnotationDbi)

###Importing transacript abundances and construction of a gene-level DESeqDataSet object from Salmon quants files###
library(tximportData)
dir <- system.file("extdata", package = "tximportData", mustWork=TRUE)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
coldata

coldata <- coldata[1:6,]
coldata$names <- coldata$SampleName
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf")
file.exists(coldata$files)

library(tximeta)
se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
dim(gse)

library(SummarizedExperiment)
data(gse)

gse
assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
round( colSums(assay(gse)) / 1e6, 1 )
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)

###DifferentialExpressionAnalysis###

library("DESeq2")
ddsTxi <- DESeqDataSet(gse, design = ~ Exposure)
ddsTxi

#keep <- rowSums(counts(ddsTxi) >= 100) >= 2 #(group has to have >5 counts)
#dds <- ddsTxi[keep,]

ddsTxi$Exposure <- factor(ddsTxi$Exposure, levels = c("Sham","TiO2")) #The text, condition treated vs untreated, tells you that the estimates are of the logarithmic fold change log2(treated/untreated).
#dds$Exposure <- relevel(dds$Exposure, ref = "Sham")

dds <- DESeq(ddsTxi)
res <- results(dds, name="Exposure_TiO2_vs_Sham")
res
summary(res)


###Log fold change shrinkage for visualization and ranking###

resultsNames(dds)
resLFC <- lfcShrink(dds, coef="Exposure_TiO2_vs_Sham", type="apeglm")
resLFC
resOrdered <- res[order(res$pvalue),] #Orders results table by smallest p-value#
summary(res)
sum(res$padj < 0.05, na.rm=TRUE) #Tells us how many adj p-values were less than 0.05


##Exporting Results##
resOrdered <- res[order(res$padj),]
head(resOrdered)
resOrderedDF <- as.data.frame(resOrdered)[1:3647, ]
write.csv(resOrderedDF, file = "path/to/file.csv")


###Exploring Results###

##MA-plot shows the log2 fold changes attributable to a given varaible over the mean of normalized counts for all the samples in the DeseqDataSet##
##Points will be colored red if the adjusted p value is less than 0.1. Points which fall out of the window are plotted as open triangles pointing either up or down## 
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2)) #MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds#
idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

#MA-plotWithTopGene#
plotMA(resLFC, ylim = c(-2,2))
topGene <- rownames(resLFC)[which.min(resLFC$padj)]
with(resLFC[topGene, ], {
  points(baseMean, log2FoldChange, col="dodgerblue", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="red")
})

#Alternative Shrinkage Estimators#
#resultsNames(dds)
# because we are interested in treated vs untreated, we set 'coef=2'
#resNorm <- lfcShrink(dds, coef=2, type="normal")
#resAsh <- lfcShrink(dds, coef=2, type="ashr")


par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
plotMA(resLFC, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAsh, xlim=xlim, ylim=ylim, main="ashr")


##Plot Counts##
#Examining the counts of reads for a single gene across the groups#
#plotCounts(dds, gene=which.min(res$padj), intgroup="Exposure")

##GPX4##
GeneofChoice <- "ENSMUSG00000075706"
plotCounts(dds, gene = GeneofChoice, intgroup=c("Exposure"), returnData = TRUE)
plotCounts(dds, gene = GeneofChoice, intgroup=c("Exposure"), returnData = FALSE)

#Customized Plotting#
# an argument returnData specifies that the function should only return a data.frame for plotting with ggplot#
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Exposure", 
                returnData=TRUE) #Here we specify the gene which had the smallest p value from the results table created above#
library("ggplot2")
ggplot(d, aes(x=Exposure, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Exposure"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Exposure, y = count)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

ggplot(geneCounts, aes(x = Exposure, y = count)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()


##MoreInfoOnResultsColums##
#Information about which variables and tests were used can be found by calling the function mcols on the results object.#
mcols(res)$description

###Data Normalization Assessment###

##Extracting transformed values##
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)

colData(vsd)

rld <- rlog(dds, blind = FALSE)
head(assay(rld), 3)

##Effects of transformations (vst and rlog) on the varaince##
# this gives log2(n + 1)
ntd <- normTransform(dds)
library("vsn")
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

##Removing batch effects##
#mat <- assay(vsd)
#mm <- model.matrix(~Exposure, colData(vsd))
#mat <- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
#assay(vsd) <- mat
#plotPCA(vsd)

##Comparing the two transformations##
library("dplyr")
library("ggplot2")

dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df$transformation <- factor(df$transformation, levels=lvls)

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)


##Sample Distances##
sampleDists <- dist(t(assay(vsd)))
sampleDists

#Heatmap of sample-to-sample distances using the vsd-transformed values#
library("pheatmap")
library("RColorBrewer")
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$Exposure, sep = " - " )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors)

#Heatmap of sample-to-sample distances using the Poisson Distance: measure of dissimilarity between counts that also takes the inherent variance strcutrue of counts into consideration#
library("PoiClaClu")
poisd <- PoissonDistance(t(counts(dds))) #Takes original count matrix (not normalized as it normalizes internally) with samples as rows instead of columns so we transpose the counts in dds)#

samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$SampleName, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)


#PCA Plot-another way to visualize sample-to-sample distances (using the vsd-transformed values)#
plotPCA(vsd, intgroup = c("Exposure"))

#Customize PCA Plot#
pcaData <- plotPCA(vsd, intgroup = c("Exposure"), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = Exposure)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")


#Generalzied version of PCA for dimesnsion reduction of non-normally distributed data such as counts and binary matrices#
library("glmpca")
gpca <- glmpca(counts(dds), L=2)
gpca.dat <- gpca$factors
gpca.dat$Exposure <- dds$Exposure

ggplot(gpca.dat, aes(x = dim1, y = dim2, color = Exposure)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA")


##MDS Plot)similar to PCA plot but usig multidimensional scaling function-useful when we don't have a matrix of data##
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = Exposure)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")

mdsPois <- as.data.frame(colData(dds)) %>%
  cbind(cmdscale(samplePoisDistMatrix))
ggplot(mdsPois, aes(x = `1`, y = `2`, color = Exposure)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with PoissonDistances")


###Data Quality Assessment by Sample Clustering and Visualization###

library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
               decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Exposure")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

cdata <- colData(dds)
pheatmap(assay(ntd),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = as.data.frame(cdata[,"Exposure"], row.names=rownames(cdata)))

pheatmap(assay(vsd),
         cluster_rows = FALSE,
         show_rownames = FALSE,
         cluster_cols = FALSE,
         annotation_col = as.data.frame(cdata[,"Exposure"], row.names=rownames(cdata)))


##Gene Clustering##

library("genefilter")
#topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 500)
#or#
topGenes <- head(order(res$padj),decreasing = TRUE, 1000)

mat  <- assay(vsd)[ topGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("Exposure")])
colnames(anno) <- c('Exposure')
rownames(anno) <- colnames(mat)

#my_colour = list(
#  cell = c(Four = "#8B4513", One = "#708090"),
#  dex = c(trt = "#FF0000", trt2 = "#F08080", trt3 = "#FFFF00", trt4 = "#FFD700",
#          trt5 = "#2E8B57", trt6 = "#00FF00", untrt = "#00BFFF", untrt2 = "#1E90FF"))

pheatmap(mat, annotation_col = anno, show_rownames = FALSE)


##Other Visualizations##
library(vidger)

vsBoxPlot(dds, d.factor = "Exposure", type = "deseq")
vsDEGMatrix(dds, d.factor = "Exposure", padj = 0.05, type = "deseq")
#vsFourWay("trt2", "trt", "untrt2", d.factor = "Exposure", dds, type = "deseq")
vsMAMatrix(dds, d.factor = "Exposure", type = "deseq", y.lim = c(-10,10))
vsScatterMatrix(dds, d.factor = "Exposure", type = "deseq")
#vsScatterPlot("untrt2", "trt2", dds, d.factor = "group", type = "deseq")
#vsVolcano("untrt2", "trt2", dds, d.factor = "group", type = "deseq", x.lim = c(-10,10), padj = 0.05)
vsVolcanoMatrix(dds, d.factor = "Exposure", type = "deseq", lfc = 2, padj = 0.05, x.lim = c(-8,8),
                title = FALSE, legend = TRUE, grid = TRUE, counts = FALSE, facet.title.size = 10)


##Volcano Plots##
library(EnhancedVolcano)
EnhancedVolcano(res, lab = rownames(res), x = "log2FoldChange", y = "padj", pCutoff = 0.05,
                xlim = c(-4.5, 4.5), ylim = c(0, 4))



##Histogram of p vaules for genes with mean normalized count larger than 2##
hist(res$pvalue[res$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")
