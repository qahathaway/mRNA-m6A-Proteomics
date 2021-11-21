
###FULL CODE###

library(RADAR)
samplenames <- c( paste0("FSham_",1:3), paste0("FExp_",1:3) )

bam_dir <- system.file('extdata',package = 'm6A')
samplenames

radar <- countReads(samplenames = samplenames,
                    gtf = "path/to/file.gtf",
                    bamFolder = bam_dir,
                    modification = "m6A",
                    fragmentLength = 297,
                    outputDir = "path/to/file/reads",
                    saveOutput = T,
                    threads = 32,
                    binSize = 25
)

summary(radar)


radar <- normalizeLibrary(radar)
sizeFactors(radar)
radar <- adjustExprLevel(radar)
variable(radar) <- data.frame( Group = c(rep("FS_",3),rep("FE_",3)) )
radar <- filterBins(radar,minCountsCutOff = 5)


radar <- diffIP_parallel(radar, thread = 32)


top_bins <- extractIP(radar,filtered = T)[order(rowMeans( extractIP(radar,filtered = T) ),decreasing = T)[1:10000],]
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )


radar <- reportResult(radar, cutoff = 0.1, Beta_cutoff = 0.5, thread = 32)


result <- results(radar)
head(result)

write.csv(result, file = "path/to/file.csv")

radar <- PrepCoveragePlot(radar, gtf = "path/to/file.gtf")
plotHeatMap(radar)
peakDistribution(radar)

RNASeqCounts <- geneExpression(radar)

IPreadsn <- extractIP(radar, normalized = TRUE)
IPreadsa <- extractIP(radar, adjusted = TRUE)
IPreadsf <- extractIP(radar, filtered = TRUE)
write.csv(IPreadsn, file = "path/to/file/IP_Reads_Normalized.csv")
write.csv(IPreadsa, file = "path/to/file/IP_Reads_Adjusted.csv")
write.csv(IPreadsf, file = "path/to/file/IP_Reads_Filtered.csv")

Inputreads <- extractInput(radar, normalized = TRUE)
write.csv(Inputreads, file = "path/to/file/Input_Reads_Normalized.csv")

##Cnot1##
plotGeneCov(radar,geneName = "ENSMUSG00000036550", center = median, libraryType = "opposite", adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000036550", center = median, libraryType = "opposite",ZoomIn = c(96457229,96457253 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000036550", center = median, libraryType = "opposite",ZoomIn = c(96457229,96457253 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000036550", libraryType = "opposite",ZoomIn = c(96457229,96457253 ), adjustExprLevel = T , split = T)

##Xrn2##
plotGeneCov(radar,geneName = "ENSMUSG00000027433", center = median, libraryType = "opposite", adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000027433", center = median, libraryType = "opposite",ZoomIn = c(146878646,146878670 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000027433", center = median, libraryType = "opposite",ZoomIn = c(146878646,146878670 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000027433", libraryType = "opposite",ZoomIn = c(146878646,146878670 ), adjustExprLevel = T , split = T)

##Gpx4##
plotGeneCov(radar,geneName = "ENSMUSG00000075706", center = median, libraryType = "opposite", adjustExprLevel = T)
plotGeneCov(radar,geneName = "ENSMUSG00000075706", center = median, libraryType = "opposite",ZoomIn = c(79891800,79892273 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000075706", center = median, libraryType = "opposite",ZoomIn = c(79891800,79892273 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000075706", libraryType = "opposite",ZoomIn = c(79891800,79892273 ), adjustExprLevel = T , split = T)
