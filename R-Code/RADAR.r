
###FULL CODE###

library(RADAR)
samplenames <- c( paste0("FSham_",1:3), paste0("FExp_",1:3) )

bam_dir <- system.file('extdata',package = 'AminaM6A')
samplenames

radar <- countReads(samplenames = samplenames,
                    gtf = "/home/john/R/x86_64-pc-linux-gnu-library/4.1/AminaM6A/extdata/Mus_musculus.GRCm39.104.gtf",
                    bamFolder = bam_dir,
                    modification = "m6A",
                    fragmentLength = 297,
                    outputDir = "/home/john/R/x86_64-pc-linux-gnu-library/4.1/AminaM6A/reads",
                    saveOutput = T,
                    threads = 32,
                    binSize = 25
)

summary(radar)


radar <- normalizeLibrary(radar)
sizeFactors(radar)
radar <- adjustExprLevel(radar)
#radar <- adjustExprLevel(radar, adjustBy = "pos")
variable(radar) <- data.frame( Group = c(rep("FS_",3),rep("FE_",3)) )
radar <- filterBins(radar,minCountsCutOff = 5)


radar <- diffIP_parallel(radar, thread = 32)
#radar <- diffIP_parallel(radar, thread = 32, exclude = "FS_3")
#radar <- diffIP_parallel(radar, thread = 32, exclude = "FE_2")


top_bins <- extractIP(radar,filtered = T)[order(rowMeans( extractIP(radar,filtered = T) ),decreasing = T)[1:10000],]
plotPCAfromMatrix(top_bins,group = unlist(variable(radar)) )

#radar <- reportResult(radar, cutoff = 0.05, Beta_cutoff = 0.5, thread = 32)
radar <- reportResult(radar, cutoff = 0.1, Beta_cutoff = 0.5, thread = 32)
#radar <- reportResult(radar, cutoff = 0.05, Beta_cutoff = 0, thread = 32)
#radar <- reportResult(radar, cutoff = 0.1, Beta_cutoff = 0, thread = 32)
#radar <- reportResult(radar, cutoff = 0.2, Beta_cutoff = 0, thread = 32)

result <- results(radar)
head(result)

write.csv(result, file = "/home/john/Documents/AminaM6A/AminaM6A/Results/Sham3_vs_Exp3_DE_ALL.csv")

radar <- PrepCoveragePlot(radar, gtf = "/home/john/R/x86_64-pc-linux-gnu-library/4.1/AminaM6A/extdata/Mus_musculus.GRCm39.104.gtf")
plotHeatMap(radar)
peakDistribution(radar)

RNASeqCounts <- geneExpression(radar)

IPreadsn <- extractIP(radar, normalized = TRUE)
IPreadsa <- extractIP(radar, adjusted = TRUE)
IPreadsf <- extractIP(radar, filtered = TRUE)
write.csv(IPreadsn, file = "/home/john/Documents/Amina_Final/m6A/IP_Reads_Normalized.csv")
write.csv(IPreadsa, file = "/home/john/Documents/Amina_Final/m6A/IP_Reads_Adjusted.csv")
write.csv(IPreadsf, file = "/home/john/Documents/Amina_Final/m6A/IP_Reads_Filtered.csv")

Inputreads <- extractInput(radar, normalized = TRUE)
write.csv(Inputreads, file = "/home/john/Documents/Amina_Final/m6A/Input_Reads_Normalized.csv")

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

##Erp29##
plotGeneCov(radar,geneName = "ENSMUSG00000029616", center = median, libraryType = "opposite", adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000029616", center = median, libraryType = "opposite",ZoomIn = c(121568813,121568820 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000029616", center = median, libraryType = "opposite",ZoomIn = c(121568813,121568820 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000029616", libraryType = "opposite",ZoomIn = c(121568813,121568820 ), adjustExprLevel = T , split = T)

##Lysmd4##
plotGeneCov(radar,geneName = "ENSMUSG00000043831", center = median, libraryType = "opposite", adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000043831", center = median, libraryType = "opposite",ZoomIn = c(66884852,66884876 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000043831", center = median, libraryType = "opposite",ZoomIn = c(66884852,66884876 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000043831", libraryType = "opposite",ZoomIn = c(66884852,66884876 ), adjustExprLevel = T , split = T)

##Flad1##
plotGeneCov(radar,geneName = "ENSMUSG00000042642", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000042642", center = median, libraryType = "opposite",ZoomIn = c(89310148,89310172 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000042642", center = mean, libraryType = "opposite",ZoomIn = c(89310148,89310172 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000042642", libraryType = "opposite",ZoomIn = c(89310148,89310172 ), adjustExprLevel = T , split = T)

##Cox20##
plotGeneCov(radar,geneName = "ENSMUSG00000026500", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000026500", center = median, libraryType = "opposite",ZoomIn = c(178149681,178149705 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000026500", center = mean, libraryType = "opposite",ZoomIn = c(178149681,178149705 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000026500", libraryType = "opposite",ZoomIn = c(178149681,178149705 ), adjustExprLevel = T , split = T)

##Gpx4##
plotGeneCov(radar,geneName = "ENSMUSG00000075706", center = median, libraryType = "opposite", adjustExprLevel = T)
plotGeneCov(radar,geneName = "ENSMUSG00000075706", center = median, libraryType = "opposite",ZoomIn = c(79891800,79892273 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000075706", center = median, libraryType = "opposite",ZoomIn = c(79891800,79892273 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000075706", libraryType = "opposite",ZoomIn = c(79891800,79892273 ), adjustExprLevel = T , split = T)

##Sod2##
plotGeneCov(radar,geneName = "ENSMUSG00000006818", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000006818", center = median, libraryType = "opposite",ZoomIn = c(13234000,13237691 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000006818", center = mean, libraryType = "opposite",ZoomIn = c(13234000,13237691 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000006818", libraryType = "opposite",ZoomIn = c(13234000,13237691 ), adjustExprLevel = T , split = T)

##Cat##
plotGeneCov(radar,geneName = "ENSMUSG00000027187", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000027187", center = median, libraryType = "opposite",ZoomIn = c(103284194,103287194 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000027187", center = mean, libraryType = "opposite",ZoomIn = c(103284194,103287194 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000027187", libraryType = "opposite",ZoomIn = c(103284194,103287194 ), adjustExprLevel = T , split = T)

##Dnmt1##
plotGeneCov(radar,geneName = "ENSMUSG00000004099", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000004099", center = median, libraryType = "opposite",ZoomIn = c(20823576,20823888 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000004099", center = mean, libraryType = "opposite",ZoomIn = c(20823576,20823888 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000004099", libraryType = "opposite",ZoomIn = c(20823576,20823888 ), adjustExprLevel = T , split = T)

##Fam102b##
plotGeneCov(radar,geneName = "ENSMUSG00000040339", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000040339", center = median, libraryType = "opposite",ZoomIn = c(108878565,108878714 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000040339", center = mean, libraryType = "opposite",ZoomIn = c(108878565,108878714 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000040339", libraryType = "opposite",ZoomIn = c(108878565,108878714 ), adjustExprLevel = T , split = T)

##Igf1r##
plotGeneCov(radar,geneName = "ENSMUSG00000005533", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000005533", center = median, libraryType = "opposite",ZoomIn = c(67882176,67882275 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000005533", center = mean, libraryType = "opposite",ZoomIn = c(67882176,67882275 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000005533", libraryType = "opposite",ZoomIn = c(67882176,67882275 ), adjustExprLevel = T , split = T)

##Ndufs4##
plotGeneCov(radar,geneName = "ENSMUSG00000021764", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000021764", center = median, libraryType = "opposite",ZoomIn = c(114424689,114424738 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000021764", center = mean, libraryType = "opposite",ZoomIn = c(114424689,114424738 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000021764", libraryType = "opposite",ZoomIn = c(114424689,114424738 ), adjustExprLevel = T , split = T)

##Ndufab1##
plotGeneCov(radar,geneName = "ENSMUSG00000030869", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000030869", center = median, libraryType = "opposite",ZoomIn = c(121686359,121686408 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000030869", center = mean, libraryType = "opposite",ZoomIn = c(121686359,121686408 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000030869", libraryType = "opposite",ZoomIn = c(121686359,121686408 ), adjustExprLevel = T , split = T)

##Ndufa11##
plotGeneCov(radar,geneName = "ENSMUSG00000002379", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000002379", center = median, libraryType = "opposite",ZoomIn = c(57028600,57028649 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000002379", center = mean, libraryType = "opposite",ZoomIn = c(57028600,57028649 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000002379", libraryType = "opposite",ZoomIn = c(57028600,57028649 ), adjustExprLevel = T , split = T)

##Cox7c##
plotGeneCov(radar,geneName = "ENSMUSG00000017778", center = median, libraryType = "opposite")
plotGeneCov(radar,geneName = "ENSMUSG00000017778", center = median, libraryType = "opposite",ZoomIn = c(86194434,86194483 ) )
plotGeneCov(radar,geneName = "ENSMUSG00000017778", center = mean, libraryType = "opposite",ZoomIn = c(86194434,86194483 ), adjustExprLevel = T )
plotGeneCov(radar,geneName = "ENSMUSG00000017778", libraryType = "opposite",ZoomIn = c(86194434,86194483 ), adjustExprLevel = T , split = T)
