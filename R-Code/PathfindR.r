library(pathfindR)

input <- read.csv(file = 'path/to/file.csv', header = TRUE)
input2 <- read.csv(file = 'path/to/file.csv', header = TRUE)

ModpathObject <- run_pathfindR(input, gene_sets = "mmu_KEGG")
ModpathObject2 <- run_pathfindR(input2, gene_sets = "mmu_KEGG")

write.csv(ModpathObject, 'path/to/file/results.csv')
write.csv(ModpathObject2, 'path/to/file/results.csv')

clustered1 <-cluster_enriched_terms(ModpathObject, use_description = TRUE)
clustered2 <-cluster_enriched_terms(ModpathObject2, use_description = TRUE)

write.csv(clustered1, 'path/to/file/results_clustered.csv')
write.csv(clustered2, 'path/to/file/results_clustered.csv')

enrichment_chart(clustered1[1:10, ], plot_by_cluster = TRUE)
enrichment_chart(clustered2[1:10, ], plot_by_cluster = TRUE)

enrichment_chart(pathObject, top_terms = 10)
enrichment_chart(pathObject2, top_terms = 10)

UpSet_plot(clustered1, input, num_terms = 10, method = "heatmap", use_description = TRUE, low = "green", mid = "black", high = "red")
UpSet_plot(clustered2, input2, num_terms = 10, method = "heatmap", use_description = TRUE, low = "green", mid = "black", high = "red")
