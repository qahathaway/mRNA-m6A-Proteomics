library(pathfindR)

input <- read.csv(file = 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_m6A_pathfindR.csv', header = TRUE)
input2 <- read.csv(file = 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_pathfindR.csv', header = TRUE)

pathObject <- run_pathfindR(input, gene_sets = "mmu_KEGG")
pathObject2 <- run_pathfindR(input2, gene_sets = "mmu_KEGG")

write.csv(pathObject, 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_m6A_pathfindR_results.csv')
write.csv(pathObject2, 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_pathfindR_results.csv')

ModpathObject <- read.csv(file = 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_m6A_pathfindR_results.csv', header = TRUE)
ModpathObject2 <- read.csv(file = 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_pathfindR_results.csv', header = TRUE)

clustered1 <-cluster_enriched_terms(ModpathObject, use_description = TRUE)
clustered2 <-cluster_enriched_terms(ModpathObject2, use_description = TRUE)

write.csv(clustered1, 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_m6A_pathfindR_results_clustered.csv')
write.csv(clustered2, 'C:\\Users\\ak0086\\Documents\\Lab\\6mA and m6A sequencing project\\Amina_Final\\mRNA_pathfindR_results_clustered.csv')

enrichment_chart(clustered1[1:10, ], plot_by_cluster = TRUE)
enrichment_chart(clustered2[1:10, ], plot_by_cluster = TRUE)

enrichment_chart(pathObject, top_terms = 10)
enrichment_chart(pathObject2, top_terms = 10)

UpSet_plot(clustered1, input, num_terms = 10, method = "heatmap", use_description = TRUE, low = "green", mid = "black", high = "red")
UpSet_plot(clustered2, input2, num_terms = 10, method = "heatmap", use_description = TRUE, low = "green", mid = "black", high = "red")