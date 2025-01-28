# Load required R packages
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}
library(clusterProfiler)
library(org.Hs.eg.db)

# Load data
data <- read.table("/Users/tianlimin/gene_summary.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Inspect data structure
str(data)
summary(data)

# Process gene names and clean data
gene_list <- unique(unlist(strsplit(as.character(data$Gene.refGene), ";")))
gene_list <- gene_list[gene_list != "NONE"]

# Convert gene names to Entrez IDs
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Perform KEGG enrichment analysis
kegg_results <- enrichKEGG(
  gene = gene_entrez$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Visualize results: Bar plot
barplot(kegg_results, 
        showCategory = 10, 
        title = "Top 10 Enriched KEGG Pathways by Gene Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Visualize results: Dot plot
dotplot(kegg_results, 
        showCategory = 10, 
        title = "Significance and Gene Count in KEGG Pathways") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Install and load pathview for pathway visualization
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}
library(pathview)

# Visualize a specific pathway (e.g., hsa05200: Cancer Pathway)
pathview(gene.data = gene_entrez$ENTREZID, pathway.id = "hsa05200", species = "hsa")

# Sort KEGG results by gene count and select top pathways
kegg_results_sorted <- kegg_results[order(kegg_results$Count, decreasing = TRUE), ]
top_pathways <- head(kegg_results_sorted, 5)$ID

# Batch visualize top pathways
for (id in top_pathways) {
  pathview(gene.data = gene_entrez$ENTREZID, pathway.id = id, species = "hsa")
}

# Display interactive table
library(DT)
datatable(as.data.frame(kegg_results), 
          options = list(pageLength = 10, autoWidth = TRUE), 
          caption = "Interactive KEGG Enrichment Results")
