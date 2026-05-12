################## RNA-Seq Data Analysis with GEO ###########################

##################### Installing and Loading Packages #####################

# Install required packages.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Sc.sgd.db")

#  org.Sc.sgd.db is installed from Bioconductor because it is not available on CRAN.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")

# clusterProfiler is installed from Bioconductor because it is not available on CRAN.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")

# pathview is installed from Bioconductor because it is not available on CRAN.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")

# EnhancedVolcano is installed from Bioconductor because it is not available on CRAN.

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# DESeq2 is installed from Bioconductor because it is not available on CRAN.   

install.packages("tidyverse")
install.packages("pheatmap")

# Load required libraries.

library(DESeq2)
library(org.Sc.sgd.db)
library(pheatmap)
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(EnhancedVolcano)

# DESeq2 is used for differential expression analysis.
# org.Sc.sgd.db provides Saccharomyces cerevisiae gene annotations.
# pheatmap is used for heatmap visualization.
# tidyverse is used for data manipulation.
# clusterProfiler and pathview are used for enrichment analysis.

######################## Loading Data ##################################

# Read gene count data from a text file.

count_data <- read.delim("GSE165308_gene_count.txt", row.names = 1)

# Select the first six columns containing count data.

count_data <- dplyr::select(count_data, 1:6 )

##################### Defining Experimental Conditions ######################

# Create sample metadata for DESeq2 analysis.

sample_info <- data.frame(condition=c( "control", "control", "control", "knockout", "knockout", "knockout"))

rownames(sample_info) <- colnames(count_data)

# Match sample names between metadata and count matrix.

## Checking Data Compatibility ##

# Verify that sample metadata and count matrix are aligned.

all(rownames(sample_info)==colnames(count_data))

##################### Creating DESeqDataSet ###############################

# Create DESeqDataSet object for DESeq2 analysis.

dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~condition)

dds

## Filtering Low Count Genes ## 

# Remove genes with low read counts.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## Setting Reference Level ##

# Set control as the reference condition.

dds$condition <- relevel(dds$condition, ref = "control") 

############# Running Differential Expression Analysis ###############

# Run DESeq2 differential expression analysis.

dds <- DESeq(dds) 

# Extract differential expression results.

res <- results(dds)

res 

########################## PCA Plot Generation ##############################

# Perform variance stabilizing transformation (VST).

vstdata <- vst(dds, blind = F) 

# Generate PCA plot using transformed data.

plotPCA(vstdata, intgroup="condition") 

## Creating MA Plot ## 

# Generate MA plot for differential expression results.

plotMA(dds, alpha=0.05) 

############## Converting Results into Dataframe ###################

# Convert DESeq2 results into a dataframe.

res_df <- as.data.frame(res) 

# Label significantly differentially expressed genes.

res_df$significant <- ifelse( res_df$padj < 0.05, "yes", "no") 

## Designing MA Plot with ggplot2 ##

# Visualize differential expression results using ggplot2 .

ggplot(res_df, aes(log(baseMean), log2FoldChange, color=significant)) +
  geom_point() 

## Retrieving Gene Symbols ##

# Map ENSEMBL IDs to gene symbols using org.Hs.eg.db.

res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = rownames(res_df),
                        keytype = "ENSEMBL",
                        column = "SYMBOL") 


## Enhanced Volcano Plot ##

# Generate volcano plot using EnhancedVolcano.

EnhancedVolcano(res_df, x="log2FoldChange", y="padj", lab = res_df$symbol) 

############### Obtaining Normalized Counts ##########################

# Extract normalized count values.

normalized_counts <- counts(dds, normalized=T) 

# Convert normalized counts into a dataframe.

normalized_counts <- as.data.frame(normalized_counts)

# Merge normalized counts with DESeq2 results.

res_normalize <- merge(res_df, normalized_counts, by=0) 

# Select top 20 significant genes.
# Remove missing values and sort genes by significance.
top_20 <- res_normalize %>% 
  drop_na(symbol) %>% 
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>% 
  arrange(padj) %>% 
  select(9:15) %>% 
  head(20)

# Assign gene symbols as row names.

rownames(top_20) <- top_20$symbol
top_20$symbol <- NULL 

# Generate heatmap using log2-transformed normalized counts.

pheatmap(log2(top_20 + 1)) 

#################### Gene Ontology and KEGG Pathway #######################

# Functional enrichment analysis is used to interpret gene expression results.

## Gene Ontology Enrichment ## 

# Filter significantly expressed genes for enrichment analysis.

res_go <- res_df %>%
  filter(significant == "yes")

# Create gene list for enrichment analysis.

gene_list <- rownames(res_go) 

## Biological Process Ontology Analysis ##

# Perform Gene Ontology enrichment analysis for Biological Process terms.

ego1 <- enrichGO(gene=gene_list,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENSEMBL",
                 ont = "BP",
                 readable = T) 

# Visualize Biological Process enrichment results. 

barplot(ego1) 

## Molecular Function Ontology Analysis ## 

# Perform Gene Ontology enrichment analysis for Molecular Function terms.

ego2 <- enrichGO(gene=gene_list,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENSEMBL",
                 ont = "MF",
                 readable = T) 

# Visualize Molecular Function enrichment results.

barplot(ego2) 

## Cellular Component Ontology Analysis ##

# Perform Gene Ontology enrichment analysis for Cellular Component terms.

ego3 <- enrichGO(gene=gene_list,
                 OrgDb = "org.Hs.eg.db",
                 keyType = "ENSEMBL",
                 ont = "CC",
                 readable = T) 

# Visualize Cellular Component enrichment results.

ego3_plot <- barplot(ego3) 

ego3_plot 

# Save Cellular Component enrichment plot.

png("ego3_plot.png", width = 1200, height = 600, res = 150)
print(ego3_plot)
dev.off()

# Visualize GO term-gene relationships using goplot().

options(ggrepel.max.overlaps = Inf)
goplot(ego1) 

## KEGG Pathway Enrichment Analysis ##

# Search for Homo sapiens organism information in KEGG.

search <- search_kegg_organism("Homo sapiens") 

search 

# Map ENSEMBL IDs to UniProt IDs.

res_go$uniprot <- mapIds(org.Hs.eg.db,
                         keys = rownames(res_go),
                         keytype = "ENSEMBL",
                         column = "UNIPROT")

# Store UniProt IDs for KEGG analysis. 

uniprot_id <- res_go$uniprot 

# Perform KEGG pathway enrichment analysis.

kegg_enrich <- enrichKEGG(gene = uniprot_id,
                          organism = "hsa",
                          keyType = "uniprot") 

head(kegg_enrich) 

# Display KEGG enrichment results.

browseKEGG(kegg_enrich, "hsa03008")

# Open KEGG pathway visualization in browser.

# Visualize KEGG pathway using pathview().

hsa03008 <- pathview(gene.data = uniprot_id,
                     species = "hsa",
                     pathway.id = "hsa03008",
                     gene.idtype = "uniprot") 





