Bulk RNA sequencing analysis of endothelial cells
================
SChear
2024-01-01

``` r
library(DESeq2)
library(ggplot2)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(clusterProfiler)
library(EnhancedVolcano)
library(fgsea)
library(magrittr)
library(tidyverse)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(edgeR)
library(matrixStats)
library(circlize)
```

``` r
setwd("/Users/sueannechear/Bioinformatics/rnaseq/iMac/quant/second_analysis/Endothelialcells github")
counts<-read.delim("cleancount_V3.csv",header=T, sep=",")
rownames(counts)  <- counts[,1]
counts<- counts[ -c(1) ]
head(counts,3)
```

Samples used in this analysis:

``` r
df<-data.frame("Abbreviation"=c("iPSC","HUVEC","iBMEC","HBMEC","BEPC","CEC"),
               "Cell"=c("Human induced pluripotent stem cells","Human umbilical vein endothelial cells","iPSC derived brain microvascular endothelial cells","Human brain microvascular endothelial cells","Bronchial epithelial cells","Colon epithelial cells"),
               "Source"=c("PRJNA754196","GSE57662","GSE129290","GSE97575","GSE85402","GSE82207"), "Replicates"=c(3,3,3,3,2,2))
print(df)
```

    ##   Abbreviation                                               Cell      Source
    ## 1         iPSC               Human induced pluripotent stem cells PRJNA754196
    ## 2        HUVEC             Human umbilical vein endothelial cells    GSE57662
    ## 3        iBMEC iPSC derived brain microvascular endothelial cells   GSE129290
    ## 4        HBMEC        Human brain microvascular endothelial cells    GSE97575
    ## 5         BEPC                         Bronchial epithelial cells    GSE85402
    ## 6          CEC                             Colon epithelial cells    GSE82207
    ##   Replicates
    ## 1          3
    ## 2          3
    ## 3          3
    ## 4          3
    ## 5          2
    ## 6          2

Create a DESeq dataset : 18 samples; 62753 genes

``` r
condition<-factor(c("HUVEC","HUVEC","HUVEC","iEC","iEC",
                    "BEPC","BEPC",
                    "iPSC","iPSC","iPSC","HBMEC","HBMEC","HBMEC",
                    "CEC","CEC","iBMEC","iBMEC","iBMEC"
))

sample<-factor(colnames(counts))

coldata<-data.frame(sample,condition)

dds<-DESeqDataSetFromMatrix(countData = counts,
                            colData = coldata,
                            design=~condition)

dds
```

    ## class: DESeqDataSet 
    ## dim: 62753 18 
    ## metadata(1): version
    ## assays(1): counts
    ## rownames(62753): ENSG00000000003 ENSG00000000005 ... ENSG00000292372
    ##   ENSG00000292373
    ## rowData names(0):
    ## colnames(18): HUVEC_1 HUVEC_2 ... iBMEC.1 iBMEC.2
    ## colData names(2): sample condition

We retain genes with at least 1 count per million (CPM) in at least two
samples. Genes remained after filtering: 20362

``` r
dds = dds[ rowSums(edgeR::cpm(counts(dds)) > 1)>=2, ]

nrow(dds)
```

    ## [1] 20362

QC for dispersion of variability in data. Fitted line trend below 1,
indicating the data is a good fit for the DESeq model.

``` r
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)
plotDispEsts(dds)
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Transform data using VST method before PCA.

``` r
vstdata<-vst(dds,blind=F)

meanSdPlot(assay(vstdata), ranks=FALSE)
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

PCA: to examine variation between samples. Gnerally,HUVEC,iEC and HBMEC
are clustered close together,indicating similar transcriptomics. PC1
clearly separates primary endothelial cells from other cell types. iPSC
and epithelial cells are more similar in transcriptomics than they are
with endothelial cells.

``` r
plotPCA(vstdata,intgroup="condition")
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

``` r
#check pc3 and pc4
plotPCA <- function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3], PC4 = pca$x[, 4], group = group, 
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    ggplot(data = d, aes_string(x = "PC3", y = "PC4", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC3: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC4: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed()
}

#print(plotPCA(vstdata,intgroup="condition"))
```

``` r
ntop <- 500
rv <- rowVars(assay(vstdata))
select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
mat <- t( assay(vstdata)[select, ] )

pca<-prcomp(mat)
pca <- as.data.frame(pca$x)
```

Which genes impact the most in PC1?

``` r
getLoadings = function(dds){
  
  mat<-assay(vstdata)
  pca = prcomp(t(mat), retx = TRUE)
  
  return(pca$rotation)
}

loadings_vstdata = getLoadings(vstdata) %>% as.data.frame()
# Annotating gene names
loadings_vstdata$symbol = mapIds(org.Hs.eg.db,
                              keys=rownames(loadings_vstdata),
                              column="SYMBOL",
                              keytype="ENSEMBL",
                              multiVals="first")
```

    ## 'select()' returned 1:many mapping between keys and columns

``` r
# show the top 10 genes from PC1
loadings_vstdata %>% 
  # select only the PCs we are interested in
  dplyr::select(symbol, PC1) %>%
  # convert to "long" format
  pivot_longer(cols = "PC1", names_to = "PC1", values_to = "loadings") %>% 
  # for PC1
  group_by(PC1) %>% 
  # arrange by descending order
  arrange(desc(abs(loadings))) %>% 
  # take the 10 top rows
  slice(1:20) %>%
  pull(symbol)
```

    ##  [1] "KRT8"    "KRT19"   "KRT18"   "SLC1A2"  "CD24"    "FAM107A" "SPARCL1"
    ##  [8] "GFAP"    "A2M"     "ADGRF5"  "CD34"    "ATP1B2"  "CDH1"    "SLC6A1" 
    ## [15] "VWF"     "CRABP2"  "PMP2"    "PECAM1"  "TOP2A"   "CLDN6"

Which genes impact the most in PC2?

``` r
# show the top 20 genes from PC2
loadings_vstdata %>% 
  # select only the PCs we are interested in
  dplyr::select(symbol, PC2) %>%
  # convert to "long" format
  pivot_longer(cols = "PC2", names_to = "PC2", values_to = "loadings") %>% 
  # for PC2
  group_by(PC2) %>% 
  # arrange by descending order
  arrange(desc(abs(loadings))) %>% 
  # take the 10 top rows
  slice(1:20) %>%
  pull(symbol)
```

    ##  [1] "MMP1"    "CD93"    "ANKRD1"  "MMRN1"   "ESM1"    "PECAM1"  "PTPRZ1" 
    ##  [8] "KIF1A"   "THBS1"   "ATP1A2"  "PTX3"    "CLEC14A" "EFEMP1"  "CDH5"   
    ## [15] "PLP1"    "PLVAP"   "SOX18"   "CLDN6"   "CDH1"    "MMRN2"

Cluster dendrogram

``` r
rv <- rowVars(assay(vstdata))
o <- order(rv,decreasing=TRUE)
dists <- dist(t(assay(vstdata)[head(o,500),]))
hc <- hclust(dists)
plot(hc, labels=vstdata$sample)
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

Correlation matrix heat map of transcript expression across all samples.
Euclidean distance is used as the similarity measure and clustering
samples based on the ‘complete’ method. The lower the numbers, the
stronger the correlation between samples.

``` r
sampleDists<-dist(t(assay(vstdata)))
sampleDistMatrix<-as.matrix(sampleDists)
colnames(sampleDistMatrix)

colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, display_numbers = TRUE,
         clustering_distance_cols=sampleDists, col=colors, fontsize_number=7)
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-14-1.png)<!-- -->

## Heatmap to visualize clustering using top 100 genes.

``` r
#get the indices of the top variable genes
topVarGenes <- head(order(rowVars(assay(vstdata)), decreasing = TRUE), 100)

#subset the data matrix to include only the top variable genes
mat  <- assay(vstdata)[ topVarGenes, ]


#center the data
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vstdata))

#map ensembl IDs to gene symbols
symbols <- mapIds(org.Hs.eg.db, keys = rownames(mat), column = "SYMBOL", keytype = "ENSEMBL")

symbols <- symbols[!is.na(symbols)]
  symbols <- symbols[match(rownames(mat), names(symbols))]
  rownames(mat) <- symbols
  keep <- !is.na(rownames(mat))
  mat <- mat[keep,]


#create a heatmap with hierarchical clustering
#heatmap_result <- pheatmap(mat, annotation_col = anno, fontsize_row=5)
```

``` r
# Perform hierarchical clustering separately
hc_rows <- hclust(dist(mat), method = "complete")

# Create a heatmap with hierarchical clustering
#heatmap_result <- pheatmap(mat, annotation_col = anno, fontsize_row = 5, clustering_distance_rows = "correlation")

heatmap_result <- pheatmap(mat, fontsize_row = 5, clustering_distance_rows = "correlation")
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-16-1.png)<!-- -->

``` r
# Extract the cluster assignments for the rows
cluster_assignments <- cutree(hc_rows, k = 2)

# Print or use the cluster assignments as needed
print(cluster_assignments)

# Create a data frame with gene symbols and cluster assignments
gene_cluster_df <- data.frame(
  GeneSymbol = names(cluster_assignments),
  ClusterAssignment = cluster_assignments
)

# Order the data frame by cluster assignments
gene_cluster_df <- gene_cluster_df[order(gene_cluster_df$ClusterAssignment, gene_cluster_df$GeneSymbol), ]

# Print or use the data frame as needed
print(gene_cluster_df)

write.csv(gene_cluster_df, "genecluster.csv",row.names=F)
```

Pathways enriched in cluster 1 :these are genes upregulated in
endothelial cells.

``` r
# Select genes in Cluster 1
genes_in_cluster1 <- gene_cluster_df$GeneSymbol[gene_cluster_df$ClusterAssignment == 1]

print(genes_in_cluster1)

entrez_ids <- mapIds(org.Hs.eg.db, keys = genes_in_cluster1, keytype = "SYMBOL", column = "ENTREZID")

GO_results <- enrichGO(gene = genes_in_cluster1, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results)


plot(barplot(GO_results, showCategory = 15,cex.names = 0.5))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Pathways enriched in cluster 2; these are genes upregulated in
non-endothelial cells.

``` r
# Select genes in Cluster 2
genes_in_cluster2 <- gene_cluster_df$GeneSymbol[gene_cluster_df$ClusterAssignment == 2]

print(genes_in_cluster2)

entrez_ids <- mapIds(org.Hs.eg.db, keys = genes_in_cluster2, keytype = "SYMBOL", column = "ENTREZID")

GO_results2 <- enrichGO(gene = genes_in_cluster2, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL", ont = "BP")

as.data.frame(GO_results2)


plot(barplot(GO_results2, showCategory = 15,cex.names = 0.5))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-18-1.png)<!-- -->

### Differential expression of gene analysis:

### iEC vs iBMEC. (padjusted value=0.05)

``` r
pds<-dds

pds$condition<-relevel(pds$condition, ref="iBMEC")

pds<-DESeq(pds)

res = results(pds, contrast=c("condition","iEC","iBMEC"), alpha=0.05)

summary(res)
```

    ## 
    ## out of 20362 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 3358, 16%
    ## LFC < 0 (down)     : 4922, 24%
    ## outliers [1]       : 1444, 7.1%
    ## low counts [2]     : 0, 0%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
sigs<-na.omit(res)

sigs.df<-as.data.frame(sigs)

columns(org.Hs.eg.db)

sigs.df$symbol<-mapIds(org.Hs.eg.db, keys=rownames(sigs.df), keytype = "ENSEMBL",column="SYMBOL")

sigs.df<-sigs.df%>%filter(!str_detect(symbol,'NA'))
```

Retrieve gene ontology terms associated with upregulated genes in iEC
samples

``` r
genes_to_test <- rownames(sigs[sigs$log2FoldChange >5,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)


plot(barplot(GO_results, showCategory = 15))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

Retrieve gene ontology terms associated with downregulated genes in iEC
samples

``` r
genes_to_test <- rownames(sigs[sigs$log2FoldChange < -2,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)


plot(barplot(GO_results, showCategory = 15))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-22-1.png)<!-- -->

### iBMEC vs iPSC. (padjusted value=0.05)

``` r
pds<-dds

pds$condition<-relevel(pds$condition, ref="iPSC")

pds<-DESeq(pds)

res = results(pds, contrast=c("condition","iBMEC","iPSC"), alpha=0.05)

summary(res)
```

    ## 
    ## out of 20362 with nonzero total read count
    ## adjusted p-value < 0.05
    ## LFC > 0 (up)       : 4182, 21%
    ## LFC < 0 (down)     : 4232, 21%
    ## outliers [1]       : 1448, 7.1%
    ## low counts [2]     : 0, 0%
    ## (mean count < 3)
    ## [1] see 'cooksCutoff' argument of ?results
    ## [2] see 'independentFiltering' argument of ?results

``` r
sigs<-na.omit(res)

sigs.df<-as.data.frame(sigs)

columns(org.Hs.eg.db)

sigs.df$symbol<-mapIds(org.Hs.eg.db, keys=rownames(sigs.df), keytype = "ENSEMBL",column="SYMBOL")

sigs.df<-sigs.df%>%filter(!str_detect(symbol,'NA'))
```

Retrieve gene ontology terms associated with upregulated genes in iBMEC
samples

``` r
genes_to_test <- rownames(sigs[sigs$log2FoldChange >5,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)


plot(barplot(GO_results, showCategory = 15))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

Retrieve gene ontology terms associated with downregulated genes in
iBMEC samples

``` r
genes_to_test <- rownames(sigs[sigs$log2FoldChange < -2,])

GO_results <- enrichGO(gene = genes_to_test, OrgDb = "org.Hs.eg.db", keyType = "ENSEMBL", ont = "BP")

as.data.frame(GO_results)


plot(barplot(GO_results, showCategory = 15))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-26-1.png)<!-- -->

### Heatmap expression of various cell marker genes

``` r
# Extract the transformed data matrix
transformed_matrix <- assay(vstdata)

# Convert the matrix to a data frame
transformed_df <- as.data.frame(transformed_matrix)

# Extract row names (Ensembl IDs)
ensembl_ids <- rownames(transformed_df)

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
symbols <- mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL")

# Add gene symbols as a new column to the data frame
transformed_df$Symbol <- symbols[match(ensembl_ids, names(symbols))]

# Remove the initial row names
rownames(transformed_df) <- NULL

#save as csv
write.csv(transformed_df, "forheatmap.csv",row.names=F,sep='\t')
```

``` r
df<-read.delim("forheatmap.csv",header=T, sep=",")
head(df)


df2<-df%>%filter(Symbol%in% c("TMEM119","P2RY12","HEXB","FCRLS","SALL1","C1Q","GPR34","OLFML3","MERTK",
  "PROS1","TYRO3","TGFBR1","CD31","NG2","PDGFRB","CD146","NESTIN","VWF",
  "ACE","ADAMTS13","PECAM1","VCAM1","ICAM1","ICAM2","CD47","SELE","SLP",
  "CDH5","NECTIN2","ESAM","LEF1","FZD3","NOTUM","APCDD1","AXIN2","DIXDC1",
  "TNFRSF19","MECOM")) 

head(df2)

df3 <- df2 %>% 
  column_to_rownames(var = "Symbol")

#set a color scheme
#colors<-colorRampPalette(rev(brewer.pal(9,"Blues")))(255)

pheatmap(df3, scale="row",display_numbers=TRUE,fontsize_number=7,color = colorRampPalette(rev(c("#D73027", "#FC8D59", 
        "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-28-1.png)<!-- -->

``` r
pheatmap(df3, cluster_cols=FALSE,  display_numbers=TRUE,fontsize_number=7,
          color = colorRampPalette(rev(c("#D73027", "#FC8D59", 
        "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-28-2.png)<!-- -->

``` r
pheatmap(df3, cluster_cols=FALSE, cluster_rows = FALSE, display_numbers=TRUE,fontsize_number=7,
          color = colorRampPalette(rev(c("#D73027", "#FC8D59", 
        "#FEE090", "#FFFFBF", "#E0F3F8", "#91BFDB", "#4575B4")))(100))
```

![](Endothelial_cells_files/figure-gfm/unnamed-chunk-28-3.png)<!-- -->
