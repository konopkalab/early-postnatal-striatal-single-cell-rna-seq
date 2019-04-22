# DATA CLUSTERING

This script performs data clustering on single cell RNA-seq UMI counts dataset in multiple steps such as filtering, normalization, scaling using Seurat pipeline (Seurat v2.3.4). Towards the end, gene markers per cluster are also determined using seurat functions. *sessionInfo* is provided towards the end of the script.



#### LOAD QC DATA

Load seurat object created during data quality check for downstream processing.

```{r}
load("STR_SEURAT_DATA_QC.RData")
```



#### FILTERING THE DATA

Filter data based on cutoffs for number of UMIs, number of genes and percent mitochondrial genes using the QC plots  and save updated seurat object.

```{r}
## set the filtering cutoffs
nUlo <- -Inf # cutoff for number of UMIs - lower bound
nUhi <- 50000 # cutoff for number of UMIs - higher bound
nGlo <- -Inf # cutoff for number of genes - lower bound
nGhi <- Inf # cutoff for number of genes - higher bound
pMlo <- -Inf # cutoff for percent mitochondrial content - lower bound
pMhi <- 0.10 # cutoff for percent mitochondrial content - higher bound

## plot the data showing filtering cutoffs
plotseuObjAll_UM <- ggplot(seuObjAll@meta.data, aes(x = nUMI, y = pMito, colour = Genotype)) + geom_point(shape = 20, size = 2, alpha = 1) + labs(title = paste("ALL_DATA", round(cor(seuObjAll@meta.data$nUMI, seuObjAll@meta.data$pMito), 2), sep = ", "), x = "Number of UMIs", y = "Percent Mito") + annotate("rect", xmin = nUlo, xmax = nUhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "blue") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = pMlo, ymax = pMhi, alpha = 0.1, fill = "green") + theme_bw()
plotseuObjAll_UG <- ggplot(seuObjAll@meta.data, aes(x = nUMI, y = nGene, colour = Genotype)) + geom_point(shape = 20, size = 2, alpha = 1) + labs(title = paste("ALL_DATA", round(cor(seuObjAll@meta.data$nUMI, seuObjAll@meta.data$nGene), 2), sep = ", "), x = "Number of UMIs", y = "Number of Genes") + annotate("rect", xmin = nUlo, xmax = nUhi, ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "red") + annotate("rect", xmin = -Inf, xmax = Inf, ymin = nGlo, ymax = nGhi, alpha = 0.1, fill = "green") + theme_bw()
plotQCseuObjAll2 <- grid.arrange(plotseuObjAll_UM, plotseuObjAll_UG, ncol = 2)
ggsave(paste(seuObjAll@project.name, "_FILT.pdf", sep = ""), plot = plotQCseuObjAll2, width = 10, height = 4, units = "in", dpi = 150)

## data filtering using cutoffs for number of UMIs and percent mitochondrial content
seuObjAllFilt <- FilterCells(object = seuObjAll, subset.names = c("nUMI", "pMito"), low.thresholds = c(nUlo, pMlo), high.thresholds = c(nUhi, pMhi))

## save filtered data along with cutoffs used
save(seuObjAllFilt, nUlo, nUhi, nGlo, nGhi, pMlo, pMhi, file = "STR_SEURAT_DATA_FILT_TEMP.RData")
```



#### REMOVING GENES FROM CHROMOSOMES X, Y AND M

Remove the genes from Chr X, Chr Y to avoid gender biases and Chr M followed by saving the updated seurat object. List of genes from Chr X, Chr Y and Chr M is created using reference annotation (Gencode vM17) used during alignment step.

```{r}
## read the list of genes from ChrX, ChrY and Chr M
mGenes <- scan("GENCODE_vM17_MM10_ChrM_GENES.txt", what = "", sep = "\n") ## 13 genes
xGenes <- scan("GENCODE_vM17_MM10_ChrX_GENES.txt", what = "", sep = "\n") ## 1086 genes
yGenes <- scan("GENCODE_vM17_MM10_ChrY_GENES.txt", what = "", sep = "\n") ## 183 genes
mxyGenes <- unique(sort(c(mGenes, xGenes, yGenes)))

## create the list of genes and cells to retain after filtering
keepGenes <- unique(sort(setdiff(row.names(seuObjAllFilt@data), mxyGenes)))
keepCells <- colnames(seuObjAllFilt@data)

## filter UMI data using genes and cells to retain
strDataFilt <- seuObjAllFilt@raw.data[keepGenes, keepCells]

## initialize the seurat object with the filtered raw (non-normalized data)
seuObjFilt <- CreateSeuratObject(raw.data = strDataFilt, project = "STR_ALL_FILT")

## update meta data for pMito
pMitoFilt <- seuObjAll@meta.data[keepCells, "pMito"]
names(pMitoFilt) <- row.names(seuObjAll@meta.data[keepCells,])
seuObjFilt <- AddMetaData(object = seuObjFilt, metadata = pMitoFilt, col.name = "pMito")

## update meta data for library
libFilt <- seuObjAll@meta.data[keepCells, "Library"]
names(libFilt) <- row.names(seuObjAll@meta.data[keepCells,])
seuObjFilt <- AddMetaData(object = seuObjFilt, metadata = libFilt, col.name = "Library")

## update meta data for genotype
genoFilt <- seuObjAll@meta.data[keepCells, "Genotype"]
names(genoFilt) <- row.names(seuObjAll@meta.data[keepCells,])
seuObjFilt <- AddMetaData(object = seuObjFilt, metadata = genoFilt, col.name = "Genotype")

## save filtered data
save(seuObjFilt, keepCells, keepGenes, file = "STR_SEURAT_DATA_FILT.RData")
```



#### NORMALIZE AND SCALE THE DATA

Normalize and scale the data using Seurat pipeline and save updated seurat object.

```{r}
## log-normalize the data using scaling factor
seuObjFiltNorm <- NormalizeData(object = seuObjFilt, normalization.method = "LogNormalize", scale.factor = 10000)

## identify ~2000 top variable genes
pdf("STR_ALL_VAR_GENES.pdf", width = 8, height = 6)
seuObjFiltNorm <- FindVariableGenes(object = seuObjFiltNorm, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.18, x.high.cutoff = 5, y.cutoff = 0.5)
dev.off()

## scale the data by regressing to remove unintended variation from number of UMIs and percent mitochondrial content per cell
seuObjFiltNorm <- ScaleData(object = seuObjFiltNorm, vars.to.regress = c("nUMI", "pMito"))

## save normalized, scaled and regressed data
save(seuObjFiltNorm, file = "STR_SEURAT_DATA_FILT_NORM.RData")
```



#### PCA ANALYSIS

Run PCA analysis to identify number of top significant PCs to be used for downstream analysis and save updated seurat object.

```{r
## run PCA analysis using identified variable genes
seuObjFiltNormPC <- RunPCA(object = seuObjFiltNorm, pc.genes = seuObjFiltNorm@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)

## examine and visualize PCA results
pdf("STR_ALL_PCA_1.pdf", width = 9, height = 21)
VizPCA(object = seuObjFiltNormPC, pcs.use = 1:50)
dev.off()

pdf("STR_ALL_PCA_2.pdf", width = 9, height = 6)
PCAPlot(object = seuObjFiltNormPC, dim.1 = 1, dim.2 = 2)
dev.off()

seuObjFiltNormPC <- ProjectPCA(object = seuObjFiltNormPC, do.print = FALSE)

pdf("STR_ALL_PCA_3.pdf", width = 9, height = 21)
PCHeatmap(object = seuObjFiltNormPC, pc.use = 1:50, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

pdf("STR_ALL_PCA_4.pdf", width = 6, height = 4)
PCElbowPlot(object = seuObjFiltNormPC, num.pc = 50)
dev.off()

## determine statistically significant principal components
seuObjFiltNormPC <- JackStraw(object = seuObjFiltNormPC, num.pc = 50, num.replicate = 100)

pdf("STR_ALL_PCA_5.pdf", width = 12, height = 8)
JackStrawPlot(object = seuObjFiltNormPC, PCs = 1:50)
dev.off()

## save PCA analysis
save(seuObjFiltNormPC, file = "STR_SEURAT_DATA_FILT_NORM_PCA.RData")
```



#### DATA CLUTERING

Using the significant number of PCs, first the number of clustering is identified using Lovain clustering approach as part of the Seurat pipeline followed by clustering the cells into tSNE and UMAP space. Clustering plots are also generated. Updated seurat object is saved.

```{r}
## run Lovain clustering algorithm on statistically significant PCs to identify clusters
numPCs <- 50 ## based on Jackstraw plot generated in earlier step
seuObjFiltNormPCClust <- FindClusters(object = seuObjFiltNormPC, reduction.type = "pca", dims.use = 1:numPCs, resolution = 1.6, print.output = 0, save.SNN = TRUE, force.recalc = TRUE) ## tweak the values of resolution for finetuned clustering

## print the detected clusters and number of cells assigned to clusters
table(seuObjFiltNormPCClust@ident)

## run non-linear dimensional reduction (tSNE)
seuObjFiltNormPCClust <- RunTSNE(object = seuObjFiltNormPCClust, dims.use = 1:numPCs, do.fast = TRUE)

pdf("STR_ALL_CLUST_TSNE.pdf", width = 7, height = 6)
TSNEPlot(object = seuObjFiltNormPCClust, do.label = TRUE, pt.size = 1)
dev.off()

## run UMAP
seuObjFiltNormPCClust <- RunUMAP(seuObjFiltNormPCClust, reduction.use = "pca", dims.use = 1:numPCs)

umap1 <- DimPlot(seuObjFiltNormPCClust, reduction.use = "umap", do.return = TRUE, do.label = TRUE)
ggsave("STR_ALL_CLUST_UMAP.pdf", plot = umap1, width = 7, height = 6, units = "in", dpi = 300)

## save post-clustering seurat object
save(seuObjFiltNormPCClust, file = "STR_SEURAT_DATA_FILT_NORM_PCA_CLUST.RData")
```



#### CLUSTER MARKERS

Run DEG test to identify cluster marker genes to detect the cell types.

```{r}
## find markers for every cluster compared to all remaining cells
str.markers <- FindAllMarkers(object = seuObjFiltNormPCClust, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

write.table(str.markers, "STR_ALL_DEG_TABLE.txt", row.names = T, col.names = T, quote = F, sep = "\t")
save(str.markers, file = "STR_SEURAT_DATA_FILT_NORM_PCA_CLUST_DEG.RData")
```



####SESSION INFO

```{r}
sessionInfo()

R version 3.4.1 (2017-06-30)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Red Hat Enterprise Linux Server 7.4 (Maipo)

Matrix products: default
BLAS/LAPACK: /cm/shared/apps/intel/compilers_and_libraries/2017.6.256/linux/mkl/lib/intel64_lin/libmkl_gf_lp64.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] reshape2_1.4.3 gridExtra_2.3  dplyr_0.7.7    Seurat_2.3.4   Matrix_1.2-14 
[6] cowplot_0.9.3  ggplot2_3.1.0 

loaded via a namespace (and not attached):
  [1] tsne_0.1-3          segmented_0.5-3.0   nlme_3.1-131       
  [4] bitops_1.0-6        bit64_0.9-7         httr_1.3.1         
  [7] RColorBrewer_1.1-2  prabclus_2.2-6      tools_3.4.1        
 [10] backports_1.1.2     irlba_2.3.2         R6_2.3.0           
 [13] rpart_4.1-11        KernSmooth_2.23-15  Hmisc_4.1-1        
 [16] lazyeval_0.2.1      colorspace_1.3-2    trimcluster_0.1-2.1
 [19] nnet_7.3-12         withr_2.1.2         tidyselect_0.2.5   
 [22] bit_1.1-14          compiler_3.4.1      htmlTable_1.12     
 [25] hdf5r_1.0.0         diptest_0.75-7      caTools_1.17.1.1   
 [28] scales_0.5.0        checkmate_1.8.5     lmtest_0.9-36      
 [31] DEoptimR_1.0-8      mvtnorm_1.0-8       robustbase_0.93-2  
 [34] ggridges_0.5.0      pbapply_1.3-4       dtw_1.20-1         
 [37] proxy_0.4-22        stringr_1.3.1       digest_0.6.18      
 [40] mixtools_1.1.0      foreign_0.8-69      R.utils_2.7.0      
 [43] base64enc_0.1-3     pkgconfig_2.0.2     htmltools_0.3.6    
 [46] bibtex_0.4.2        htmlwidgets_1.2     rlang_0.3.0.1      
 [49] rstudioapi_0.8      bindr_0.1.1         jsonlite_1.5       
 [52] zoo_1.8-4           ica_1.0-2           mclust_5.4.1       
 [55] gtools_3.8.1        acepack_1.4.1       R.oo_1.22.0        
 [58] magrittr_1.5        modeltools_0.2-22   Formula_1.2-3      
 [61] lars_1.2            Rcpp_0.12.19        munsell_0.5.0      
 [64] reticulate_1.10     ape_5.2             R.methodsS3_1.7.1  
 [67] stringi_1.2.4       gbRd_0.4-11         MASS_7.3-47        
 [70] flexmix_2.3-14      gplots_3.0.1        Rtsne_0.13         
 [73] plyr_1.8.4          grid_3.4.1          parallel_3.4.1     
 [76] gdata_2.18.0        crayon_1.3.4        doSNOW_1.0.16      
 [79] lattice_0.20-35     splines_3.4.1       SDMTools_1.1-221   
 [82] knitr_1.20          pillar_1.3.0        igraph_1.2.2       
 [85] fpc_2.1-11.1        codetools_0.2-15    stats4_3.4.1       
 [88] glue_1.3.0          metap_1.0           latticeExtra_0.6-28
 [91] data.table_1.11.8   png_0.1-7           Rdpack_0.9-0       
 [94] foreach_1.4.4       tidyr_0.8.1         gtable_0.2.0       
 [97] RANN_2.6            purrr_0.2.5         kernlab_0.9-26     
[100] assertthat_0.2.0    class_7.3-14        survival_2.41-3    
[103] tibble_1.4.2        snow_0.4-3          iterators_1.0.10   
[106] bindrcpp_0.2.2      cluster_2.0.6       fitdistrplus_1.0-9 
[109] ROCR_1.0-7
```



Last updated: 4/17/2019.

