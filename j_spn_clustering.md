# SPN CLUSTERING

This script performs normalization, scaling on subset of data corresponding to neuronal cells from primary dataset which is already processed through Seurat pipeline, filtered, normalized and clustered. Since the primary dataset is already filtered after QC, no filtered is needed here. Since, data being used for subset is raw UMI count data, normalization and scaling is required.



####DATA NORMALIZATION

```{r}
## load seurat object created after data subset
load("STR_SPN_SEURAT_DATA.RData")

# log-normalize the data
aaxNorm <- NormalizeData(object = aax, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify ~2000 variable genes
pdf("AAX_SELECTED_VAR_GENES_woXYMT.pdf", width = 8, height = 6)
aaxNorm <- FindVariableGenes(object = aaxNorm, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.2, x.high.cutoff = 2.5, y.cutoff = 0.5)
dev.off()

# scale the data by regressing out for un-intended variation (nUMI per cell, batch effects, cell cycle, mito.genes, etc.)
aaxNorm <- ScaleData(object = aaxNorm, vars.to.regress = c("nUMI", "pMito"))

# save RData file
save(aaxNorm, file = "STR_SPN_SEURAT_DATA_NORM.RData")
```



#### PCA ANALYSIS

```{r}
# PCA using identified variable genes on regressed (scaled) data
aaxPC <- RunPCA(object = aaxNorm, pc.genes = aaxNorm@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)

# examine and visualize PCA results a few different ways
pdf("STR_SPN_SEURAT_PCA_1.pdf", width = 9, height = 21)
VizPCA(object = aaxPC, pcs.use = 1:50)
dev.off()

pdf("STR_SPN_SEURAT_PCA_2.pdf", width = 7, height = 6)
PCAPlot(object = aaxPC, dim.1 = 1, dim.2 = 2)
dev.off()

aaxPC <- ProjectPCA(object = aaxPC, do.print = FALSE)

pdf("STR_SPN_SEURAT_PCA_3.pdf", width = 9, height = 21)
PCHeatmap(object = aaxPC, pc.use = 1:50, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

pdf("STR_SPN_SEURAT_PCA_4.pdf", width = 6, height = 4)
PCElbowPlot(object = aaxPC, num.pc = 50)
dev.off()

# determine statistically significant principal components, takes long time
aaxPC <- JackStraw(object = aaxPC, num.pc = 50, num.replicate = 100)

pdf("STR_SPN_SEURAT_PCA_5.pdf", width = 12, height = 21)
JackStrawPlot(object = aaxPC, PCs = 1:50)
dev.off()

## save RData file
save(aaxPC, file = "STR_SPN_SEURAT_DATA_NORM_PCA.RData")
```



####DATA CLUSTERING

```{r}
## based on PCA analysis, identify the number of PCs for clustering
numPCs <- 42
aaxNormClust <- FindClusters(object = aaxPC, reduction.type = "pca", dims.use = 1:numPCs, resolution = 1.4, print.output = 0, save.SNN = TRUE)

## print the number of detected clusters and number of cell per cluster
table(aaxNormClust@ident)

# clustering using tSNE
aaxNormClust <- RunTSNE(object = aaxNormClust, dims.use = 1:numPCs, do.fast = TRUE, perplexity = 50)

pdf("STR_SPN_SEURAT_CLUST_TSNE.pdf", width = 7, height = 6)
TSNEPlot(object = aaxNormClust, do.label = TRUE, pt.size = 0.5)
dev.off()

## clustering using UMAP
aaxNormClust <- RunUMAP(aaxNormClust, reduction.use = "pca", dims.use = 1:numPCs)

umap1 <- DimPlot(aaxNormClust, reduction.use = "umap", do.return = TRUE, do.label = TRUE)
ggsave("STR_SPN_SEURAT_CLUST_UMAP.pdf", plot = umap1, width = 7, height = 6, units = "in", dpi = 300)

save(aaxNormClust, file = "STR_SPN_SEURAT_DATA_NORM_PCA_CLUST.RData")
```



####DEG ANALYSIS

```{r}
# find markers for every cluster compared to all remaining cells
aax.markers <- FindAllMarkers(object = aaxNormClust, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

write.table(aax.markers, "STR_SPN_SEURAT_DATA_NORM_PCA_CLUST_DEG_TABLE.txt", row.names = T, col.names = T, quote = F, sep = "\t")
save(aax.markers, file = "STR_SPN_SEURAT_DATA_NORM_PCA_CLUST_DEG.RData")
```



#### SESSION INFO

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

