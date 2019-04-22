# ADDITIONAL CODE

In addition to main analysis, code chunks for few additional analyses, plots is covered here.



#### EXPORT RAW DATA WITH META INFOTMATION

Raw UMI count data across all cells and their corresponding clusters and addtional meta information is exported for further use such as creating a subset of data for neuronal cells.

```{r}
## load seurat object with UMAP clustering information
load("STR_SEURAT_DATA_FILT_NORM_PCA_CLUST.RData")

data.temp <- as.matrix(seuObjFiltNormPCClust@raw.data)
expData <- as.data.frame(t(data.temp))

metaTemp <- as.data.frame(seuObjFiltNormPCClust@meta.data)
colnames(metaTemp) <- c("nGenes", "nUMI", "Genotype", "Library", "pMito", "Cluster")
metaTemp$NewGeno <- rep(0, nrow(metaTemp))
metaTemp[metaTemp$Genotype == "CTL","NewGeno"] <- 1
metaTemp[metaTemp$Genotype == "D1CKO","NewGeno"] <- 2
metaTemp[metaTemp$Genotype == "D2CKO","NewGeno"] <- 3
metaTemp[metaTemp$Genotype == "DDCKO","NewGeno"] <- 4
metaTemp$Cluster <- sprintf("%02d", as.numeric(metaTemp$Cluster))

umapProj <- seuObjFiltNormPCClust@dr$umap@cell.embeddings
colnames(umapProj) <- c("UMAP_1", "UMAP_2")

metaData <- merge(umapProj, metaTemp, by = "row.names", all = TRUE)
row.names(metaData) <- metaData$Row.names
metaData$Row.names <- NULL

allData <- merge(metaData, expData, by = "row.names")
row.names(allData) <- allData$Row.names
colnames(allData)[colnames(allData) == "Row.names"] <- "CellBarcode"

## save RData file
save(allData, file = "STR_SEURAT_DATA_FILT_NORM_PCA_CLUST_RAWDATA.RData")
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

