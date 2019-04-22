# DATA QC

Starting with the collapsed and genotypes combined UMI counts table, this script performs basic data quality check using functions available as part of *Seurat* R package v2.3.4.



#### CREATE SEURAT OBJECT

Load *RData* file containing a collapsed and genotypes combined data frame for UMI counts across all libraries/samples and their technical replicates followed by creating a seurat object without filtering for anything yet. Also update meta data within seurat object such as genotypes, library information.

```{r}
## load collapsed and genotypes combined data frame
load("STR_Collapsed_Combined_Genotypes.RData")

## create seurat object
seuObjAll <- CreateSeuratObject(raw.data = allData, project = "STR_ALL")

## extract meta data from rownames
metaTemp <- as.data.frame(matrix(unlist(strsplit(rownames(seuObjAll@meta.data), "_")), ncol = 4, byrow = TRUE))
row.names(metaTemp) <- row.names(seuObjAll@meta.data)
colnames(metaTemp) <- c("Source", "Genotype", "Library", "CellBarcode")

## update genotypes in seurat object meta data
geno <- metaTemp$Genotype
names(geno) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = geno, col.name = "Genotype")

## update library information in seurat object meta data
lib <- metaTemp$Library
names(lib) <- row.names(metaTemp)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = lib, col.name = "Library")

## save seurat object
save(seuObjAll, file = "STR_SEURAT_DATA.RData")
```



#### UMI HISTOGRAM

Seurat object can now be used to create violin plot or box plot or histogram for number of UMIs and number of genes.

```{r}
## violin plot for number of UMIs across genotype
plotVio1 <- ggplot(seuObjAll@meta.data, aes(x = factor(Genotype), y = nUMI, fill = Genotype)) + geom_violin(scale = "width", trim = FALSE) + labs(title = paste(seuObjAll@project.name, "nUMI", sep = " | "), x = "Libraries", y = "Number of UMIs") + theme_bw() + theme(legend.position = "none")
ggsave(paste(seuObjAll@project.name, "_UMI_VIOLIN.pdf", sep = ""), plot = plotVio1, width = 6, height = 6, units = "in", dpi = 150)	

## violin plot for number of genes across genotype
plotVio2 <- ggplot(seuObjAll@meta.data, aes(x = factor(Genotype), y = nGene, fill = Genotype)) + geom_violin(scale = "width", trim = FALSE) + labs(title = paste(seuObjAll@project.name, "nGenes", sep = " | "), x = "Libraries", y = "Number of Genes") + theme_bw() + theme(legend.position = "none")
ggsave(paste(seuObjAll@project.name, "_GENES_VIOLIN.pdf", sep = ""), plot = plotVio2, width = 6, height = 6, units = "in", dpi = 150)	

## histogram for number of UMIs across libraries
plotHist1 <- ggplot(seuObjAll@meta.data, aes(x = log10(nUMI), fill = Library, alhpa = 0.5)) + geom_histogram(bins = 100) + labs(title = paste(seuObjAll@project.name, "nUMI", sep = " | "), x = "Log10(Number of UMIs)", y = "Number of Cells") + theme_bw() + scale_fill_discrete(name = "Libraries")
ggsave(paste(seuObjAll@project.name, "_LIB_HIST.pdf", sep = ""), plot = plotHist1, width = 7, height = 6, units = "in", dpi = 150)

## histogram for number of UMIs across all cells
plotHist2 <- ggplot(seuObjAll@meta.data, aes(x = log10(nUMI))) + geom_histogram(bins = 100) + labs(title = paste(seuObjAll@project.name, "nUMI", sep = " | "), x = "Log10(Number of UMIs)", y = "Number of Cells") + theme_bw()
ggsave(paste(seuObjAll@project.name, "_HIST.pdf", sep = ""), plot = plotHist2, width = 6, height = 6, units = "in", dpi = 150)
```



#### CALCULATE PERCENT MITOCHONDRIAL CONTENT

Further, UMI counts data across all genes is used to identify mitochondrial genes to calculate percent mitochondrial content per cell and meta data of seurat object is also updated accordingly. Now seurat object can be used to create a violin plot for percent mitochondrial content.

```{r}
## identify mitochondrial genes and calculate percent mitochondrial content per cell
mito.genes <- grep(pattern = "^mt-", x = rownames(x = seuObjAll@data), value = TRUE)
percent.mito <- Matrix::colSums(seuObjAll@data[mito.genes, ])/Matrix::colSums(seuObjAll@data)
seuObjAll <- AddMetaData(object = seuObjAll, metadata = percent.mito, col.name = "pMito")

## violin plot for percent mitochondrial content across genotypes
plotVio3 <- ggplot(seuObjAll@meta.data, aes(x = factor(Genotype), y = pMito, fill = Genotype)) + geom_violin(scale = "width", trim = FALSE) + labs(title = paste(seuObjAll@project.name, "pMito", sep = " | "), x = "Libraries", y = "Percent Mito") + theme_bw() + theme(legend.position = "none")  + scale_y_continuous(labels = scales::percent)
ggsave(paste(seuObjAll@project.name, "_MITO_VIOLIN.pdf", sep = ""), plot = plotVio3, width = 6, height = 6, units = "in", dpi = 150)
```



#### ADDITITIONAL QC PLOTS

Plotting percent mitochondrial content and number of genes per cell across number of UMIs can further facilitate filtering.

```{r}
## violin plots for number of UMIs, genes and percent mitochondrial content across all cells
plot_nU <- VlnPlot(object = seuObjAll, features.plot = "nUMI", point.size.use = 0) + theme(legend.position = "none")
plot_nG <- VlnPlot(object = seuObjAll, features.plot = "nGene", point.size.use = 0) + theme(legend.position = "none")
plot_pM <- VlnPlot(object = seuObjAll, features.plot = "pMito", point.size.use = 0) + theme(legend.position = "none") + scale_y_continuous(labels = scales::percent)
plotQC1 <- grid.arrange(plot_nU, plot_nG, plot_pM, ncol = 3)
ggsave(paste(seuObjAll@project.name, "_QC_1.pdf", sep = ""), plot = plotQC1, width = 12, height = 4, units = "in", dpi = 300)

## plotting pMito against nUMI & nGenes against nUMI
pdf(paste(seuObjAll@project.name, "_QC_2.pdf", sep = ""), width = 12, height = 6)
par(mfrow = c(1, 2))
GenePlot(object = seuObjAll, gene1 = "nUMI", gene2 = "pMito", cex.use = 1)
GenePlot(object = seuObjAll, gene1 = "nUMI", gene2 = "nGene", cex.use = 1)
dev.off()
```



#### SAVE UPDATED SEURAT OBJECT

```{r}
save(seuObjAll, file = "STR_SEURAT_DATA_QC.RData")
```



#### R SESSION INFO

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
[1] gridExtra_2.3 Seurat_2.3.4  Matrix_1.2-14 cowplot_0.9.3 ggplot2_3.1.0

loaded via a namespace (and not attached):
  [1] tsne_0.1-3          segmented_0.5-3.0   nlme_3.1-131       
  [4] bitops_1.0-6        bit64_0.9-7         httr_1.3.1         
  [7] RColorBrewer_1.1-2  prabclus_2.2-6      tools_3.4.1        
 [10] backports_1.1.2     irlba_2.3.2         R6_2.3.0           
 [13] rpart_4.1-11        KernSmooth_2.23-15  Hmisc_4.1-1        
 [16] lazyeval_0.2.1      colorspace_1.3-2    trimcluster_0.1-2.1
 [19] nnet_7.3-12         withr_2.1.2         tidyselect_0.2.5   
 [22] bit_1.1-14          compiler_3.4.1      htmlTable_1.12     
 [25] hdf5r_1.0.0         labeling_0.3        diptest_0.75-7     
 [28] caTools_1.17.1.1    scales_0.5.0        checkmate_1.8.5    
 [31] lmtest_0.9-36       DEoptimR_1.0-8      mvtnorm_1.0-8      
 [34] robustbase_0.93-2   ggridges_0.5.0      pbapply_1.3-4      
 [37] dtw_1.20-1          proxy_0.4-22        stringr_1.3.1      
 [40] digest_0.6.18       mixtools_1.1.0      foreign_0.8-69     
 [43] R.utils_2.7.0       base64enc_0.1-3     pkgconfig_2.0.2    
 [46] htmltools_0.3.6     bibtex_0.4.2        htmlwidgets_1.2    
 [49] rlang_0.3.0.1       rstudioapi_0.8      bindr_0.1.1        
 [52] jsonlite_1.5        zoo_1.8-4           ica_1.0-2          
 [55] mclust_5.4.1        gtools_3.8.1        acepack_1.4.1      
 [58] dplyr_0.7.7         R.oo_1.22.0         magrittr_1.5       
 [61] modeltools_0.2-22   Formula_1.2-3       lars_1.2           
 [64] Rcpp_0.12.19        munsell_0.5.0       reticulate_1.10    
 [67] ape_5.2             R.methodsS3_1.7.1   stringi_1.2.4      
 [70] gbRd_0.4-11         MASS_7.3-47         flexmix_2.3-14     
 [73] gplots_3.0.1        Rtsne_0.13          plyr_1.8.4         
 [76] grid_3.4.1          parallel_3.4.1      gdata_2.18.0       
 [79] crayon_1.3.4        doSNOW_1.0.16       lattice_0.20-35    
 [82] splines_3.4.1       SDMTools_1.1-221    knitr_1.20         
 [85] pillar_1.3.0        igraph_1.2.2        fpc_2.1-11.1       
 [88] reshape2_1.4.3      codetools_0.2-15    stats4_3.4.1       
 [91] glue_1.3.0          metap_1.0           latticeExtra_0.6-28
 [94] data.table_1.11.8   png_0.1-7           Rdpack_0.9-0       
 [97] foreach_1.4.4       tidyr_0.8.1         gtable_0.2.0       
[100] RANN_2.6            purrr_0.2.5         kernlab_0.9-26     
[103] assertthat_0.2.0    class_7.3-14        survival_2.41-3    
[106] tibble_1.4.2        snow_0.4-3          iterators_1.0.10   
[109] bindrcpp_0.2.2      cluster_2.0.6       fitdistrplus_1.0-9 
[112] ROCR_1.0-7
```



Last updated: 04/18/2019.

