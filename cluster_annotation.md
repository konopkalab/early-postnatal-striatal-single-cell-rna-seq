# ANNOTATING THE CLUSTERS

This script covers the two approaches used to annotate and classify cells into different cell-types. Annotation is based on a published adult mouse striatal single cell dataset. Genes corresponding to each identified cluster for published dataset were identified using Seurat pipeline followed by (i) custom developed script to perform hypergeomtric test for overlap with cell-type cluster marker genes from published dataset and (ii) expression weighted cell-type enrichment (EWCE) analysis with the same published dataset as reference.



#### IDENTIFICATION OF CLUSTER MARKERS FOR PUBLISHED REFERENCE DATASETS

The published reference dataset (Saunders et al., 2018, Cell 174, 1015–1030, http://dropviz.org) is composed of cell types from 9 regions of adult mouse brain. Here, only striatal dataset was used to identify cluster-specific gene markers using Seurat (v2.3.4) pipeline without filtering the dataset and keeping the original published cell-type annotation intact.

```{r}
## load reference dataset along with cluster/cell-type annotation provided
dge <- loadSparseDge("F_GRCm38.81.P60Striatum.raw.dge.txt.gz")
clusterOutcome <- readRDS("F_GRCm38.81.P60Striatum.cell_cluster_outcomes.RDS")
clusterAssign <- readRDS("F_GRCm38.81.P60Striatum.cluster.assign.RDS")
subclusterAssign <- readRDS("F_GRCm38.81.P60Striatum.subcluster.assign.RDS")

## load striatal subset of meta data provided
metaSTR <- read.table("metaSTR2.txt", header = T, sep = "\t")

## create a subset digital gene expression dataset matching to striatal meta data
expData <- as.data.frame(as.matrix(dge))
cells2useData <- clusterOutcome[clusterOutcome$cluster %in% unique(sort((metaSTR$cluster))),]
saunders.data <- expData[,colnames(expData) %in% row.names(cells2useData)]

## create seurat object and update metadata with available cell-type association
saunders <- CreateSeuratObject(raw.data = saunders.data, project = "SAUNDERS_STR")

saunders.meta.cluster <- cells2useData$cluster
names(saunders.meta.cluster) <- row.names(cells2useData)
saunders <- AddMetaData(object = saunders, metadata = saunders.meta.cluster, col.name = "Cluster")

saunders.meta.subcluster <- cells2useData$subcluster
names(saunders.meta.subcluster) <- row.names(cells2useData)
saunders <- AddMetaData(object = saunders, metadata = saunders.meta.subcluster, col.name = "SubCluster")

## save seurat object for reference striatal dataset
save(saunders, saunders.data, cells2useData, file = "SAUNDERS_STR_DATA.RData")

## normalize and scale the data
saundersNorm <- NormalizeData(object = saunders, normalization.method = "LogNormalize", scale.factor = 10000)

## identify ~2000 variable genes
pdf("SAUNDERS_STR_VAR_GENES.pdf", width = 8, height = 6)
saundersNorm <- FindVariableGenes(object = saundersNorm, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.3, x.high.cutoff = 4, y.cutoff = 0.5)
dev.off()

## scale the data by regressing for number of UMIs per cell
saundersNorm <- ScaleData(object = saundersNorm, vars.to.regress = "nUMI")

## save normalized, scaled and regressed dataset
save(saundersNorm, file = "SAUNDERS_STR_NORM.RData")

## replace the identities in seurat object by published cluster identities
saundersNorm@ident <- saundersNorm@meta.data$Cluster
names(saundersNorm@ident) <- row.names(saundersNorm@meta.data)
saundersNorm@meta.data$res.0.8 <- saundersNorm@meta.data$Cluster

## identify gene markers per cluster and save the data
saunders.markers <- FindAllMarkers(object = saundersNorm, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
write.table(saunders.markers, "SAUNDERS_STR_DEG_TABLE.txt", row.names = T, col.names = T, quote = F, sep = "\t")
save(saunders.markers, file = "SAUNDERS_STR_DEG.RData")

## merge published annotation with cluster-specific gene markers dataset
## this file can further be used to identify cell-types in our striatal dataset
CellTypesPerCluster <- read.table("SaundersSTR_CellTypeClusterAnnotation.txt", header = T, sep = "\t")
degSaunders <- merge(saunders.markers, CellTypesPerCluster, by = "cluster")
write.table(degSaunders, "SAUNDERS_STR_DEG_TABLE_ANNOTATED.txt", row.names = F, col.names = T, quote = F, sep = "\t")

## session info
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
[1] DropSeq.util_2.0  data.table_1.11.8 gridExtra_2.3     dplyr_0.7.7      
[5] Seurat_2.3.4      Matrix_1.2-14     cowplot_0.9.3     ggplot2_3.1.0    

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
 [85] fpc_2.1-11.1        reshape2_1.4.3      codetools_0.2-15   
 [88] stats4_3.4.1        glue_1.3.0          metap_1.0          
 [91] latticeExtra_0.6-28 png_0.1-7           Rdpack_0.9-0       
 [94] foreach_1.4.4       tidyr_0.8.1         gtable_0.2.0       
 [97] RANN_2.6            purrr_0.2.5         kernlab_0.9-26     
[100] assertthat_0.2.0    class_7.3-14        survival_2.41-3    
[103] tibble_1.4.2        snow_0.4-3          iterators_1.0.10   
[106] bindrcpp_0.2.2      cluster_2.0.6       fitdistrplus_1.0-9 
[109] ROCR_1.0-7
```



#### PREPARE REFERENCE DEG TABLE FOR HYPERGEOMETRIC TEST

Hypergeometric test used in the analysis requires the reference dataset to be in a specific format.

```{r}
## load DEG table generated by Seurat DEG analysis with cluster and cell-type association
saundersDEG <- read.table("SAUNDERS_STR_DEG_TABLE_ANNOTATED.txt", header = T, sep = "\t")

## filter DEG table
saundersSIGt <- saundersDEG[saundersDEG$p_val_adj <= 0.05,]
saundersSIG <- saundersSIGt[saundersSIGt$pct.1 >= 0.5 | saundersSIGt$pct.2 >= 0.5,]

## create list of data frames for each cell-type & associated gene markers
degSaunders <- vector("list", length(unique(sort(saundersSIG$celltype))))
names(degSaunders) <- unique(sort(saundersSIG$celltype))

for (i in 1:length(unique(sort(saundersSIG$celltype))))
  {
  print(paste("Extracting Genes for", unique(sort(saundersSIG$celltype))[i] , sep = " "))
  tempDEG <- saundersSIG[saundersSIG$celltype == unique(sort(saundersSIG$celltype))[i],]
  filtDEG <- tempDEG[,c("gene", "celltype")]
  degSaunders[[i]] <- filtDEG
  }

## save deg dataset
save(degSaunders, file = "Reference_DEG_for_Hypergeometric.RData")


## session info
sessionInfo()

R version 3.4.1 (2017-06-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.3

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] RColorBrewer_1.1-2  pheatmap_1.0.12     ggplot2_3.1.0       reshape2_1.4.3      gmp_0.5-13.2       
[6] multtest_2.38.0     Biobase_2.42.0      BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0       compiler_3.5.1   pillar_1.3.1     plyr_1.8.4       bindr_0.1.1      tools_3.5.1     
 [7] tibble_2.0.1     gtable_0.2.0     lattice_0.20-38  pkgconfig_2.0.2  rlang_0.3.1      Matrix_1.2-15   
[13] rstudioapi_0.9.0 bindrcpp_0.2.2   withr_2.1.2      stringr_1.3.1    dplyr_0.7.8      stats4_3.5.1    
[19] tidyselect_0.2.5 glue_1.3.0       R6_2.3.0         survival_2.43-3  purrr_0.3.0      magrittr_1.5    
[25] scales_1.0.0     MASS_7.3-51.1    splines_3.5.1    assertthat_0.2.0 colorspace_1.4-0 stringi_1.2.4   
[31] lazyeval_0.2.1   munsell_0.5.0    crayon_1.3.4 
```



#### HYPERGEOMETRIC TEST FOR OVERLAP

This part of the code perform overlap between gene markers for every cluster in our dataset against gene markers for every cell-type associated gene markers in reference dataset. Based on the number of overlapping genes, hypergeometric test assigns a p-value for the overlap using number of genes expressed in our dataset as background. P-values are further corrected using Benjamini–Hochberg procedure to yield adjusted p-values which are used to generate a heatmap to visualize the enrichment of gene markers in our cluster against reference cell-types. Thus our clusters are also annotated based on the significant enrichment across reference cell-types.

```{r}
## custom function to perform hypergeometric test
enrich_pvalue <- function(N, A, B, k)
 {
 m <- A + k
 n <- B + k
 i <- k:min(m,n)
 as.numeric( sum(chooseZ(m,i)*chooseZ(N-m,n-i))/chooseZ(N,n) )
 }

## load DEG table generated by Seurat pipeline for cluser-specific markers for our dataset
degtable <- list.files(path = ".", pattern = 'DEG_TABLE')
exp <- read.table(paste("./", degtable, sep = ""), sep="\t", header=T)
tab <- exp[exp$p_val_adj <= 0.05,] ## filter DEG table
tab <- tab[tab$pct.1 >= 0.75 | degDataSigt$pct.2 >= 0.75,] ## filter DEG table
tab <- tab[c(7,6)] ## use only selected columns for gene name and cluster association
tab$cluster <- paste("Cluster_",tab$cluster,sep="")
colnames(tab) <- c("Gene","DEFINITION")
tab$DEFINITION <- factor(tab$DEFINITION, levels = c(paste("Cluster", unique(sort(exp$cluster)), sep = "_")))
Genes <- as.data.frame(table(tab$DEFINITION))

## load a list of data frames having gene name and associated cell types generated above
load("Reference_DEG_for_Hypergeometric.RData")
GeneSets <- degSaunders

## perform overlaps and save overlapping data and number of overlapping genes
ln <- length(GeneSets)
cl <- length(Genes$Var1)
TEMP <- list()
INT <- list()
for (i in 1:ln)
 {
 TEMP[[i]] <- tab[tab$Gene %in% GeneSets[[i]]$gene,]
 INT[[i]] <- as.data.frame(table(TEMP[[i]]$DEFINITION))
 }
names(INT) <- names(GeneSets)
names(TEMP) <- names(GeneSets)
save(TEMP, file <- "GeneSets_Intersections.RData")
save(INT, file <- "GeneSets_Intersections_Number.RData")

## perform hypergeometric test using custom function above
PVAL_Brain <- data.frame()
for (i in 1:cl)
 {
 for (j in 1:ln)
  {
  PVAL_Brain[i,j] <- enrich_pvalue(7500,  Genes[i,2] - INT[[j]][i,2], length(GeneSets[[j]]$gene) - INT[[j]][i,2], INT[[j]][i,2]) ## using background as number of expressed genes in our dataset
  }
 }
names(PVAL_Brain) <- names(GeneSets)
rownames(PVAL_Brain) <- Genes$Var1

## calculate fdr or adjusted p-value for overlaps using Benjamini–Hochberg procedure
ADJ_Brain <- apply(PVAL_Brain, 2, function(x){p.adjust(x,"BH")})

## save p-values and adjusted p-value tables
write.table(PVAL_Brain, "ENRICH_TABLES/Pvalue_BrainExp_GeneSets_overlap.txt", sep="\t", quote=F)
write.table(ADJ_Brain, "ENRICH_TABLES/Adj_BrainExp_GeneSets_overlap.txt", sep="\t", quote=F)

# generate heatmap for overlap associated adjusted p-value
ADJ_Brain[ADJ_Brain > 0.05] <- 1
LOG <- as.data.frame(-log10(ADJ_Brain))
LOG$DEFINITION <- rownames(LOG)
rownames(LOG) <- NULL
df <- melt(LOG)
colnames(df) <- c("DEFINITION","Assoc","Pval")
colfunc <- colorRampPalette(c("white", "red2"))
df$DEFINITION <- factor(df$DEFINITION, levels = c(paste("Cluster", unique(sort(exp$cluster)), sep = "_")))

heatmap1 <- ggplot(df, aes(x = DEFINITION, y = Assoc, fill = Pval)) + 
									geom_tile(colour = "black")+geom_text(label = signif(df$Pval,3), size = 2.7,fontface="bold") + 
									scale_fill_gradientn(colours = colfunc(30),name = "-Log10(P)") + 
                  scale_x_discrete(expand = c(0, 0)) + 
                  scale_y_discrete(expand = c(0, 0)) + 
                  coord_equal() + 
                  theme_bw() + 
                  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
                  theme(strip.text.x = element_text(size=30, color = "black",family="serif",face="bold"),
                        strip.background = element_rect(colour="black", fill="white"))+ 
									scale_y_discrete(name="") +
                  scale_x_discrete(name="") + 
                  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                  theme(panel.margin = unit(8, "lines")) +
                  theme(axis.text=element_text(color = "black",family="serif",face="bold",size=16)) + 
                  theme(panel.border=element_rect(color="black"))+theme(legend.key.size =  unit(0.2, "in")) + 
                  theme(legend.position= "right") + 
                  theme(legend.text = element_text(size=12,color = "black",family="serif",face="bold") +
                  theme(legend.title = element_text(size=20,color = "black",family="serif",face="bold")))
ggsave(filename="ENRICH_PLOTS/ADJ_Brain_HEATMAP_GeneSets.pdf", plot=heatmap1, width=24, height=12, units="in")

## session info
sessionInfo()

R version 3.4.1 (2017-06-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.3

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] RColorBrewer_1.1-2  pheatmap_1.0.12     ggplot2_3.1.0       reshape2_1.4.3      gmp_0.5-13.2       
[6] multtest_2.38.0     Biobase_2.42.0      BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.0       compiler_3.5.1   pillar_1.3.1     plyr_1.8.4       bindr_0.1.1      tools_3.5.1     
 [7] tibble_2.0.1     gtable_0.2.0     lattice_0.20-38  pkgconfig_2.0.2  rlang_0.3.1      Matrix_1.2-15   
[13] rstudioapi_0.9.0 bindrcpp_0.2.2   withr_2.1.2      stringr_1.3.1    dplyr_0.7.8      stats4_3.5.1    
[19] tidyselect_0.2.5 glue_1.3.0       R6_2.3.0         survival_2.43-3  purrr_0.3.0      magrittr_1.5    
[25] scales_1.0.0     MASS_7.3-51.1    splines_3.5.1    assertthat_0.2.0 colorspace_1.4-0 stringi_1.2.4   
[31] lazyeval_0.2.1   munsell_0.5.0    crayon_1.3.4
```



#### EXPRESSION WEIGHTED CELL_TYPE ENRICHMENT (EWCE) ANALYSIS

This is another method to annotate clusters with cell-types based on gene markers across our and reference dataset. This method is publicly available as a R-package (https://github.com/NathanSkene/EWCE). First, prepare the reference dataset for EWCE analysis followed by EWCE analysis using our dataset.

```{r}
## load reference dataset generated during Seurat analysis described above
load("SAUNDERS_STR_DATA.RData")
expData <- as.data.frame(t(saunders.data))

## load cell-type association
CellTypesPerCluster <- read.table("SaundersSTR_CellTypeClusterAnnotation.txt", header = T, sep = "\t")

## re-arrange meta data
cells2useData$cell_id <- row.names(cells2useData)
metaData <- merge(cells2useData, CellTypesPerCluster, by = "cluster")
row.names(metaData) <- metaData$cell_id
colnames(metaData) <- c("cluster", "subcluster", "reason", "cell_id", "level1class")
metaData$level2class <- metaData$cluster

## combine expression data with meta data
dataSaunders <- merge(metaData, expData, by = "row.names")
row.names(dataSaunders) <- dataSaunders$Row.names
dataSaunders$Row.names <- NULL
dataSaunders$reason <- NULL
saundersData <- as.data.frame(t(dataSaunders))

## save the intermediate dataset
save(saundersData, file = "SAUNDERS_STR_DATA2.RData")

## convert data to EWCE compatible format
## this part of the code is based on the EWCE R-package
mRNA <- saundersData
expr <- mRNA[6:nrow(mRNA), ]
cell_type <- unlist(mRNA[4,])
subcell_type = sprintf("%s_%s", cell_type, unlist(mRNA[5,]))
matching_cell_type <- as.character(unique(data.frame(cell_type = cell_type, subcell_type = subcell_type))[, "cell_type"])
expr2 <- as.numeric(as.matrix(expr))
expr3 <- matrix(as.numeric(as.matrix(expr)), nrow = nrow(expr), ncol = ncol(expr))
rownames(expr3) <- row.names(mRNA)[6:nrow(mRNA)]

count <- 0
for (sct in unique(subcell_type)) {
    count <- count + 1
    sub_expr <- expr3[, subcell_type == sct]
    sct_expr <- data.frame(temp = apply(sub_expr, 1, mean))
    colnames(sct_expr) = sct
    if (count == 1) {
        all_scts = sct_expr
    	}
    else {
        all_scts = cbind(all_scts, sct_expr)
    	}
	}
rownames(all_scts) = rownames(expr3)

keepGenes = rownames(all_scts)[apply(all_scts, 1, max) > 0.2]
all_scts = all_scts[keepGenes, ]
cTs = unique(cell_type)
geneList = unique(rownames(all_scts))
count = 0
for (gs in geneList) {
    count = count + 1
    exp1 = unlist(all_scts[gs, ])
    exp2 = exp1/sum(exp1)
    exp3 = data.frame(e = exp2, cell = matching_cell_type)
    exp4 = aggregate(exp3$e, sum, by = list(exp3$cell))
    exp5 = data.frame(t(exp4[, 2]))
    colnames(exp5) = as.character(exp4[, 1])
    rownames(exp5) = gs
    if (count == 1) {
        cell_dists = exp5
	    }
    else {
        cell_dists = rbind(cell_dists, exp5)
    	}
	}

scTs = unique(subcell_type)
geneList = unique(rownames(all_scts))
count = 0
for (gs in geneList) {
    count = count + 1
    exp1 = unlist(all_scts[gs, ])
    exp2 = data.frame(t(exp1/sum(exp1)))
    rownames(exp2) = gs
    if (count == 1) {
        subcell_dists = exp2
    	}
    else {
        subcell_dists = rbind(subcell_dists, exp2)
    	}
	}

celltype_data <- list(all_scts = all_scts, cell_dists = cell_dists, subcell_dists = subcell_dists)

## save EWCE compatible dataset
save(celltype_data, file = "SAUNDERS_STR_DATA3.RData")

## load reference mouse-human gene homolog data
data("mouse_to_human_homologs")

## EWCE functions to perform bootstrap enrichment test for overlap
getOverlapEWCE1 <- function(clusterGenes)
  {
  m2h <- unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
  mouse.hits <- unique(clusterGenes)
  mouse.bg <- unique(setdiff(m2h$MGI.symbol, mouse.hits))
  reps <- 1000
  full_results1 <- bootstrap.enrichment.test(sct_data=celltype_data, mouse.hits=mouse.hits, mouse.bg=mouse.bg, reps=reps, sub=TRUE)
  return(full_results1$results)
  }

getOverlapEWCE2 <- function(clusterGenes)
  {
  m2h <- unique(mouse_to_human_homologs[,c("HGNC.symbol","MGI.symbol")])
  mouse.hits <- unique(clusterGenes)
  mouse.bg <- unique(setdiff(m2h$MGI.symbol, mouse.hits))
  reps=1000
  full_results2 <- bootstrap.enrichment.test(sct_data=celltype_data, mouse.hits=mouse.hits, mouse.bg=mouse.bg, reps=reps, sub=FALSE)
  return(full_results2$results)
  }


## load the DEG table for our dataset from Seurat pipeline
degData <- read.table("STR_ALL_DEG_TABLE.txt", sep = "\t", header = T)
degDataSigt <- degData[degData$p_val_adj <= 0.05,]
degDataSig <- degDataSigt[degDataSigt$pct.1 >= 0.75 | degDataSigt$pct.2 >= 0.75,]

## perform EWCE analysis
minorAnnotation <- vector("list", length(unique(sort(degDataSig$cluster))))
names(minorAnnotation) <- paste("Cluster", unique(sort(degDataSig$cluster)), sep = "_")
majorAnnotation <- vector("list", length(unique(sort(degDataSig$cluster))))
names(majorAnnotation) <- paste("Cluster", unique(sort(degDataSig$cluster)), sep = "_")

for(i in 0:max(unique(sort(degDataSig$cluster))))
	{
	print(paste("Processing Cluster ", i, sep = ""))
	degDataSigClu <- degDataSig[degDataSig$cluster == i,]
	degDataSigClu2 <- degDataSigClu[order(degDataSigClu$avg_logFC, decreasing = T),]
	deGenes <- degDataSigClu2$gene
	j <- i + 1

	ewceMinor <- getOverlapEWCE1(deGenes)
	ewceMinor1 <- ewceMinor[order(ewceMinor$fold_change, decreasing = T),]
	minorAnnotation[[j]] <- ewceMinor1
	plotMinor <- ewce.plot(minorAnnotation[[j]], mtc_method="BH")
	# plotMinor + labs(title = "Minor Cell Types") + theme(plot.margin = unit(c(1,1,1.5,1.2), "cm"))

	ewceMajor <- getOverlapEWCE2(deGenes)
	ewceMajor1 <- ewceMajor[order(ewceMajor$fold_change, decreasing = T),]
	majorAnnotation[[j]] <- ewceMajor1
	plotMajor <- ewce.plot(majorAnnotation[[j]], mtc_method="BH")
	# plotMajor + labs(title = "Major Cell Types") + theme(plot.margin = unit(c(1,1,1.5,1.2), "cm"))
	
	plotAnnotation <- grid.arrange(plotMajor, plotMinor, ncol = 1)
	ggsave(paste("STR_ALL_ANNOTATION_EWCE_CLUSTER", i, ".pdf", sep = ""), plotAnnotation, width = 12, height = 12, units = "in", dpi = 150)
	}

save(majorAnnotation, minorAnnotation, file = "STR_ALL_ANNOTATION_EWCE.RData")


## session info
sessionInfo()

R version 3.4.1 (2017-06-30)
Platform: x86_64-apple-darwin15.6.0 (64-bit)
Running under: macOS  10.14.3

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] grid      parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] gridExtra_2.3       readxl_1.2.0        limma_3.38.3        cowplot_0.9.4       EWCE_1.2.0         
 [6] RColorBrewer_1.1-2  pheatmap_1.0.12     ggplot2_3.1.0       reshape2_1.4.3      gmp_0.5-13.2       
[11] multtest_2.38.0     Biobase_2.42.0      BiocGenerics_0.28.0

loaded via a namespace (and not attached):
 [1] progress_1.2.0       tidyselect_0.2.5     purrr_0.3.0          splines_3.5.1        lattice_0.20-38     
 [6] colorspace_1.4-0     stats4_3.5.1         blob_1.1.1           survival_2.43-3      XML_3.98-1.16       
[11] rlang_0.3.1          pillar_1.3.1         glue_1.3.0           withr_2.1.2          DBI_1.0.0           
[16] bit64_0.9-7          bindrcpp_0.2.2       bindr_0.1.1          plyr_1.8.4           stringr_1.3.1       
[21] cellranger_1.1.0     munsell_0.5.0        gtable_0.2.0         memoise_1.1.0        IRanges_2.16.0      
[26] biomaRt_2.38.0       AnnotationDbi_1.44.0 Rcpp_1.0.0           scales_1.0.0         S4Vectors_0.20.1    
[31] bit_1.1-14           hms_0.4.2            digest_0.6.18        stringi_1.2.4        dplyr_0.7.8         
[36] tools_3.5.1          bitops_1.0-6         magrittr_1.5         lazyeval_0.2.1       RCurl_1.95-4.11     
[41] tibble_2.0.1         RSQLite_2.1.1        crayon_1.3.4         pkgconfig_2.0.2      MASS_7.3-51.1       
[46] Matrix_1.2-15        prettyunits_1.0.2    httr_1.4.0           assertthat_0.2.0     rstudioapi_0.9.0    
[51] R6_2.3.0             compiler_3.5.1
```



#### REFERENCES

1. Saunders et al., 2018, Cell 174, 1015-1030. <https://doi.org/10.1016/j.cell.2018.07.028>
2. Skene et al., 2016, Front. Neurosci. 10:16. <https://doi.org/10.3389/fnins.2016.00016>



