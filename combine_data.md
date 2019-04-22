# COMBINE DATA

This early post-natal striatal single cell RNA-seq data was generated for four genotypes (i) CTL, (ii) D1cKO, (iii) D2cKO & (iv) DDcKO. Each genotype has 4 biological replicates (4 mice). Single cell RNA-seq library was generated for each sample using 10X Genomics v2 microfluidics platform. And each library was sequenced for more than once (technical replicates) to increase the reads depth. Purpose of this script is to collapse the technical replicates for each library for each genotype followed by combining all libraries from all genotypes into a single UMI counts matrix for downstream processing. The code chunks need to be re-run changing the library/genotype name in order to process all libraries for all genotypes.



#### COLLAPSE TECHINCAL REPLICATES PER LIBRARY/SAMPLE

###### Read Reference Genes ID/Names List (Gencode vM17)

Raw UMI counts matrix for each sequencing run of every library is not identical due to many reasons. The goal of this part is to match all technical replicates using a reference list of gene names/symbols from Gencode vM17 reference annotation. 

```{r}
refGenes <- read.table("gencode.vM17.protein_coding_gene_id_symbol.txt", sep = "\t", header = T)
refGenesSymbol <- as.data.frame(unique(sort(refGenes$GeneSymbol)))
colnames(refGenesSymbol) <- "genes"
row.names(refGenesSymbol) <- refGenesSymbol$genes
```



###### Read and update UMI count tables

Read UMI Counts table for each technical run for each library/sample and merge the table with reference gene symbols. For the genes missing in the counts table, replace NAs add zeros. The data frames for each technical replicate are organized in a list *dataList* and list of cell barcodes from each technical replicate are organized in a list *libCells*.

```{r}
## reading the list of count tables for a specific library
## and merging with reference genes list from reference annotation
nameLIB <- "AA1" # library name
filesList <- list.files(path = ".", pattern = "*AA1_Counts.tsv.gz")

dataList <- vector("list", length(filesList))
names(dataList) <- gsub("_Counts.tsv.gz", "", filesList)

libCells <- vector("list", length(filesList))
names(libCells) <- gsub("_Counts.tsv.gz", "", filesList)

for(i in 1:length(dataList))
 {
 techRun <- read.table(paste(names(dataList)[[i]], "_Counts.tsv.gz", sep = ""), header = T, sep = "\t", row.names = 1)
 print(paste("Technical Run", i, "Data:", ncol(techRun), "Cells", sep = " "))

 techRun2 <- merge(refGenesSymbol, techRun, by = "row.names", all.x = T)
 row.names(techRun2) <- techRun2$Row.names
 techRun2$Row.names <- NULL
 techRun2$genes <- NULL

 techRun2[is.na(techRun2)] <- 0

 dataList[[i]] <- techRun2
 libCells[[i]] <- colnames(techRun2)
 }
```



###### Identify and extract common cells for each library

Using the list of lists for cell barcodes (*libCells*) and list of data frames (*dataList*) for each library created in earlier step, identify cell barcodes that are common to all technical replicates and extract data from each technical replicate corresponding to common cell barcodes into another list of data frames *dataCOMMON*. To generate a plot of common cell barcodes across technical replicates, script uses *UpSetR* R package.

```{r}
## identify common cell barcodes & generate a corresponding plot using UpSetR package
pdf(paste("AA_10X", nameLIB, "Common_CellBarcodes.pdf", sep = "_"), width = 6, height = 4)
libPlot <- upset(fromList(libCells), order.by = "freq")
dev.off()

libCombos <- Reduce(c, lapply(2:length(libCells), function(x) combn(1:length(libCells), x, simplify=FALSE) ))
libIntersect <- lapply(libCombos, function(x) Reduce(intersect, libCells[x]) )

commonCellsLib <- libIntersect[[max(length(libIntersect))]]
print(paste("Common Cells:", length(commonCellsLib), sep = " "))

## write the list of common cell barcodes for each library
write.table(commonCellsLib, paste("AA_10X", nameLIB, "Common_CellBarcodes.txt", sep = "_"), row.names = F, col.names = F, quote = F, sep = "\t")

## fetch the UMI counts data for common cells
dataCOMMON <- vector("list", length(dataList))
names(dataCOMMON) <- names(dataList)

for(j in 1:length(dataCOMMON))
 {
 dataCOMMON[[j]] <- dataList[[j]][,colnames(dataList[[j]]) %in% commonCellsLib]
 colnames(dataCOMMON[[j]]) <- paste(names(dataCOMMON)[[j]], colnames(dataCOMMON[[j]]), sep = "_")
 dataCOMMON[[j]]$genes <- row.names(dataCOMMON[[j]])
 }
```



###### Collapse common cells across technical replicate

UMI counts are then collapsed for cell barcodes that are common to all technical replicates for each library. Collapsed UMI count data for all technical replicates per library are stored in RData for further use.

```{r}
## collapse the UMI counts for common cells by adding the counts
combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , dataCOMMON)
row.names(combinedData) <- combinedData$genes
combinedData$genes <- NULL

combinedMeta <- as.data.frame(matrix(unlist(strsplit(colnames(combinedData), "_")), ncol = 4, byrow = TRUE))
row.names(combinedMeta) <- colnames(combinedData)
colnames(combinedMeta) <- c("Genotype", "Run", "Library", "CellBarcode")

newDataTemp <- as.data.frame(t(combinedData))
newData <- merge(newDataTemp, combinedMeta, by = "row.names")
row.names(newData) <- newData$Row.names
newData$Row.names <- NULL
cellBC <- list(newData$CellBarcode)
newData$Genotype <- NULL
newData$Library <- NULL
newData$CellBarcode <- NULL
newData$Run <- NULL

newDataAggrTemp <- aggregate(newData, by = cellBC, FUN = "sum")
row.names(newDataAggrTemp) <- newDataAggrTemp$Group.1
newDataAggrTemp$Group.1 <- NULL

newDataAggr <- as.data.frame(t(newDataAggrTemp))
colnames(newDataAggr) <- paste(nameLIB, colnames(newDataAggr), sep = "_")
newDataAggr$genes <- row.names(newDataAggr)

## save collapsed data frame for each library
save(newDataAggr, file = paste(nameLIB, "Collapsed_Data.RData", sep = "_"))
```



#### COMBINE LIBRARIES/SAMPLES PER GENOTYPE

In the next step, collapsed UMI counts per library/sample for each genotype is combined into a single data frame. Combined genotype data is stored in RData for further use.

```{r}
## combine all libraries for one genotype
listRep <- list("AA1", "AA4", "AA7", "AA8")
nameCat <- "CTL"
print(paste("Processing", length(listRep), "replicates for", nameCat, sep = " "))
names(listRep) <- listRep

dataCOMBINED <- vector("list", length(listRep))
names(dataCOMBINED) <- names(listRep)

for(i in 1:length(dataCOMBINED))
 {
 print(paste("Reading data for ", names(dataCOMBINED)[i], sep = ""))
 dtlist <- load(paste(listRep[[i]], "_Collapsed_Data.RData", sep = ""))
 dataCOMBINED[[i]] <- newDataAggr
 rm(list = dtlist)
 rm(dtlist)
 }

combinedData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , dataCOMBINED)
row.names(combinedData) <- combinedData$genes
combinedData$genes <- NULL
colnames(combinedData) <- paste(nameCat, colnames(combinedData), sep = "_")
print(dim(combinedData))

## save combined data frame
save(combinedData, file = paste(nameCat, "_Collapsed_Combined.RData", sep = ""))
```



#### COMBINE LIBRARIES/SAMPLES FROM MULTIPLE GENOTYPES

As a last step, UMI count tables corresponding to each genotype are combined to create a massive data frame will all libraries/samples. This combined data frame is stored in RData and becomes the starting point for data quality check and filtering followed by downstream clustering and differential gene expression analysis.

```{r}
## combine all genotypes together
listGeno <- list("CTL", "D1CKO", "D2CKO", "DDCKO")
nameCat <- "STR"
names(listGeno) <- listGeno

dataCOMBINED <- vector("list", length(listGeno))
names(dataCOMBINED) <- names(listGeno)

for(i in 1:length(dataCOMBINED))
 {
 print(paste("Reading data for ", names(dataCOMBINED)[i], sep = ""))
 dtlist <- load(paste(listGeno[[i]], "_Collapsed_Combined.RData", sep = ""))
 combinedData$genes <- row.names(combinedData)
 dataCOMBINED[[i]] <- combinedData
 rm(list = dtlist)
 rm(dtlist)
 }

allData <- Reduce(function(x, y) { merge(x, y, all = TRUE, by = "genes") } , dataCOMBINED)
row.names(allData) <- allData$genes
allData$genes <- NULL
colnames(allData) <- paste(nameCat, colnames(allData), sep = "_")

## save combined data frame for all genotypes
save(allData, file = paste(nameCat, "_Collapsed_Combined_Genotypes.RData", sep = ""))
```



The final combined matrix has a total of 63,984 cells.

| Genotype  | Cells      |
| --------- | ---------- |
| CTL       | 14,603     |
| D1cKO     | 17,206     |
| D2cKO     | 10,001     |
| DDcKO     | 22,174     |
| **Total** | **63,984** |



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
[1] UpSetR_1.3.3

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.19     crayon_1.3.4     dplyr_0.7.7      assertthat_0.2.0
 [5] grid_3.4.1       plyr_1.8.4       R6_2.3.0         gtable_0.2.0    
 [9] magrittr_1.5     scales_0.5.0     ggplot2_3.1.0    pillar_1.3.0    
[13] rlang_0.3.0.1    lazyeval_0.2.1   bindrcpp_0.2.2   labeling_0.3    
[17] glue_1.3.0       purrr_0.2.5      munsell_0.5.0    compiler_3.4.1  
[21] pkgconfig_2.0.2  colorspace_1.3-2 tidyselect_0.2.5 bindr_0.1.1     
[25] gridExtra_2.3    tibble_1.4.2  
```




