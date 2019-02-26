# DEMULTIPLEXING
..* demultiplex the raw sequencing data in the form of BCL files
..* resulting files are FASTQ.GZ files
..* uses Illumina's bcl2fastq and 10X Genomics CellRanger software program
```shell
cellranger mkfastq --run=<path to BCL files> --samplesheet=<path to sample sheet>
```

# QUALITY CHECK
..* run quality check on the fastq.gz files
..* uses FASTQC software program
```shell
ls *.fastq.gz | sed "s/.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; fastqc %.fastq.gz'
```

## WHITELIST
## estimate and create a whitelist of real cell barcodes from the quality filtered fastq.gz files
## uses UMI Tools software program

ls *R1*.gz | sed "s/_R1.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; umi_tools whitelist --stdin=%_R1.fastq.gz --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --expect-cells=<number of expected cells> --plot-prefix=%_Expect_Whitelist --log=%_Whitelist_log.out --error=%_Whitelist_log.err --stdout=%_Whitelist.txt'

## EXTRACT READS
## extract the reads corresponding to estimated whitelist of real cell barcodes
## also appends the cell-barcode and umi information from R1 fastq files to read names of R2 fastq files
## uses UMI Tools software program

ls *R1*.gz | sed "s/_R1.fastq.gz//g" | xargs -I % -n 1 -P 48 sh -c 'echo %; umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCNNNNNNNNNN --stdin=%_R1.fastq.gz --stdout=%_R1_extracted.fastq.gz --read2-in=%_R2.fastq.gz --read2-out=%_R2_extracted.fastq.gz --filter-cell-barcode --whitelist=%_Whitelist.txt --log=%_Extract_log.out --error=%_Extract_log.err'

## READS ALIGNMENT
## align the extracted reads with reference genome and annotation
## uses STAR software program

## Prior alignment, reference genome index need to be built
## This is one-time only step and can be used multipli times for same genome build.
## Genome build: Mouse genome GRCm38.p6 (MM10)
## Annotation: Gencode vM17

STAR --runMode genomeGenerate --runThreadN 48 --genomeDir <genome index directory> --genomeFastaFiles <reference genome fasta file> --sjdbGTFfile <reference annotation gtf file> --sjdbOverhang 100


## Once the genome index is ready, run the alignment as below
for FQFILE in `ls *R2_extracted*.gz`
 do
  prefx=`echo ${FQFILE} | sed "s/_R2_extracted.fastq.gz//g"`
  echo "Processing" ${prefx}

  STAR --runThreadN 48 \
       --genomeDir <path the reference genome index> \
       --readFilesIn ${FQFILE} \
       --readFilesCommand zcat \
       --sjdbGTFfile <path to reference annotation gtf file> \
       --outFilterType BySJout  \
       --outFilterMismatchNoverReadLmax 0.04 \
       --outFilterMultimapNmax 10 \
       --alignSJoverhangMin 10 \
       --alignSJDBoverhangMin 1 \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMismatchNmax 5 \
       --twopassMode Basic \
       --outFileNamePrefix ${prefx}_STAR_

 done
