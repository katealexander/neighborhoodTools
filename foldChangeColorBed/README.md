# Rationale
Create a BED file from Deseq2 output file of RNA-seq data to visualize increasing (yellow) and decreasing (blue) genes on genome browser track

# Design
From a Deseq2 output file, generates a BED file with a color column and header that can be directly copy-pasted onto UCSC Genome Browser's "Add custom track" feature.

# Outputs
A bed file for each dataset in the "datasets" directory will be output into the "beds" directory

# Software
RStudio with: dplyr and scales

# Needed inputs
- A gtf file with the locations of genes (e.g. [hg38.ncbiRefSeq.gtf](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/))
- One or more deseq2 output files, stored in directory "datasets/"

# Running the script
Open in RStudio, update first line of code to your working directory. Then run.

# Considerations
We designed the script to accomodate different names for the Deseq2 columns. If the inputed file does not have one of the names, and error message will occur. Additional names can be added in this portion of the script:

    ```## alternative names for columns in expression table; they will be replaced with common name
    namesForAdjP <- c("padj", "p.adj", "adj.P.Val")
    namesForPvalue <- c("pvalue", "P.Value")
    namesForLog2FC <- c("log2FoldChange", "logFC")
    namesForSymbol <- c("Symbol", "geneID")```

# R session info
> sessionInfo()
R version 4.4.1 (2024-06-14)
Platform: aarch64-apple-darwin20
Running under: macOS 15.6

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/New_York
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] scales_1.4.0 dplyr_1.1.4 

loaded via a namespace (and not attached):
 [1] RColorBrewer_1.1-3 R6_2.6.1           tidyselect_1.2.1   farver_2.1.2       magrittr_2.0.4     glue_1.8.0         tibble_3.3.0       pkgconfig_2.0.3   
 [9] generics_0.1.4     lifecycle_1.0.4    cli_3.6.5          vctrs_0.6.5        withr_3.0.2        compiler_4.4.1     rstudioapi_0.17.1  tools_4.4.1       
[17] pillar_1.11.1      rlang_1.1.6    
