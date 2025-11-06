setwd("~/Desktop/foldChangeColorBed/")

library(dplyr)
library(scales)

## make a directory if one does not already exist
if (!dir.exists("beds")) dir.create("beds")

genes <- read.delim(file = "hg38.ncbiRefSeq.gtf", header = F, sep = "\t")
genes <- genes[genes$V3 == "transcript",]
genes$gene <- sub(".*gene_name ([^;]+);.*", "\\1", genes$V9) 
colnames(genes) <- c("chr", "dataset", "type", "start", "stop","anno1","strand", "anno2", "info", "gene")
genes$TSS <- ifelse(genes$strand == "+", genes$start, genes$stop)
chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")
genes <- genes[genes$chr %in% chromosomes, ]

datasets <- list.files(path = "datasets", full.names = TRUE) #this will be a list of our datasets

## alternative names for columns in expression table; they will be replaced with common name
namesForAdjP <- c("padj", "p.adj", "adj.P.Val")
namesForPvalue <- c("pvalue", "P.Value")
namesForLog2FC <- c("log2FoldChange", "logFC")
namesForSymbol <- c("Symbol", "geneID")

# Make a color ramp from blue -> white -> red
color_fun <- colorRampPalette(c("#2166ac", "#f7f7f7", "#E1AD01"))

for(data in datasets) {
  ## Load the gene expression data
  #data <- "datasets/GSE272396.THP1IfnTreated.12h.tsv"
  dataTable <- read.delim(file = data, sep = "\t", header = TRUE) #loading the dataset
  dataTable <- dataTable[,colSums(is.na(dataTable)) != nrow(dataTable)]
  dataTable <- na.omit(dataTable)
  nameOfExperiment <- sub(pattern = ".*datasets/(.*)\\.tsv$", replacement = "\\1", x = data)
  
  ###### if the column names are not as expected, replace them with common name in expression table ###### 
  # if the common name is not in the expression table colnames, but an alternative name is, replace the alternative name with the common name
  # otherwise create an error message that it could not be found
  # check padj
  if (!"adj.P.Val" %in% colnames(dataTable)){
    if (any(namesForAdjP %in% colnames(dataTable))) {
      # find the position
      positionIndex <- which(colnames(dataTable) %in% namesForAdjP)
      # rename the column at that position to the common name
      names(dataTable)[positionIndex] <- "adj.P.Val"
    } else {
      errorMessage <- paste0("Could not find the adjusted p-value column name in: ", nameOfExperiment)
      print(errorMessage)
      next
    }
  }
  
  # check p-value
  if (!"P.Value" %in% colnames(dataTable)){
    if (any(namesForPvalue %in% colnames(dataTable))) {
      # find the position
      positionIndex <- which(colnames(dataTable) %in% namesForPvalue)
      # rename the column at that position to the common name
      names(dataTable)[positionIndex] <- "P.Value"
    } else {
      errorMessage <- paste0("Could not find the p-value column name in: ", nameOfExperiment)
      print(errorMessage)
      next
    }
  }
  
  # check fold change
  if (!"logFC" %in% colnames(dataTable)){
    if (any(namesForLog2FC %in% colnames(dataTable))) {
      # find the position
      positionIndex <- which(colnames(dataTable) %in% namesForLog2FC)
      # rename the column at that position to the common name
      names(dataTable)[positionIndex] <- "logFC"
    } else {
      errorMessage <- paste0("Could not find the logFC column name in: ", nameOfExperiment)
      print(errorMessage)
      next
    }
  }
  
  # check symbol
  if (!"Gene.symbol" %in% colnames(dataTable)){
    if (any(namesForSymbol %in% colnames(dataTable))) {
      # find the position
      positionIndex <- which(colnames(dataTable) %in% namesForSymbol)
      # rename the column at that position to the common name
      names(dataTable)[positionIndex] <- "Gene.symbol"
    } else {
      errorMessage <- paste0("Could not find the gene symbol column name in: ", nameOfExperiment)
      print(errorMessage)
      next
    }
  }
  
  dataTable <- dataTable[dataTable$adj.P.Val < 0.05,]
  df <- merge(genes, dataTable, by.x = "gene", by.y = "Gene.symbol", all = F)
  df <- df[!duplicated(df$gene), ]
  
  # Map fold changes to the color ramp:
  logFC <- df$logFC
  # Clip logFC to a range
  logFC <- pmax(pmin(logFC, 0.6), -0.6)
  # Scale to 0–1
  scaled <- rescale(logFC, from=c(-0.6,0.6))
  
  # Get colors as hex
  cols <- color_fun(100)  # 100-color palette
  hex_colors <- cols[as.integer(scaled*99)+1]
  
  # Convert hex to RGB strings for UCSC
  rgb_mat <- col2rgb(hex_colors)
  rgb_str <- apply(rgb_mat, 2, function(x) paste(x, collapse=","))
  
  # Build BED
  bed <- df %>%
    mutate(start = TSS,
           end = TSS + 50,
           cleanExp = gsub("\\s+", "_", nameOfExperiment),
           cleanSym = gsub("\\s+", "_", gene),
           name = paste0(cleanExp, "_", cleanSym),
           score = 0,
           strand = ".",
           x = start,
           y = end,
           color = rgb_str) %>%
    select(chr, start, end, name, score, strand, x, y, color)
  
  # Write header + data
  safe_experiment <- gsub("\\s+|/|\\|", "_", nameOfExperiment)  # replace spaces, slashes, pipes
  fileName <- paste0("beds/", safe_experiment, "_foldchange_gradient.bed")
  cat(
    paste0(
      'track name="', safe_experiment, '_FoldChangeGradient" ',
      'description="Fold change gradient for ', nameOfExperiment, '" ',
      'itemRgb="On"\n'
    ),
    file = fileName
  )
  write.table(bed, fileName,
              quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE)
  
}



