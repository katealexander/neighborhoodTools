#!/usr/bin/env python3

## Calculate statistics for changes in gene expression of gene neighborhoods

import os
import sys
import pybedtools
from pybedtools import BedTool, Interval
import numpy as np
import pandas as pd
import random

#python3 identifyChangingNeighborhoods.v2.py hg38.chrom.sizes hg38.ncbiRefSeq.gtf deseq2File.txt binSizes
#for binSizes, comma-separate for multiple

## output file
outputFileName = 'neighborhoodResults.txt'

## number of permutations for permutation test
nperm = 1000

## p value cutoff for calling significant bins
cutoff = 0.001

## define bin sizes
binSizes = [int(sys.argv[4])]

## define experiment name
experimentName = sys.argv[3].split("/")[-1].split(".tsv")[0]
print("Working on", experimentName)

## chromosomes to include
chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"] # human chromosomes

## defining which columns of expression file have the info
# gene expression header
with open(sys.argv[3], "r") as f:
    header = f.readline().strip().split("\t")

fcName = "log2FoldChange"
fcIndex = header.index(fcName) #fold change index

padjName = "padj"
padjIndex = header.index(padjName) #adjusted p value index

geneName = "Symbol"
geneIndex = header.index(geneName) #gene symbol index

## get chromosome sizes into dictionary
chrSizes = {}
with open(sys.argv[1], "r") as f:
    for line in f:
        line = line.strip()
        chr = line.split("\t")[0]
        size = int(line.split("\t")[1])
        chrSizes[chr] = size
   
## Make sliding bins to iterate through based on genome size and binSizes. Creates a list of lists
print("Making sliding genomic bins. Bin sizes:", binSizes)
bins = []
for size in binSizes:
    slide = int(int(size)/4) # slide by 25% of the binSize
    for chromosome in chromosomes:
        chrSize = chrSizes[chromosome]
        start = 0
        stop = int(size)
        while stop <= chrSize:
            bin = [chromosome, start, stop] # add bin to bins
            bins.append(bin)
            start = start + slide
            stop = stop + slide
            

## load exisiting dataframe if it has already been computed
precomputedFileName = "precomputedDataframes/" + sys.argv[2][:-4] + "." + ".".join(str(i) for i in binSizes) + ".csv"
if os.path.exists(precomputedFileName):
    print("Loading precomputed dataframe...")
    genesInBinsDF = pd.read_csv(precomputedFileName)
    geneNames = genesInBinsDF.columns.tolist()
    nGenes = len(geneNames)
else:
    ## get genes into bed file format
    print("Getting genes into bed file format...")
    geneCoordsList = []
    geneNames = []
    with open(sys.argv[2], "r") as f:
        for line in f:
            line = line.strip()
            annotation = line.split("\t")[2]
            chr = line.split("\t")[0]
            geneName = line.split("\t")[8].split('"')[1]
            if annotation == "transcript" and chr in chromosomes and geneName not in geneNames:
                start = int(line.split("\t")[3])
                stop = int(line.split("\t")[4])
                score = line.split("\t")[5]
                strand = line.split("\t")[6]
                geneCoordsList.append([chr, start, stop, geneName, score, strand])
                geneNames.append(geneName)
    geneCoordsBed = BedTool(geneCoordsList)
            
    ## make bin by gene dataframe. Rows are bin indexes, columns are genes. Gene gets a 1 if it's in the bin
    print("Making dataframe for which genes belong in which genomic bin...")
    nGenes = len(geneNames)
    nBins = len(bins)
    zeroArray = np.zeros((nBins, nGenes))
    genesInBinsDF = pd.DataFrame(zeroArray, columns=geneNames)
    i = 0
    for bin in bins:
        binInterval = Interval(bin[0], bin[1], bin[2])
        interList = geneCoordsBed.all_hits(binInterval)
        names = [interval.name for interval in interList]
        if names:
            genesInBinsDF.loc[i, names] = [1] * len(names)

        i += 1
    ## save dataframe
    genesInBinsDF.to_csv(precomputedFileName, index=False)



## empty dataframe with all genes as rows and real + # permutations as columns
zeroArray = np.zeros((nGenes, nperm + 1))
upDownDF = pd.DataFrame(zeroArray, index=geneNames)

## set real up and down genes
print("Reading gene expression file to get increasing and decreasing genes...")
with open(sys.argv[3], "r") as f:
    header = f.readline().strip()
    for line in f:
        line.strip()
        fc = line.split("\t")[fcIndex]
        padj = line.split("\t")[padjIndex]
        gene = line.split("\t")[geneIndex]
        if gene in geneNames and padj != "NA" and fc != "NA":
            padj = float(padj)
            fc = float(fc)
            if padj < 0.05: #significant
                if fc < 0 : #down
                    upDownDF.loc[gene, 0] = -1
                else: #up
                    upDownDF.loc[gene, 0] = 1

## replace remaining columns with bootstrapped values from column 0 of upDownDF
print("Generating dataframe with permutations that match the experimental number of increasing and decreasing genes. Number of permutations: ", nperm)
realUpDown = upDownDF[0].values
nUp = upDownDF[0].value_counts()[1]
nDown = upDownDF[0].value_counts()[-1]
for i in range(1,nperm+1):
    randomUpIndexes = random.sample(range(0, len(upDownDF)), nUp)
    remainingIndexes = list(set(range(0, len(upDownDF))) - set(randomUpIndexes))
    randomDownIndexes = random.sample(remainingIndexes, nDown)
    upDownDF.iloc[randomUpIndexes, i] = 1
    upDownDF.iloc[randomDownIndexes, i] = -1
    upDownDF.iloc[randomUpIndexes, i] = 1
    upDownDF.iloc[randomDownIndexes, i] = -1

## Matrix multiplication.
## Each row is a bin, the first column is the sum of up/down genes for the real data, remaining is for permutations
print("Calculating the sum of up and down genes per genomic bin...")
genesPerBinDF = genesInBinsDF @ upDownDF

## Get the increasing or decreasing bins
upIndexes = []
downIndexes = []
for i, row in genesPerBinDF.iterrows():
    nAbove = sum(1 for value in row[1:] if value >= row[0])
    nBelow = sum(1 for value in row[1:] if value <= row[0])
    pAbove = (nAbove + 1) / (nperm + 1)
    pBelow = (nBelow + 1) / (nperm + 1)
    if pAbove < cutoff:
        upIndexes.append(i)
    elif pBelow < cutoff:
        downIndexes.append(i)
upBins = [bins[i] for i in upIndexes]
downBins = [bins[i] for i in downIndexes]

## sort, merge, and save significant bins
sigUpBinsBed = pybedtools.BedTool(upBins)
sigDownBinsBed = pybedtools.BedTool(downBins)
upSorted = sigUpBinsBed.sort()
downSorted = sigDownBinsBed.sort()
upMerged = upSorted.merge()
downMerged = downSorted.merge()

upBedName = 'beds/' + experimentName + "_" + ".".join(str(i) for i in binSizes) + '_neighborhoods_up.bed'
downBedName = 'beds/' + experimentName + "_" + ".".join(str(i) for i in binSizes) +'_neighborhoods_down.bed'

upMerged.saveas(upBedName)
downMerged.saveas(downBedName)

## calculations for output file
print("Generating outputs...")
variances = genesPerBinDF.var() # variance of as a measure of how neighborhoody our dataset is
nAbove = sum(1 for value in variances[1:] if value >= variances[0]) # count how many of permutation bins have higher variance than actual
pVariance = (nAbove + 1)/(nperm + 1)
genomeSize = sum(chrSizes.values())
upSize = sum(interval.length for interval in upMerged)
downSize = sum(interval.length for interval in downMerged)

## write a header if an output is not already present
outFileHeader = 'experimentName\tVariance\tpVariance\tnumUpGenes\tnumDownGenes\tpercentGenomeUpNeighborhood\tpercentGenomeDnNeighborhood\tnumPermutations\tbinSizes\n'
if not os.path.exists(outputFileName):
    with open(outputFileName, 'w') as f:
        f.write(outFileHeader)

## append to existing file
with open(outputFileName, 'a') as f:
    toWrite = "\t".join([
        experimentName,
        f"{variances[0]:.3f}",
        f"{pVariance:.3e}",
        str(nUp),
        str(nDown),
        f"{(upSize / genomeSize) * 100:.3f}",
        f"{(downSize / genomeSize) * 100:.3f}",
        str(nperm),
        ",".join(map(str, binSizes)) if isinstance(binSizes, (list, tuple)) else str(binSizes)
    ]) + "\n"
    f.write(toWrite)
