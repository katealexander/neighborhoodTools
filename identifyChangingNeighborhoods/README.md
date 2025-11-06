# Rationale
Detect increasing and decreasing gene neighborhoods from Deseq2 output files of RNA-seq data.

# Design
This analysis detects genomic neighborhoods that contain significantly more increasing or decreasing genes in a pairwise comparison than expected by random chance. It calculates significance using a permutation test of the same number of randomly selected genes. After identifying significantly increasing or decreasing genomic bins, it uses pybedtools to merge overlapping bins and returns these as bed files. 

# Outputs
1. Bed files with significantly up an down gene neighborhoods.
2. A tab-delimited text file, "neighborhoodResults.txt" that provides summaries for all datasets analysed. 

# Software
Python 3 with: numpy, pybedtools, sys, os, pandas, and random

To get python package:

```pip install pybedtools```

# Needed inputs
- A gtf file with the locations of genes (e.g. [hg38.ncbiRefSeq.gtf](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz))
- A chromosome sizes file (e.g. "hg38.chrom.sizes")
- One or more deseq2 output files, stored in directory "datasets/"

# Before running
1. Navigate to the directory containing the script, "identifyChangingNeighborhoods.v2.py" in terminal
    ```cd path/to/directory```
2. Ensure the following directories are present: beds, and precomputedDataframes. If not, add the directories:
    ```mkdir beds```
    ```mkdir precomputedDataframes```

# Running the script
To run the script for one dataset:

```python3 identifyChangingNeighborhoods.v2.py hg38.chrom.sizes hg38.ncbiRefSeq.gtf datasets/DATASETNAME.tsv 2500000```

The last argument is the bin sizes. This can be just one bin size, as shown for 1Mb size neighborhoods. Or, if desired, it can be multiple, for example: 1000000,5000000

To run the script on all the datasets in the datasets directory:

```for dataset in datasets/*; do python3 identifyChangingNeighborhoods.v2.py hg38.chrom.sizes hg38.ncbiRefSeq.gtf $dataset 2500000; done```

# Considerations
- You will want to ensure the dataset of interest has enough differential genes to make the analysis interpretable. For example, neighborhoods will likely be more challenging to detect with 50 differential genes versus 2000.
- Smaller bin sizes will make things take longer. Bin size should be selected based on what is seen in your specific dataset. Load the bed files onto a genome browser and compare with visualizations of increasing and decreasing genes (see foldChangeColorBed)
- The script stores a dataframe based on the transcript locations used and bin sizes. This dataframe consists of 1s and 0s demarking which gene belongs in which bin. If you reuse the same binSizes and same transcript (gtf) file, the script will re-load this dataframe rather than regenerate it. If you edit any part of the script that impacts the bins variable, make sure to delete these files and re-generate them. 

# Software details
This is the exact version of python and all installed packages at the time of development.
Python 3.12.5
Package                Version
---------------------- -----------
cellpose               3.1.1.1
certifi                2025.4.26
charset-normalizer     3.4.2
docutils               0.21.2
fastremap              1.15.0
filelock               3.16.1
fsspec                 2024.12.0
idna                   3.10
imagecodecs            2024.12.30
Jinja2                 3.1.5
joblib                 1.5.2
llvmlite               0.43.0
MarkupSafe             3.0.2
mpmath                 1.3.0
natsort                8.4.0
ncls                   0.0.70
networkx               3.4.2
numba                  0.60.0
numpy                  1.26.4
opencv-python-headless 4.11.0.86
packaging              24.2
pandas                 2.3.3
pip                    25.3
pybedtools             0.12.0
Pygments               2.19.1
PyQt6                  6.8.0
PyQt6-Qt6              6.8.1
PyQt6_sip              13.9.1
pyqtgraph              0.13.7
pyranges               0.1.4
pysam                  0.23.3
python-dateutil        2.9.0.post0
pytz                   2025.2
QtPy                   2.4.2
requests               2.32.3
roifile                2024.9.15
scipy                  1.15.1
setuptools             75.8.0
six                    1.17.0
sorted_nearest         0.0.41
statistics             1.0.3.5
superqt                0.7.1
sympy                  1.13.1
tabulate               0.9.0
tifffile               2025.1.10
torch                  2.5.1
tqdm                   4.67.1
typing_extensions      4.12.2
tzdata                 2025.2
urllib3                2.4.0
