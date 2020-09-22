# teratoma-analysis-code

## Processed Seurat objects and metadata



## Reproducing the core analysis starting from the processed counts matrices and genotype dictionaries

1. Clone this repository into your local directory: `git clone https://github.com/yanwu2014/teratoma-analysis-code.git`

2. Download the processed data from GEO: GSE156170
    1. GSE156170_teratoma_merged_counts.tar.gz
    2. GSE156170_external_reference_data.tar.gz

3. Untar the processed data files which should create a `Counts/` directory and a `Reference_Data` directory
````
tar xvf GSE156170_teratoma_merged_counts.tar.gz
tar xvf GSE156170_external_reference_data.tar.gz
````

4. Run the R scripts within each Figure directory to reproduce the analysis used to create that figure. The scripts need to be run in a certain order reflected in their numbering. For example run `01_human_clustering.R` before `02_human_cluster_mapping.R` in Figure1. Scripts with the same number can be run in any order. The scripts in Figure 1 need to be run first as they will generate the clustering results needed for the rest of the analysis.

## Computing the editing rates for the screens in Figure 4/S4 using the amplicon Fastq files

## Generating the genotype dictionaries in Figure 2 and Figure4/S4 from the lentiviral barcode/gRNA amplicon Fastq files

## Generating the counts matrices from the raw 10X cDNA fastq files
