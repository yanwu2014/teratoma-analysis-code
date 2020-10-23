# teratoma-analysis-code (under construction)

## Processed Seurat objects and metadata

The processed Seurat object for the H1 teratomas and the corresponding cell type annotations and metadata can be found on our [FTP server](shorturl.at/gwI08)


## Reproducing the core analysis starting from the processed counts matrices and genotype dictionaries

1. Clone this repository into your local directory: `git clone https://github.com/yanwu2014/teratoma-analysis-code.git`

2. Download the processed data from [GEO: GSE156170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156170)
    1. [GSE156170_teratoma_merged_counts.tar.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156170&format=file&file=GSE156170%5Fteratoma%5Fmerged%5Fcounts%2Etar%2Egz)
    2. [GSE156170_external_reference_data.tar.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156170&format=file&file=GSE156170%5Fexternal%5Freference%5Fdata%2Etar%2Egz)

3. Untar the processed data files which should create a `Counts/` directory and a `Reference_Data` directory
````
tar xvf GSE156170_teratoma_merged_counts.tar.gz
tar xvf GSE156170_external_reference_data.tar.gz
````

4. Run the R scripts within each Figure directory to reproduce the analysis used to create that figure. The scripts need to be run in a certain order reflected in their numbering. For example run `01_human_clustering.R` before `02_human_cluster_mapping.R` in Figure1. Scripts with the same number can be run in any order. The scripts in Figure 1 need to be run first as they will generate the clustering results needed for the rest of the analysis.


## Generating the genotype dictionaries in Figure 2 and Figure 4/S4 from the lentiviral barcode/gRNA amplicon Fastq files

The genotype dictionaries map either CRISPR-Cas9 gRNAs or lentiviral barcodes to single cells. We provide the processed genotype dictionaries in both this github repository and the supplementary files at [GEO: GSE156170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156170). The genotype dictionaries end in `pheno_dict.csv.gz`. 

### Reproducing the genotyping for Figure 4

To reproduce the genotype dictionaries for the embryonic lethal screen and replicate screen in Figure 4 we can download the original gRNA barcode Fastq files and reprocess them using either the files in the `Genotyping` directory. We'll walk through an example using the gRNA fastq files from one of the 10X runs from the embryonic lethal screen which you can download at [GSM4725940](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4725934):

1. Download and install [PicardTools](https://broadinstitute.github.io/picard/) if not already installed. These scripts assume the `picard.jar` file is at `$HOME/PicardTools/picard.jar` but you can install PicardTools anywhere as long as you edit the `PicardToolsPath` parameter in the genotyping scripts
2. Clone https://github.com/yanwu2014/genotyping-matrices into your home directory (or wherever you like as long as you edit the `GenotypingMatricesPath` in your genotyping scripts)
3. Download the gRNA fastq files from GEO and move to the `Genotyping/Lethal_Screen/` directory
4. Edit the fastqName parameter in the `01_run_genotyping.sh` script to match the files that were just downloaded. The full fastq file names will be assumed to be `[fastqName]_R1_001.fastq.gz` and `[fastqName]_R2_001.fastq.gz` in the rest of the shell script.
5. Run `01_run_genotyping.sh` which should generate a genotype dictionary
6. Download the rest of the gRNA fastq files (GSM4725934 - GSM4725939, GSM4725946 - GSM4725951) and change the cellBarcodePath and outputFileName to generate the remaining genotype dictionaries. For example if you want to generate the genotype dictionary for the embryonic lethal screen teratoma 2 10X replicate 1, set `cellBarcodePath=dm-ter-screen2-1_cell_barcodes.tsv` and `outputFileName=dm-ter-screen2-1_pheno_dict.csv`
7. To merge the genotype dictionaries you can use the `02_merge_pheno_dicts.R` script which has the usage: `Rscript 02_merge_pheno_dicts.R [output_merged_pheno_dict.csv] [input_pheno_dict_1.csv] [input_pheno_dict_2].csv ...`

### Reproducing the genotyping for Figure S4

To reproduce the genotype dictionaries for the neural disase screen in FigureS4 simply download the appropriate gRNA fastq files (GSM4725956 - GSM4725959) from [GEO: GSE156170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156170) and move them to the `Genotyping/Neural_Disease_Screen/` directory instead and use the same steps as for the embryonic lethal screen.

### Reproducing the single cell lentiviral barcoding for Figure 2

To reproduce the lentiviral barcoding genotype dictionaries again download the appropriate barcode fastq files (GSM4725916 - GSM4725918) from [GEO: GSE156170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156170) and move to the `Lentiviral_Barcoding/` directory

### Reproducing the pre and post injection barcode fraction in Figure 2

1. Download the pre and post teratoma injection gDNA fastq files (GSM4725919 - GSM4725924) from [GEO: GSE156170](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE156170) and move to the `Lentiviral_Barcoding/` directory.
2. Run `01_count_bcs.py` on each gDNA fastq file with the usage: `python 01_count_bcs.py [gDNA_fastq_file.fastq.gz] [min_barcode_reads] > [output_file]`. For the paper we set `min_barcode_reads = 1` but feel free to play around with the parameter to see how it affects the number of barcodes you get out. The `02_compute_barcode_frac.py` script by default expects barcode counts files in the format: `dm-ter-bc[1-3]_[pre/post]_gDNA_bc_counts.txt` so we recommend using that naming convention. 
3. If you used the naming convention recommended in the previous step you can simply run `python 02_compute_barcode_frac.py`. Otherwise edit `02_compute_barcode_frac.py` so that it lines up with your barcode counts file names.


## Computing the editing rates for the screens in Figure 4/S4 using the amplicon Fastq files

1. Install CRISPResso https://github.com/lucapinello/CRISPResso
2. Download the [GSE156170_editing_rate_data.tar.gz](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156170&format=file&file=GSE156170%5Fediting%5Frate%5Fdata%2Etar%2Egz) and extract.
3. Run `python 01_run_editing_analysis.py` which is a wrapper script that runs CRISPResso on all of the amplicon fastq files in the directory
4. Run `python 02_prase_editing_rates.py` which extracts the CRISPResso output and formats it into a single tab separated output


## Generating the counts matrices from the raw 10X cDNA fastq files

1. Download the merged hg19/mm10 genome reference from the CellRanger site: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
2. Download and install CellRanger https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
3. Download the fastq file for the 10X run you want to analyze (i.e. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4725909)
4. Run CellRanger using `cellranger count` on the fastq files using default settings (adjusting the local cores and local memory usage as appropriate for your system)
