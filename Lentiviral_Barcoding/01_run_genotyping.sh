#!/bin/bash

## Path to cell barcodes output by CellRanger 
cellBarcodePath=dm-ter-bc1_barcodes.tsv

## Lentiviral barcode Fastq File Prefix
fastqName=DM_Barcoded-1_cDNA

## Output barcode dictionary name
outputFileName=dm-ter-bc1_pheno_dict.csv

## Path to PicardTools installation
PicardToolsPath=$HOME/PicardTools/picard.jar 

## Path to genotyping-matrices directory
GenotypingMatricesPath=$HOME/genotyping-matrices/

## Convert from fastq to bam	
java -Xmx8g -jar $PicardToolsPath FastqToSam F1=$fastqName\_R2.fastq.gz O=$fastqName\_R2.bam SM=$fastqName SORT_ORDER=unsorted TMP_DIR=tmp/

## Tag Read2 bam with Read1 cell/molecule barcodes
python $GenotypingMatricesPath/convert_gRNA.py $fastqName\_R1.fastq.gz $fastqName\_R2.bam 
samtools index $fastqName\_R2.tagged.bam

## Parse lentiviral barcodes and store in python dictionary
## Uses a slightly modified version of the parseCellBarcodes.py script in genotyping-matrices
python parseLentiBarcodes.py $cellBarcodePath $fastqName\_R2.tagged.bam 20 'GGCTGTTACGCG' 'CTACTGAC'
 
## Genotype cells: outputs a genotype dictionary in csv format
## Uses a modified version of the genotyping script in genotyping-matrices
python2 createBCdict.py $fastqName\_R2_cell_barcodes.pickle $outputFileName $barcodeToGene 0.05 0.05 0.0025
