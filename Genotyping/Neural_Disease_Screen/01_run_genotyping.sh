#!/bin/bash

## Cell Barcodes from 10X CellRanger 
cellBarcodePath=dm-ter-screen1-1_cell_barcodes.tsv

## gRNA Fastq File Prefix
fastqName=DM_sgRNA_amp_M1_1_S54_L005

## Output dictionary name
outputFileName=dm-ter-screen1-1_pheno_dict.csv

## File mapping gRNA guide sequences to genes
barcodeToGene=dm-ter_screen_library.txt

## Sequences directly before and after the guide sequence
handle1='AACACCG'
handle2='GTTTTAGAGCTA'

## Path to PicardTools
PicardToolsPath=$HOME/PicardTools/picard.jar

## Path to genotyping-matrices
GenotypingMatricesPath=$HOME/genotyping-matrices/

## Convert from fastq to bam	
java -Xmx8g -jar $PicardToolsPath FastqToSam F1=$fastqName\_R2_001.fastq.gz O=$fastqName\_R2.bam SM=$fastqName SORT_ORDER=unsorted TMP_DIR=tmp/

## Tag Read2 bam with Read1 cell/molecule barcodes
python $GenotypingMatricesPath/convert_gRNA.py $fastqName\_R1_001.fastq.gz $fastqName\_R2.bam 
samtools index $fastqName\_R2.tagged.bam

## Parse cell barcodes and store in python dictionary
python $GenotypingMatricesPath/parseCellBarcodes.py $cellBarcodePath $fastqName\_R2.tagged.bam 20 $handle1 $handle2 
## Plot distribution of reads per UMI in order to set appropriate thresholds
python $GenotypingMatricesPath/plotUMIdist.py $fastqName\_R2_cell_barcodes.pickle 0.01

## Genotype cells: outputs a genotype dictionary in csv format
python2 $GenotypingMatricesPath/genotypeCells.py $fastqName\_R2_cell_barcodes.pickle $outputFileName $barcodeToGene 0.1 0.1 0.001
