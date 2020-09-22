import os
import sys
import subprocess as sp

amplicons_file = sys.argv[1]

with open(amplicons_file) as f:
    for line in f:
        fields = str.strip(line).split("\t")
        gRNA = fields[0]

        gRNA_r1 = gRNA + "_R1.fastq.gz"
        gRNA_r2 = gRNA + "_R2.fastq.gz"

        out_dir = gRNA + "_editing_rates/"
        if os.path.exists(out_dir):
            if os.path.exists(gRNA_r1) and os.path.exists(gRNA_r2):
                in_dir = "CRISPResso_on_" + gRNA + "_R1_" + gRNA + "_R2/"
            else:
                in_dir = "CRISPResso_on_" + gRNA + "_R1/"

            quant_file = out_dir + in_dir + "Quantification_of_editing_frequency.txt"
            
            if os.path.exists(quant_file):
                lines = []
                with open(quant_file) as f:
                    for line in f:
                        lines.append(line)
                unmodified = lines[1].split(":")[1].split(" ")[0]
                total_reads = lines[6].split(":")[1].split(" ")[0]
                rate = 1 - float(unmodified)/float(total_reads)
                print gRNA + "\t" + str(rate) + "\t" + total_reads
 
