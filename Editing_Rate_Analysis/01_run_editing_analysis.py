import sys
import subprocess as sp
import os

amplicons_file = sys.argv[1]

with open(amplicons_file) as f:
    for line in f:
        fields = str.strip(line).split("\t")
        gRNA = fields[0]
        amplicon = fields[1]
        gRNA_seq = fields[2]

        gRNA_r1 = gRNA + "_R1.fastq.gz"
        gRNA_r2 = gRNA + "_R2.fastq.gz"

        if os.path.isfile(gRNA_r1) and os.path.isfile(gRNA_r2):
            out_dir = gRNA + "_editing_rates"
            cmd = "CRISPResso -r1 " + gRNA_r1 + " -r2 " + gRNA_r2 + " -a " + amplicon + \
                    " -o " + out_dir + " -p 12 -g " + gRNA_seq + " --ignore_substitutions " + \
                    " --trim_sequences"
            sp.call(cmd, shell = True)
        
        elif os.path.isfile(gRNA_r1) and not os.path.isfile(gRNA_r2):
            out_dir = gRNA + "_editing_rates"
            cmd = "CRISPResso -r1 " + gRNA_r1 + " -a " + amplicon + \
                    " -o " + out_dir + " -p 12 -g " + gRNA_seq + " --ignore_substitutions " + \
                    " --trim_sequences"
            sp.call(cmd, shell = True)

        else:
            print gRNA + " fastq files not found"

            
