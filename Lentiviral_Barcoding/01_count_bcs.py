import pysam as ps
#import spacerCalling as sc
from hammingSearch import approxHammingSearch
from dedupBarcodes import exactExists,argmin,collapseBarcodesExact
from collections import Counter
from itertools import islice
import sys


BC_START_HANDLE = 'GGCTGTTACGCG'
BC_END_HANDLE = 'CTACTGAC'


def getBarcode(read, barcode_length, bc_start_handle, bc_end_handle,
               rev_comp = False):
    bc_start_handle = bc_start_handle.upper()
    bc_end_handle = bc_end_handle.upper()

    read_seq = read.sequence
    if rev_comp:
        read_seq = reverse_complement(read_seq)

    # ensure that the sequence right before the barcode is what we expect
    bcStart = list(approxHammingSearch(bc_start_handle, read_seq))
    if len(bcStart) < 1: 
        return False
    left_pointer = argmin(bcStart) + len(bc_start_handle)
    
    bcEnd = list(approxHammingSearch(bc_end_handle, read_seq))
    if len(bcEnd) < 1: 
        return False
    right_pointer = argmin(bcEnd)

    barcode = read_seq[left_pointer:right_pointer]

    # Ensure the read covers the entire barcode
    if (len(barcode) < barcode_length - 2) or (len(barcode) > barcode_length + 2):
        return False

    return barcode



complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(seq):    
    bases = list(seq) 
    bases = reversed([complement.get(base,base) for base in bases])
    bases = ''.join(bases)
    return bases


def main():
    fastq_file = sys.argv[1]
    min_reads = int(sys.argv[2])

    guideCounts = Counter()
    with ps.FastxFile(fastq_file) as fh:
        for read in fh:
            gbc = getBarcode(read, barcode_length = 20, bc_start_handle = BC_START_HANDLE,
                             bc_end_handle = BC_END_HANDLE, rev_comp = True)
            if gbc:
                guideCounts[gbc] += 1
     
    guideCounts = {k:v for k,v in guideCounts.items() if v >= min_reads}
    
    for k,v in guideCounts.items():
        print k + "\t" + str(v)



if __name__ == "__main__": main()
