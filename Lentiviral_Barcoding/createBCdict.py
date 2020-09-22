# This is a wrapper script for genotyping cells
import sys
import pandas as pd
from spacerCalling import *
import numpy as np

barcodeFile = sys.argv[1]
outputFile = sys.argv[2]

readThreshFrac = float(sys.argv[3])
umiThreshFrac = float(sys.argv[4])

minFrac = float(sys.argv[5])

with open(barcodeFile) as f:
    cellBarcodes = cp.load(f)

all_bcs = {}
for cell in cellBarcodes.values():
    molTags = cell.molTags
    for bc in molTags:
        all_bcs[bc] = bc

genotype_dict = callBarcodes(cellBarcodes, all_bcs, readThreshFrac,
                             umiThreshFrac, minUMIFrac = minFrac,
                             minGBCreads = 1, minGBCumis = 1)

with open(outputFile, 'w') as f:
    for genotype,cells in genotype_dict.items():
        outLine = genotype + ',\"' + ",".join(cells) + '\"' + '\n'
        f.write(outLine)

noPlasmids = 0
gbc_dist = Counter()
for cell in cellBarcodes.values():
    if cell.type == 'usable':
        gbc_dist[len(cell.genotype)] += 1
    elif cell.type == 'noGRNA':
        gbc_dist[0] += 1
    else:
        noPlasmids += 1

print 'Total_Cells\t' + str(len(cellBarcodes))
for k,v in gbc_dist.items():
    print str(k) + "\t" + str(v)

