## Compute fraction of barcodes post vs pre injection

samples = ["dm-ter-bc1", "dm-ter-bc2", "dm-ter-bc3"]
for sample in samples:
    sample_pre_file = sample + "_pre_gDNA_bc_counts.txt"
    sample_post_file = sample + "_post_gDNA_bc_counts.txt"
    
    sample_pre_barcodes = set()
    with open(sample_pre_file) as f:
        for line in f:
            bc = line.split("\t")[0]
            sample_pre_barcodes.add(bc)
    
    num_post = 0
    num_novel = 0
    with open(sample_post_file) as f:
        for line in f:
            bc = line.split("\t")[0]
            if bc in sample_pre_barcodes:
                num_post += 1
            else:
                n_reads = int(line.split("\t")[1])
                if n_reads > 10:
                    num_novel += 1

    num_pre = len(sample_pre_barcodes)
    barcode_frac = float(num_post)/num_pre
    
    print sample + "\t" + str(num_pre) + "\t" + str(num_post) + "\t" + str(barcode_frac) + "\t" + str(num_novel)
