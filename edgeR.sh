# Code to identify differential peaks between two conditions using edgeR
# Step 1: Merge all the bam files from the analysis - this is used to identify regions of detectable chromatin regions

# Step 2: Call peaks on the above file using MACS2 with a loose threshold of 0.1
macs2 callpeak --name Maob_CD4_allsamples -t ../position_adjust/Maob_cd4_combined_shift.bed --outdir Maob_CD4_allsamples/ --format BED --nomodel --shift 100 --extsize 200 -B --SPMR -q 0.1

# Collect the summit values and extend 200 bp on either side to generate an saf file. These are constant length regions which will be used to count

# Use the above regions to perform individual read counting. Combine all counts into an excel file

# Generate the targets file and perform edgeR!
