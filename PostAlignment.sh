# Collect alignment statistics chromosome-wise. This is useful to get an average estimate of mitochondrial reads in the data. 
# Typical estimates are cell specific, but can range from 5% to 70%
samtools idxstats TCCTGAGC.sort.bam | head -n 22 > file.sort.chromosomestats.txt

# Collect reads that map to chromosomes under study. Remove sex chromosomes and MT
samtools view -b TCCTGAGC.sort.bam 1 10 11 12 13 14 15 16 17 18 19 2 20 3 4 5 6 7 8 9 >  file.sort.filtered.bam
samtools index file.sort.filtered.bam

# Remove PCR Duplicates using PICARD. Needs high memory.
picard MarkDuplicates I=file.sort.filtered.bam  O=file.sort.filtered.rmdup.bam M=file.duplicatemetrics.txt REMOVE_DUPLICATES=TRUE CREATE_INDEX=true
samtools idxstats file.sort.filtered.rmdup.bam | head -n 22 > file.chromosomestats.txt

#Convert the BAM files to BED
module load bedtools
bedtools bamtobed -i file.sort.filtered.rmdup.bam > file.bed

# Read start site shifting. 
# Adjust the reads start site to represent the center of transposon binding events. 
# Tn5 transposase binds as a dimer and inserts 2 adapters separated by 9 bp.
# To fix this, shift the + strand, offset by +4bp and - strand, offset by -5bp (Adey et al 2010)
awk 'BEGIN {OFS = "\t"} ; {if ($6 == "+") print $1, $2 + 4, $3 + 4, $4, $5, $6; else print $1, $2 - 5, $3 - 5, $4, $5, $6}' input.bed > input_shift.bed

#This BED file should be used for all downstream analyses.
