# Nucleosome position analysis. For this analysis merge the samples by group i.e. have just 1 bam file per condition

# Reference file index
samtools idx Macaca_mulatta.MMUL_1.dna.toplevel.fa
# Merge bam files per condition, remember to use the file that is PCR duplicate, MT and X chromosome removed, and sorted and indexed
samtools merge control.bam ../../duplicate_removal/AGGCAGAA.sort.filtered.rmdup.bam ../../duplicate_removal/CGTACTAG.sort.filtered.rmdup.bam

# Generate index of the merged files
samtools index control.bam

#Call broad peaks using MACS2
macs2 callpeak --name Control -t control.bam --outdir Control --format BAM -nomodel --shift 100 --extsize 200 --broad

#Run NucleoATAC
nucleoatac run --bed ../mergebams/Control/Control_peaks.broadPeak.bed --bam ../mergebams/control.bam --fasta ../Reference/Macaca_mulatta.MMUL_1.dna.toplevel.fa --out Control
