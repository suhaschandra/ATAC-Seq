# Install zimba
R
install.packages(c("R.oo"))
install.packages(c("R.utils", "quantreg","doParallel","doMC","foreach"))  # for R>3.0
# only version that works with R 3.0:
# get it from here: https://code.google.com/p/zinba/issues/detail?id=69
install.packages("zinba_2.03.1.tar.gz", repos=NULL)

# Make bed files from bams
# system("bedtools bamtobed -i /data/mapped/sample.bam > ~/zinba/reads/sample.bed")
# system("bedtools bamtobed -i /data/mapped/control.bam > ~/zinba/reads/control.bed")

# Calling peaks with ZINBA
library(zinba)

root = "/zinba"

sampleName = "sample"
controlName = "control"

# dirs
reads = file.path(root, "reads", fsep = .Platform$file.sep)
mappability = file.path(root, "mappability", fsep = .Platform$file.sep)
alignability = file.path(root, "alignability", fsep = .Platform$file.se)
alignCount = file.path(root, "aligncount", fsep = .Platform$file.sep)
output = file.path(root, "output", sampleName, fsep = .Platform$file.sep)

dir.create(reads)
dir.create(mappability)
dir.create(alignability)
dir.create(alignCount)
 
# static files & variables
genome = file.path(root, "genome", "hg19.2bit", fsep = .Platform$file.sep)

treatment = file.path(reads, paste0(sampleName, ".bed"), fsep = .Platform$file.sep)
control = file.path(reads, paste0(controlName, ".bed"), fsep = .Platform$file.sep)

alignCountOutput = file.path(alignCount, paste0(sampleName, ".basecount"), fsep = .Platform$file.sep)
averageFragmentLength=80

# Prepare for peak calling
generateAlignability(
    mapdir=mappability,  #mappability directory from unpacked mappability files
    outdir=alignability,  #directory for processed files, used later in analysis
    athresh=1,  #number of hits per read allowed during mapping process
    extension=50,  #average fragment library length
    twoBitFile=genome  #path to downloaded genome build file in .2bit format
)

basealigncount(  # only required for peak refinement, can be skipped otherwise
    inputfile=treatment,  #mapped sample reads
    outputfile=alignCountOutput,  # output path
    extension=averageFragmentLength,  #average fragment library length
    filetype="bed",  #either "bed", "bowtie", or "tagAlign"
    twoBitFile=genome  #path to downloaded genome build file in .2bit format
)

# Run the zinba pipeline
zinba(
    refinepeaks=1,  #refine peaks? 1 for yes, 0 for no
    seq=treatment,  #path to mapped experimental reads
    input=control,  #path to mapped input reads if available (default is "none")
    filetype="bed",  #either 'bed', 'bowtie', or 'tagAlign'
    threshold=0.05,  #FDR threshold, default is 0.05
    align=alignability,  #path to alignability directory
    numProc=32,  #number of CPUs to use, must be less than max available   (default 1)
    twoBit=genome,  #path to genome build in .2bit format
    outfile=output,  #prefix for outputted files
    extension=averageFragmentLength,  #average fragment library length (size selected)
    # OPTIONAL PARAMETERS #
    basecountfile=alignCountOutput,  #path to basecount file if refinepeaks is 1
    broad=FALSE,  #broad setting, TRUE or FALSE (default)
    printFullOut=1,  #print original data with enrichment estimates, 1 for yes (more space required), 0 for no (default)
    interaction=TRUE,  #whether or not to considering interaction during model selection, TRUE (default) or FALSE
    mode="peaks",  #either "peaks" for peak calling (default) or "CNV" for calling likely amplified CNV regions for reads in "seq" (input reads are best)
    FDR=TRUE  #either TRUE (default) or FALSE. If false, then uses posterior probability to threshold peaks using 1-threshold
)
