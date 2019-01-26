"""
Snakemake pipeline for ChIP-Seq peak calling in duplicated regions
"""

# -----------
# CONFIG FILE - reads in settings and parameters from config file
# -----------

REFERENCE = config["reference"]
CHROMSIZE = config["chromsize"]
SAMPLES = config["samples"]
CONTROL = config["control"]
FILENAME = config["filename"]
ADAPTERS = config["adapters"]
INDEXWINDOW = config["indexwindow"]
CROPSIZE = config["cropsize"]

if config["broad"]:
  PEAK_TYPE = "--broad"
else:
  PEAK_TYPE = ""

# ---------------
# WORKFLOW SET-UP - unsure
# ---------------

TARGETS = [expand("results/qc/{smp}_fastqc.html", smp = SAMPLES+CONTROL),
           expand("results/qc_trimmed/{smp}.trimmed_fastqc.html", smp = SAMPLES+CONTROL),
           "results/peak_calling/"+FILENAME+"_logLR.bw"]

if config["csem"]:
  TARGETS += ["results/peak_calling/"+FILENAME+"_csem_logLR.bw"]

# --------
# PIPELINE
# --------

rule all:
  input:
    TARGETS

# Takes in fastq data and generates a report on the quality of the sequencing reads
rule qc:
  input:
    "data/{smp}.fastq.gz"
  output:
    "results/qc/{smp}_fastqc.html"
  params:
    "results/qc"
  conda:
    "envs/fastqc.yaml"
  shell:
    """
    fastqc -o {params} {input}
    """

# trimmomatic "trims up" Illumina seq data by:
#     - removing adaptors
#     - removing low quality bases
#     - remove reads less than 36 bp long
rule trimming:
  input:
    "data/{smp}.fastq.gz"
  output:
    "results/data_trimmed/{smp}.trimmed.fastq.gz"
  params:
    "results/data_trimmed/{smp}.trimmed.out"
  conda:
    "envs/trimmomatic.yaml"
  threads: 5
  shell:
    """
    trimmomatic SE -threads {threads} -phred33 {input} {output} ILLUMINACLIP:{ADAPTERS}:2:30:10 SLIDINGWINDOW:4:20 MINLEN:{CROPSIZE} &> {params}
    """

# Generates a report about the quality of the trimmed data
rule qc_trimmed:
  input:
    "results/data_trimmed/{smp}.trimmed.fastq.gz"
  output:
    "results/qc_trimmed/{smp}.trimmed_fastqc.html"
  params:
    "results/qc_trimmed"
  conda:
    "envs/fastqc.yaml"
  shell:
    """
    fastqc -o {params} {input}
    """

# Unsure about this rule
rule index:
  input:
    REFERENCE
  output:
    temp(REFERENCE+".index")
  shell:
    """
    mrsfast --index {REFERENCE} --ws {INDEXWINDOW}
    """

# Generates all possible alignments for the reads
# Allows for multiple reads
rule mrsfast:
  input:
    "results/data_trimmed/{smp}.trimmed.fastq.gz",
    REFERENCE+".index"
  output:
    aligned=temp("results/alignments/{smp}.trimmed.aligned.sam"),
    unaligned=temp("results/alignments/{smp}.trimmed.unaligned")
  params:
    "results/alignments/{smp}.trimmed.aligned.out"
  threads: 5
  shell:
    """
    mrsfast --search {REFERENCE} --seqcomp --seq {input} -e 2 -o {output.aligned} -u {output.unaligned} --crop {CROPSIZE} --threads {threads} &> {params}
    """

# Converts output of MrsFast analysis from a SAM file to a BAM file
# SAM file = human readable, BAM file = binary equivalent
rule sam2bam:
  input:
    "results/alignments/{smp}.trimmed.aligned.sam"
  output:
    temp("results/alignments/{smp}.trimmed.aligned.bam")
  conda:
    "envs/samtools.yaml"
  shell:
    """
    samtools view -Sb {input} > {output}
    """

# alignments generated from MrsFast are in random order
# this sorts them into genome order
# required for input into CSEM
rule sorting:
  input:
    "results/alignments/{smp}.trimmed.aligned.bam"
  output:
    "results/alignments/{smp}.trimmed.aligned.sort-n.bam"
  conda:
    "envs/samtools.yaml"
  shell:
    """
    samtools sort -n {input} > {output}
    """

# CSEM -- Decides the most likely alignment based on variants between paralogs
# gives probability of each alignment being correct
rule csem:
  input:
    "results/alignments/{smp}.trimmed.aligned.sort-n.bam"
  output:
    "results/allocations/{smp}.trimmed.aligned.sort-n.csem.bam"
  params:
    "results/allocations/{smp}.trimmed.aligned.sort-n.csem"
  threads: 5
  shell:
    """
    run-csem --bam -p {threads} {input} 200 {params}
    """

# Unsure about this rule
rule sampling:
  input:
    "results/allocations/{smp}.trimmed.aligned.sort-n.csem.bam"
  output:
    temp("results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.header"),
    temp("results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.sam"),
    "results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.bam"
  params:
    "results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling"
  conda:
    "envs/samtools.yaml"
  shell:
    """
    bash scripts/sampling.sh {input} {params}
    """

# MACS performs peak calling: decides which reads are actually peaks
# outputs the visualization of the peaks
# Calls MACS on alignments that haven't been put through CSEM yet
rule macs:
  input:
    samples = expand("results/alignments/{smp}.trimmed.aligned.sort-n.bam", smp = SAMPLES),
    control = expand("results/alignments/{smp}.trimmed.aligned.sort-n.bam", smp = CONTROL)
  output:
    "results/peak_calling/"+FILENAME+"_control_lambda.bdg",
    "results/peak_calling/"+FILENAME+"_treat_pileup.bdg"
  params:
    samples = " ".join(expand("results/alignments/{smp}.trimmed.aligned.sort-n.bam", smp = SAMPLES)),
    control = " ".join(expand("results/alignments/{smp}.trimmed.aligned.sort-n.bam", smp = CONTROL)),
    outdir = "results/peak_calling"
  conda:
    "envs/macs2.yaml"
  shell:
    """
    macs2 callpeak -t {params.samples} -c {params.control} -f BAM -g hs --outdir {params.outdir} -n {FILENAME} -B --nomodel --extsize 400 -q 1e-2 {PEAK_TYPE}
    """

# Call MACS on the CSEM output?
rule macs_csem:
  input:
    samples = expand("results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.bam", smp = SAMPLES),
    control = expand("results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.bam", smp = CONTROL)
  output:
    "results/peak_calling/"+FILENAME+"_csem_control_lambda.bdg",
    "results/peak_calling/"+FILENAME+"_csem_treat_pileup.bdg"
  params:
    samples = " ".join(expand("results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.bam", smp = SAMPLES)),
    control = " ".join(expand("results/allocations/{smp}.trimmed.aligned.sort-n.csem.sampling.bam", smp = CONTROL)),
    outdir = "results/peak_calling"
  conda:
    "envs/macs2.yaml"
  shell:
    """
    macs2 callpeak -t {params.samples} -c {params.control} -f BAM -g hs --outdir {params.outdir} -n {FILENAME}_csem -B --nomodel --extsize 400 -q 1e-2 {PEAK_TYPE}
    """

# generate fold-enrichment and logLR track for the peaks
# fold-enrichment = normalization method, divides CHIP signals by no Ab control signals
# logLR track =
rule bdgcmp:
  input:
    sample = "results/peak_calling/{exp}_treat_pileup.bdg",
    control = "results/peak_calling/{exp}_control_lambda.bdg"
  output:
    "results/peak_calling/{exp}_logLR.bdg"
  params:
    "results/peak_calling/"
  conda:
    "envs/macs2.yaml"
  shell:
    """
    macs2 bdgcmp -t {input.sample} -c {input.control} --outdir {params} --o-prefix {wildcards.exp} -m logLR -p 0.00001
    """

# Converts peaks from BED file to a bigwig file
rule bigwig:
   input:
     "results/peak_calling/{exp}_logLR.bdg"
   output:
     "results/peak_calling/{exp}_logLR.bw"
   shell:
     """
     bedGraphToBigWig {input} {CHROMSIZE} {output}
     """
