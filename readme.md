# Snakemake file to call peaks in ChIP-Seq experiments with allocation step

## 1) Required software

Besides `snakemake`, this pipeline uses the following tools:
- fastqc
- trimmomatic
- samtools
- bowtie (short read) bowtie2 (long read) 
- macs2
- csem
- rsem (long chip pipeline only) 
- bedGraphToBigWig
- phantompeakqualtools

Most tools are automatically installed by `snakemake` through Conda. CSEM and BedGraphToBigWig must be installed locally (they are not yet available in Conda repositories). 

## 2) Folder structure

The following structure is recommended:

```
|__ adapters
    |__ sequences.fa 
|__ config.yaml # Mandatory: described in section 3
|__ reads 
    |__ sample1.fastq 
    |__ sample2.fastq
    |__ control1.fastq
    |__ ...
|__ envs
    |__ fastqc.yaml
    |__ macs2.yaml
    |__ samtools.yaml
    |__ trimmomatic.yaml
|__ readme.md
|__ reference
    |__ masked_reference.fa
    |__ chromosome_size.txt
|__ scripts
    |__ sampling.sh or long_sampling.sh
|__ Snakefile
```

The pipeline will generate `results` directory containing intermediate files and peaks called. 

## 3) Config file

Config file must contain the following information:

```
# Analysis name

filename: example

# Input files

samples:
  - sample1
  - sample2
  - ...

control:
  - control1
  - ...
  
samples_size: # from phantompeakqualtools
  - 485
  - 485

control_size: # from phantompeakqualtools
  - 665
  - 655

adapters: adapters/sequences.fa

# Reference path

reference: reference/masked_reference.fa
chromsize: reference/chromosome_sizes.txt

# Peak calling options

broad: boolean

# Workflow options

csem: boolean
```

In order to obtain the values for samples_size and control_size in the config file, we need to predict fragment size for each dataset using the alignments obtained with bwt-default and PhantomPeakQualTools. To do this, we need to run the Snakefile until sam2bam rule and then run phantompeaktools. After obtaining the sizes from phantompeaktools, update the samples_size and control_size values in the config file and then run the snakemake as below. 

For example:

```
snakemake -s Snakefile --configfile GM12878_H3K27ac.yaml -p -j 10 --until sam2bam
Rscript /share/dennislab/programs/phantompeakqualtools/run_spp.R -c=results/alignments/ENCFF001CUR.trimmed.bwt-df.bam
```

## 4) Execution

Export locally installed software to your PATH:
```
export PATH="/share/dennislab/programs/share_path:/share/dennislab/programs/csem-2.4:$PATH"
```

Dry-run:
```
snakemake --configfile config.yaml --use-conda -n
```

Obtain workflow graphs:
```
snakemake --configfile config.yaml --dag | dot -Tpng > dag.png 
snakemake --configfile config.yaml --rulegraph | dot -Tpng > rulegraph.png
```

Run interactively:
```
snakemake --configfile config.yaml --use-conda -j 20 -p
```

