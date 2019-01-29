# Snakemake file to call peaks in ChIP-Seq experiments with allocation step

## 1) Required software

Besides `snakemake`, this pipeline uses the following tools:
- fastqc
- trimmomatic
- samtools
- mrsfast
- macs2
- csem
- bedGraphToBigWig

Most tools are automatically installed by `snakemake` through Conda. MrsFast, CSEM and BedGraphToBigWig must be installed locally (they are not yet available in Conda repositories). 

## 2) Folder structure

The following structure is recommended:

```
|__ adapters
    |__ sequences.fa 
|__ config.yaml # Mandatory: described in section 3
|__ data 
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
    |__ sampling.sh
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

adapters: adapters/sequences.fa

# Reference path

reference: reference/masked_reference.fa
chromsize: reference/chromosome_sizes.txt

# mrsFast options

indexwindow: integer
cropsize: integer

# Peak calling options

broad: boolean

# Workflow options

csem: boolean
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

