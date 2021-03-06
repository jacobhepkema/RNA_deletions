# RNA_deletions

This repository contains scripts accompanying [a recent biorxiv preprint](https://www.biorxiv.org/content/10.1101/2021.09.15.460454v1).

The required data to run the scripts is available here:
* ALKBH5: SRA accession SRR9646144 (use `fasterq-dump SRR9646144` in directory `data` to get FASTQ files). `fasterq-dump` is a tool [available in the SRA toolkit](https://github.com/ncbi/sra-tools).
* ChAdOx1 HEK293: SRA accession SRR13320597 (use `fasterq-dump SRR13320597` in directory `data` to get FASTQ files)
* S protein constructs: download [here (data coming soon)](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA762212)

# Workflows

To generate the BAM files for downstream analysis, run the following:

## ALKBH5 workflow

The following will run minimap2 and convert the SAM to a BAM file (requires minimap2, samtools):

```
./ALKBH5.sh
```
Outputs will be in the `outputs/ALKBH5` directory. The script uses `-t 8` by default, edit the file to use a different number of processors.

Then, to run the downstream analysis, you need to run 
```
Rscript ALKBH5.R
```
from an environment that contains R & packages `GenomicRanges`, `Biostrings`, `Rsamtools`, `BSgenome`, and `stringr`. 
If you need this environment you can get a Singularity image from [this repository (coming soon)](https://github.com/jacobhepkema/RNA_deletions_R). If you pull the Singularity image, you can then run (after installing [Singularity](https://sylabs.io/guides/2.6/user-guide/index.html)):
```
singularity exec image_name.sif Rscript ALKBH5.R
```
where you replace `image_name.sif` with the filename of the container.

## ChAdOx1 workflow 

Similar to above, first run this script to run minimap2 & samtools:

```
./ChAdOx1.sh
```
Outputs will be in the `outputs/ChAdOx1` directory. The script uses `-t 8` by default, edit the file to use a different number of processors.

Then, to run the downstream analysis, run
```
Rscript ChAdOx1.R
```
from the same environment [described above](#alkbh5-workflow), the Singularity environment will work here too.


## S protein constructs workflow

This uses a small Nextflow pipeline. Download and install Nextflow [using the instructions on their website](https://www.nextflow.io/).
Then, run the following in an environment that contains minimap2, samtools:

```
nextflow S_protein.nf
```
If you want to run it locally, remove the `nextflow.config` file, or edit the `nextflow.config` file [accordingly](https://www.nextflow.io/docs/latest/config.html). If you want to run in manually, check the contents of the processes in `S_protein.nf` - there should be enough information in that to run minimap2 and samtools from the command line. 

Then, to run the downstream analysis, run
```
Rscript S_protein.R
```
from the same environment [described above](#alkbh5-workflow), the Singularity environment will work here too.


# Explanation of scripts

The scripts in `deletion_utils.R` are made to identify small deletions/introns using the CIGAR strings of the alignment BAM file. 

Usually, this requires Oxford Nanopore direct RNA or cDNA reads aligned using e.g. [minimap2](https://github.com/lh3/minimap2) as such:
```
minimap2 -ax splice -t8 --cs=long --sam-hit-only template.fa reads.fastq > mapped.sam
```

The steps for using the functions are usually as follows:

1. Read BAM file using [R package Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) and convert to GRanges using function `BAM_to_granges()`
2. Identify deletions in reads using `get_deletion_df_no_overlaps()`. 
3. Optionally, retain only deletions whose sequences are diverse by calculating sequence diversity with `seq_diverse()` (this filters out many deletions that might be due to poor mapping)
4. Optionally, get information about read span on template using `get_spans()`.

### `get_deletion_df_no_overlaps()`

Briefly, this method works by cycling through aligned reads, considering all CIGAR string positions that are not I, S or H (no insertions, no soft/hard clipping) so that you only look at positions mapped to the template (mapped & introns & deletions). It will assess all CIGAR strings that are D or N, and check if these events are 1) not the first or last event in what's left of the CIGAR string and 2) above a specified minimum length (very small deletions are likely noise)

There is a bit of a different workflow if you only want to consider deletions that are contained within a specified BED file (e.g. deletions within exons), for this you use the following function:

### `get_deletion_df()`

Briefly, this method works like the other function, but will only consider reads that overlap other regions (e.g. exons) (calculated using the `findOverlaps()` function from R package `GenomicRanges`), and it will make sure there are other-region-mapping reads flanking the deletion (on the _same_ 'other region').


### Used packages

`stringr`, `Biostrings`, `GenomicRanges`, `Rsamtools`, optionally also e.g. `BSgenome.Hsapiens.UCSC.hg38` if you aligned to `hg38`. 
