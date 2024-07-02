# Nanopore-mods-workflow

## Updated July 2024
This repository documents my current workflow for calling RNA modifications from Nanopore direct RNA-sequencing data (using the current SQK-RNA004 [kit](https://store.nanoporetech.com/us/direct-rna-sequencing-kit-004.html)).

### Basecall
Download the latest release of [Dorado](https://github.com/nanoporetech/dorado/releases) basecaller (I'm currently using v0.7.0 - which calls m6A and pseudoUridine).

Then - run Dorado basecaller on POD5 files. Note: Dorado does not support the old version direct RNA-sequencing kits (e.g., SQK-RNA001).

```
/dorado-0.7.0-linux-x64/bin/dorado basecaller \
-v hac,m6A \
--min-qscore 10 \
--emit-moves \
path/to/pod5_files/ \
--estimate-poly-a \
--device "cuda:0" > file.sam
```

Explanation of parameters used:

-v hac, m6a: This specifies the model speed used to basecall POD5 files (there are: fast, hac, and sup. [See here.]([url](https://github.com/nanoporetech/dorado#model-selection-foreword)) Additionally, the "m6a" will let the model include the necessary SAM headers to retain modification information.

You could plug in "psU" instead of m6A to call pseudoUridine.

--min-qscore 10: The minimum -log10 quality score as described by Illumina [here.](https://help.basespace.illumina.com/files-used-by-basespace/quality-scores)

--emit-moves: Not super necessary to have, but this options retains the raw signal information of the sequencing read, [see here.](https://github.com/nanoporetech/dorado/issues/110)

--estimate-poly-a: Estimate the polyA tail length for each read

--device "cuda:0": Specifies the GPU to basecall.

### Convert SAM to fastq
Dorado outputs reads in SAM format, but I prefer to have reads in fastq format.

```
samtools fastq -@ 5 -T "*" file.sam > file.fastq
```

The important part of this step is:

-T "*" which transfers all of the SAM tags to the fastq headers (modification information). 

### Map using Minimap2
Map to genome using Minimap2:

```
minimap2 -ax splice \
-uf -k14 \
-G 10000 \
-y --secondary=no \
genome.fasta file.fastq | \
samtools sort -o file_mapped.bam -
```

Some of these parameters, I'm getting from [here.](https://github.com/nanoporetech/dorado/issues/145)

-ax splice: Preset mode for spliced long reads (spliced because we're aligning (mostly) spliced poly-A RNAs to the genome, not transcriptome).

-uf -k14: Settings recommended by Heng Li (creator of Minimap2 & SamTools) for noisy direct RNA-seq reads. -k 14 specifies the k-mer size, I'm not sure why 14 is optimal here. -uf **sort of** specifies the strandedness of the dataset by specifying the orientation of how to find canonical splice sites (the "f" of uf specifies transcript strand).

-G 10000: Maximum intron length - use a shorter value for plants.

-y: Transfer fastq tags to output BAM tags (necessary for retaining modification information

--secondary=no: Long reads shouldn't have too many secondary aligments - and if they do - the modification calls may not be too inspiring.

### Call modifications using ModKit
Modkit is Nanopore's software for extracting modification information from BAM files and generating BED files with necessary statistics, [see here.](https://github.com/nanoporetech/modkit)

```
modkit pileup \
--filter-threshold A:0.8 --mod-thresholds a:0.99 \
--motif RRACH 2 \
--ref genome.fasta file_mapped.bam \
m6a_pileup.bed
```

The pileup subcommand runs through all sites or motifs and generates a probability that the site is modified.

Filter threshold and mod threshold are parameters that I struggle to explain succinctly, but they deal with the probability thresholds that are necessary for calling modifications at the base in question. Note that filter threshold specifies a capital "A" while mod threshold specifies a lower-case "a". This is the nomenclature that ModKit uses to distinguish an unmodied A from a modified A, respectively. [See here](https://github.com/nanoporetech/modkit/issues/198) where I discuss with the authors on how to arrive at these values.

--motif RRACH 2 specifies to only scan RRACH motifs in the genome file provided and the "2" is the 0-based offset to the base in question (the A at the center of RRACH is 2 bases offset).

## For library prep:
SQK-RNA004 allows you to use total or poly-A enriched RNA for library prep. We wanted to whether total RNA input would match poly-A enriched RNA in terms of quantification.

The following plot shows isoform-level "counter per million" from a matched sample that was split between total RNA for library prep and poly-A enrichment before library prep.

Lower abundant transcripts naturally have more variability so I would call this pretty good.

![cpm_comparison](https://github.com/kylepalos/Nanopore-mods-workflow/assets/56089443/f5caf43a-b9bd-416d-8c57-2ccf1ddfc18c)
