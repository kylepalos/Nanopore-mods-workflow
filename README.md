# Nanopore-mods-workflow

## Updated April 2025
This repository documents my current workflow for calling RNA modifications from Nanopore direct RNA-sequencing data (using the current SQK-RNA004 [kit](https://store.nanoporetech.com/us/direct-rna-sequencing-kit-004.html)).

# STEP 1:
## Basecall
Download the latest release of [Dorado](https://github.com/nanoporetech/dorado/releases) basecaller (I'm currently using v0.8.0 - which calls m6A, pseudoUridine, inosine, and m5C).

Then - run Dorado basecaller on POD5 files. Note: Dorado does not support the old version direct RNA-sequencing kits (e.g., SQK-RNA001).

```
/dorado-0.8.0-linux-x64/bin/dorado basecaller \
-v hac,inosine_m6A,m5C,pseU  \
--min-qscore 10 \
--emit-moves \
path/to/pod5_files/ \
--estimate-poly-a \
--device "cuda:0" > file.bam
```

Explanation of parameters used:

-v hac,inosine_m6A,m5C,pseU: This specifies the model speed used to basecall POD5 files (there are: fast, hac, and sup. [See here.]([url](https://github.com/nanoporetech/dorado#model-selection-foreword)) Additionally, the following modification options will detect and add the necessary SAM headers to retain modification information.

--min-qscore 10: The minimum -log10 quality score as described by Illumina [here.](https://help.basespace.illumina.com/files-used-by-basespace/quality-scores)

--emit-moves: Not super necessary to have, but this options retains the raw signal information of the sequencing read, [see here.](https://github.com/nanoporetech/dorado/issues/110)

--estimate-poly-a: Estimate the polyA tail length for each read

--device "cuda:0": Specifies the GPU to basecall.


## Automating Dorado for several input projects across different directories.
I often don't just have a single sample to basecall, but an entire experiment (i.e., replicates and condition - several samples.)

I have a script to automate your Dorado run with input directories and output file names in the **scripts** folder.


## Convert BAM to fastq
Dorado can perform the genome alignment while basecalling but I prefer to customize the parameters myself. Therefore, convert the BAM to fastq.

```
samtools fastq -@ 5 -T "*" file.bam | gzip > file.fastq.gz
```

The important part of this step is:

-T "*" which transfers all of the header tags to the fastq headers (modification information). 

**NOTE:** You should skip this step and just use ```--emit-fastq``` when running Dorado. I will update my workflow with these changes in the future.


# STEP 2:
## Map using Minimap2
Map to genome using Minimap2:

```
minimap2 -ax splice \
-uf -k14 \
-G 10000 \
-y --secondary=no \
genome.fasta file.fastq.gz | \
samtools sort -o file_mapped.bam -
```

Some of these parameters, I'm getting from [here.](https://github.com/nanoporetech/dorado/issues/145)

-ax splice: Preset mode for spliced long reads (spliced because we're aligning (mostly) spliced poly-A RNAs to the genome, not transcriptome).

-uf -k14: Settings recommended by Heng Li (creator of Minimap2 & SamTools) for noisy direct RNA-seq reads. -k 14 specifies the k-mer size, I'm not sure why 14 is optimal here. -uf **sort of** specifies the strandedness of the dataset by specifying the orientation of how to find canonical splice sites (the "f" of uf specifies transcript strand).

-G 10000: Maximum intron length - use a shorter value for plants.

-y: Transfer fastq tags to output BAM tags (necessary for retaining modification information)

--secondary=no: Long reads shouldn't have too many secondary aligments - and if they do - the modification calls may not be too inspiring.


# STEP 3:
## Call modifications using ModKit
Modkit is Nanopore's software for extracting modification information from BAM files and generating BED files with necessary statistics, [see here.](https://github.com/nanoporetech/modkit)

```
modkit pileup \
--mod-thresholds a:0.99 --mod-thresholhds m:0.99 --mod-thresholhds 17802:0.99 --mod-thresholds 17596:0.99 \
file_mapped.bam \
mod_pileup.bed
```

The pileup command extracts modified sites in [bedmethyl](https://nanoporetech.github.io/modkit/intro_pileup.html#description-of-bedmethyl-output) format.

The mod thresholds essentially deal with the probability thresholds assigned to the basecallers confidence in the modification of the site.

a = m6A, m = m5C, 17802 = pseudouridine, 17596 = inosine.

If you want more information on the thresholds, see an github issue I submitted [here](https://github.com/nanoporetech/modkit/issues/198) as well as Modkit's [documentation](https://nanoporetech.github.io/modkit/filtering_details.html).


# Threshold considerations:
This is current as of April, 2025.

I use a mod threshold for every modification type == 0.99 because, for my data, this seems to be the best balance of stringency and comprehensiveness. See plot below:

![mod_probs](https://github.com/user-attachments/assets/dd47e2d4-d617-4620-a0d6-514bc489a114)

This is a histogram of modification probabilities for all 4 RNA mod types currently available with Dorado. As you can see, there is a spike in fraction of calls at very high probabilities - therefore, setting a threshold at 0.99 (the vertical black line) should capture a good fraction of the data while reducing false positives (except for m5C currently.)


## For library prep:
SQK-RNA004 allows you to use total or poly-A enriched RNA for library prep. We wanted to whether total RNA input would match poly-A enriched RNA in terms of quantification.

The following plot shows isoform-level "counter per million" from a matched sample that was split between total RNA for library prep and poly-A enrichment before library prep.

Lower abundant transcripts naturally have more variability so I would call this pretty good.

![cpm_comparison](https://github.com/kylepalos/Nanopore-mods-workflow/assets/56089443/f5caf43a-b9bd-416d-8c57-2ccf1ddfc18c)
