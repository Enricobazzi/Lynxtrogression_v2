# Lynxtrogression_v2

## Detecting Genomic Introgressions Among Lynx Lineages

In this repository you can find all of the workflows and scripts I wrote and ran in order to detect introgressed windows in the genomes of *Lynx pardinus* and *Lynx lynx* populations.

### Selecting Samples

The genomic dataset we used for the introgression scans was comprised of all of the individuals we have sampled from the following populations:

* *Lynx pardinus* :
  * Sierra Morena (N=29) = **lpa**

* *Lynx lynx* :
  * Western clade Kirov + Urals (N=20) = **wel**
  * Eastern clade Yakutia + Primorsky krai (N=19) = **eel**
  * Southern clade Caucasus (N=12 or N=9) = **sel**

The [populations table](data/lp_ll_introgression_populations.txt) has information on how to convert sample names to their population

| sample | population |
|:------:|:----------:|
| sm     | lpa        |
| ki     | wel        |
| ur     | wel        |
| ya     | eel        |
| vl     | eel        |
| ca     | sel        |

A [table](data/samples_table.xlsx) with each sample's species, population, original study, sequencing technology, median read depth can be found in the data folder. There I also marked which samples were used in which analysis.

### Alignment to reference genomes

I align the sequencing reads to the newly generated *Lynx rufus* reference genome ([mLynRuf2.2](https://denovo.cnag.cat/lynx_rufus)). Description of how sequencing reads were processed and aligned are found in [reads_qc_and_alignment](reads_qc_and_alignment.md).

### Variant Calling from aligned reads

Genotypes of all the samples were extracted from the individual alignments and joint into a unique VCF file as described in [variant_calling](variant_calling.md)
