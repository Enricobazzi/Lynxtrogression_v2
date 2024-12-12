# Lynxtrogression_v2

## Detecting Genomic Introgressions Among Lynx Lineages

In this repository you can find all of the workflows and scripts I wrote and ran in order to detect introgressed windows in the genomes of *Lynx pardinus* and *Lynx lynx* populations.

### Selecting Samples

I describe how samples are selected in [sample_selection](sample_selection.md)

### Alignment to reference genomes

I align the sequencing reads to the newly generated *Lynx rufus* reference genome ([mLynRuf2.2](https://denovo.cnag.cat/lynx_rufus)). Description of how sequencing reads were processed and aligned are found in [reads_qc_and_alignment](reads_qc_and_alignment.md).

### Variant Calling from aligned reads

Genotypes of all the samples were extracted from the individual alignments and joint into a unique VCF file as described in [variant_calling](variant_calling.md)

### Variant Filtering

I describe the steps I take for filtering the variant dataset in [variant_filtering](variant_filtering.md).

### Demographic inference

I reconstructed the demographic history between population pairs using the software [GADMA2](https://github.com/ctlab/GADMA). All the necessary steps are described in [demographic_inference](demographic_inference.md).

### Introgression scans

I used an approach based on deep convolutional neural networks, described by [Ray et al. 2024](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010657) and implemented in [introNets](https://github.com/SchriderLab/introNets), to identify windows with introgression along the genomes of all the lynx populations (*eel not included for now*). Data preparation, scripts and description of downstream analyses and results are included in [introgression_scans_withM](introgression_scans_withM.md).


### Conda environment for analyses

Unless specified this environment is used and created as follows:
```
conda create -n lynxtrogression_v2 python=3.9
conda activate lynxtrogression_v2
pip install demes
pip uninstall ruamel.yaml
pip install "ruamel.yaml<0.18.0"
pip install pandas
conda install -c conda-forge dadi
pip install -U scikit-learn
pip install -U matplotlib
pip install seaborn
```
