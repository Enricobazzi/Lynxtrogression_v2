import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import numpy as np

def read_bed(intro):
    ibed = f"data/introgression_scans/bed_files/{intro}_intro.merged.t_dist.bed"
    bed = pd.read_csv(ibed, sep='\t', names=["chrom", "start", "end", "t_dist"])
    bed["length"] = bed["end"] - bed["start"]
    bed["n_genes"] = pd.read_csv(f"data/introgression_scans/bed_files/{intro}_intro.merged.n_genes.bed", 
                                 sep='\t', names=["chrom", "start", "end", "n_genes"])["n_genes"]
    bed["gene_overlap"] = pd.read_csv(f"data/introgression_scans/bed_files/{intro}_intro.merged.gene_overlap.bed", 
                             sep='\t', names=["chrom", "start", "end", "gene_overlap"])["gene_overlap"]
    bed["cds_overlap"] = pd.read_csv(f"data/introgression_scans/bed_files/{intro}_intro.merged.cds_overlap.bed", 
                             sep='\t', names=["chrom", "start", "end", "cds_overlap"])["cds_overlap"]
    bed["diversity"] = pd.read_csv(f"data/introgression_scans/bed_files/{intro}_intro.merged.diversity.bed",
                                sep='\t', names=["chrom", "start", "end", "diversity"])["diversity"]
    return bed

def read_ran_bed(intro: str, n: int):
    ibed = f"data/introgression_scans/{intro}_randomwins/random_windows.{n}.t_dist.bed"
    bed = pd.read_csv(ibed, sep='\t', names=["chrom", "start", "end", "t_dist"])
    bed["length"] = bed["end"] - bed["start"]
    bed["n_genes"] = pd.read_csv(f"data/introgression_scans/{intro}_randomwins/random_windows.{n}.n_genes.bed", 
                                 sep='\t', names=["chrom", "start", "end", "n_genes"])["n_genes"]
    bed["gene_overlap"] = pd.read_csv(f"data/introgression_scans/{intro}_randomwins/random_windows.{n}.gene_overlap.bed", 
                                 sep='\t', names=["chrom", "start", "end", "gene_overlap"])["gene_overlap"]
    bed["cds_overlap"] = pd.read_csv(f"data/introgression_scans/{intro}_randomwins/random_windows.{n}.cds_overlap.bed", 
                                 sep='\t', names=["chrom", "start", "end", "cds_overlap"])["cds_overlap"]
    bed["diversity"] = pd.read_csv(f"data/introgression_scans/{intro}_randomwins/random_windows.{n}.diversity.bed",
                                sep='\t', names=["chrom", "start", "end", "diversity"])["diversity"]
    return bed

lpa_to_wel = read_bed("lpa_to_wel")
lpa_to_sel = read_bed("lpa_to_sel")
wel_and_sel_to_lpa = read_bed("wel_and_sel_to_lpa")
wel_to_lpa = read_bed("wel_to_lpa")
sel_to_lpa = read_bed("sel_to_lpa")

# draw boxplots of introgressed windows length
plt.figure(figsize=(10, 5))

boxprops = dict(linestyle='-', linewidth=1, color='black')
medianprops = dict(linestyle='-', linewidth=1, color='red')

plt.subplot(1, 2, 1)
bp1 = plt.boxplot([wel_to_lpa["length"], sel_to_lpa["length"], wel_and_sel_to_lpa["length"]], 
            labels=["ELw to ILa", "ELs to ILa", "ELw+ELs to ILa"], showfliers=False,
            boxprops=boxprops, medianprops=medianprops, patch_artist=False)
plt.ylabel("Length (bp)")
plt.title("Iberian lynx\nIntrogressed windows length")

plt.subplot(1, 2, 2)
bp2 = plt.boxplot([lpa_to_wel["length"], lpa_to_sel["length"]], 
            labels=["ILa to ELw", "ILa to ELs"], showfliers=False,
            boxprops=boxprops, medianprops=medianprops, patch_artist=False)
plt.ylabel("Length (bp)")
plt.title("Eurasian lynx\nIntrogressed windows length")

plt.tight_layout()
plt.savefig("plots/introgression_scans/introgressed_windows_length.boxplots.pdf", bbox_inches='tight')

