import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

ibeds=[
    'data/introgression_scans/bed_files/wel_and_sel_to_lpa_intro.merged.bed',
    'data/introgression_scans/bed_files/lpa_to_wel_intro.merged.bed',
    'data/introgression_scans/bed_files/lpa_to_sel_intro.merged.bed'
]

dfs = []

for ibed in ibeds:
    bed = pd.read_csv(ibed, sep='\t', names=["chrom", "start", "end"])
    intro = ibed.split("/")[-1].split(".")[0]
    chrs = bed.chrom.unique()
    for n, chrom in enumerate(chrs[::-1]):
        wins = bed[bed.chrom == chrom]
        nin = np.array([wins['end'] - wins['start']]).sum()
        ntot = max(wins["end"])
        pc = f'{round(nin/ntot*100, 2)}'
        df=pd.DataFrame({
            "introgression": intro,
            "chrom": chrom,
            "chrom_length": ntot,
            "introgression_pc": pc
        }, index=[0])
        dfs.append(df)

df = pd.concat(dfs)
df["introgression_pc"] = df["introgression_pc"].astype(float)
df["chrom_length"] = df["chrom_length"].astype(float)
df["introgression"] = df["introgression"].astype(str)
df["chrom"] = df["chrom"].astype(str)

# calculate correlation coefficients
corrs = df.groupby("introgression").apply(lambda x: x["chrom_length"].corr(x["introgression_pc"]))
# add a long name to the introgressions in corrs
corrs.index = ["Southern Eurasian lynx", "Western Eurasian lynx", "Iberian lynx"]
print(corrs)

fig, ax = plt.subplots(figsize=(6, 6))
sns.scatterplot(data=df, x="chrom_length", y="introgression_pc", hue="introgression", ax=ax)
# add tendency line for each introgression
for intro in df["introgression"].unique():
    sub = df[df["introgression"] == intro]
    sns.regplot(data=sub, x="chrom_length", y="introgression_pc", scatter=False, ax=ax)
# add correlation coefficients to the legend so that 
# the legend shows the correlation coefficient next to the name of each introgression
handles, labels = ax.get_legend_handles_labels()
labels = [
    f'Iberian lynx = {str(round(corrs["Iberian lynx"], 2))}',
    f'Western Eurasian lynx = {str(round(corrs["Western Eurasian lynx"], 2))}',
    f'Southern Eurasian lynx = {str(round(corrs["Southern Eurasian lynx"], 2))}',
]

ax.legend(handles, labels, title="Correlation (r)", loc="upper right")
ax.set_xlabel("Chromosome length")
ax.set_ylabel("Introgression amount (%)")
plt.savefig('plots/introgression_scans/intro_genomes/chr_len_vs_introgression.pdf')
