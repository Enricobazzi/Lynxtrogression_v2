import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list
import seaborn as sns
import matplotlib.pyplot as plt
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Plot a window matrix from a vcf file")
    parser.add_argument("--ivcf", help="The input vcf file")
    parser.add_argument("--chr", help="The chromosome to extract")
    parser.add_argument("--start", type=int, help="Start position")
    parser.add_argument("--end", type=int, help="End position")
    parser.add_argument("--oplot", type=str, help="The output plot file")
    return parser.parse_args()

args = parse_args()

ivcf = args.ivcf
chr = args.chr
start = int(args.start)
end = args.end
oplot = args.oplot

names = []
sequences = {}
positions = []

# parse the vcf file
with open(ivcf) as file:
    for line in file:
        if line.startswith("#CHROM"):
            line = line.strip().split()
            for name in line[9:]:
                names.append(f'{name}')
        elif not line.startswith("#"):
            line = line.strip().split()
            if line[0] != chr:
                continue
            pos, ref, alt = line[1], line[3], line[4]
            if int(pos) < start or int(pos) > end:
                continue
            positions.append(pos)
            for i in range(9, len(line)):
                name1 = f'{names[i-9]}.1'
                name2 = f'{names[i-9]}.2'
                if name1 not in sequences:
                    sequences[name1] = []
                if name2 not in sequences:
                    sequences[name2] = []
                if line[i] == "0|0":
                    sequences[name1].append(0)
                    sequences[name2].append(0)
                elif line[i] == "0|1":
                    sequences[name1].append(0)
                    sequences[name2].append(1)
                elif line[i] == "1|0":
                    sequences[name1].append(1)
                    sequences[name2].append(0)
                elif line[i] == "1|1":
                    sequences[name1].append(1)
                    sequences[name2].append(1)

# build dataframe from sequences with a column for individuals and one column for each position
df = pd.DataFrame(sequences)
df = df.T

# filter for allele frequency > 0.05
# df = df.loc[:, (df.sum(axis=0) / len(df) > 0.05) & (df.sum(axis=0) / len(df) < 0.95)]

# === Compute Clustering ===
distance_matrix = pdist(df.values, metric='hamming')
linkage_matrix = linkage(distance_matrix, method='average')
leaf_order = leaves_list(linkage_matrix)
ordered_df = df.iloc[leaf_order]
ordered_names = ordered_df.index.tolist()

# === Plot Layout ===
fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 4, 0.5], wspace=0.05)

# === Plot Dendrogram (Left) ===
ax0 = fig.add_subplot(gs[0])
dendro = dendrogram(linkage_matrix, orientation='left', labels=ordered_names, ax=ax0, color_threshold=0)
ax0.invert_yaxis()
ax0.axis('off')

# === Force rendering to get leaf coordinates ===
fig.canvas.draw()

# Map from reordered leaf index to actual Y position on the plot
# dendro['leaves'] gives row indices from the original DataFrame
# ax0.collections[0] is the dendrogram lines collection
# Get the y positions directly from tick positions after plot is drawn
y_positions = ax0.get_yticks()
leaves = dendro['leaves']  # these are indices into your original DataFrame

for i, leaf_index in enumerate(leaves):
    sample_name = df.index[leaf_index]
    y = y_positions[i]

    # Population color
    if 'sm' in sample_name.lower():
        color = '#8d56a4'
    elif 'ca' in sample_name.lower():
        color = '#1b9773'
    elif 'ki' in sample_name.lower() or 'ur' in sample_name.lower():
        color = '#003c82'
    else:
        color = 'gray'

    # Plot circle at left edge of leaf line
    ax0.plot(0.01, y, 'o', markersize=8, color=color, markeredgecolor='k')

# === Heatmap (Middle) ===
ax1 = fig.add_subplot(gs[1])
sns.heatmap(ordered_df, cmap=["white", "black"], cbar=False, ax=ax1,
            xticklabels=False, yticklabels=False)
ax1.set_xlabel("")
ax1.set_ylabel("")

# === Labels + Colored Circles (Right) ===
ax2 = fig.add_subplot(gs[2])
ax2.set_xlim(0, 1)
ax2.set_ylim(0, len(ordered_names))
ax2.axis('off')

for i, name in enumerate(ordered_names):
    y = len(ordered_names) - i - 0.5
    ax2.text(0.1, y, name, va='center', fontsize=9)

plt.savefig(oplot, bbox_inches='tight')