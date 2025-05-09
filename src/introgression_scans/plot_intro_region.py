import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

def segment_max_intervals(df):
    """
    Given a DataFrame with overlapping genomic intervals (chrom, start, end, value),
    this function returns non-overlapping segments with the maximum value
    from any overlapping intervals.

    Parameters:
        df (pd.DataFrame): Input DataFrame with columns ["chrom", "start", "end", "value"]

    Returns:
        pd.DataFrame: Non-overlapping intervals with max values
    """
    # Ensure sorted input
    df = df.sort_values(by=["chrom", "start"]).reset_index(drop=True)

    result_rows = []

    # Process each chromosome separately
    for chrom in df["chrom"].unique():
        sub_df = df[df["chrom"] == chrom]
        
        # Step 1: Get unique breakpoints
        breaks = sorted(set(sub_df["start"]).union(sub_df["end"]))

        # Step 2: Loop through consecutive pairs
        for i in range(len(breaks) - 1):
            seg_start = breaks[i]
            seg_end = breaks[i + 1]

            # Step 3: Find overlapping intervals
            overlaps = sub_df[(sub_df["start"] < seg_end) & (sub_df["end"] > seg_start)]
            if not overlaps.empty:
                max_val = overlaps["value"].max()
                result_rows.append([chrom, seg_start, seg_end, max_val])

    return pd.DataFrame(result_rows, columns=["chrom", "start", "end", "value"])

chrom = "mLynRuf2.2_ChrB2_rc"
start = 29_600_000
end = 32_250_000

# load data from predictions table
table1 = pd.read_csv('data/introgression_scans/lpa-wel_predictions_withM/predictions_table.tsv', sep='\t')
# rename columns: p_ab -> welab, p_ba -> welba, p_bi -> welbi and remove p_none
table1 = table1.rename(columns={"p_ab": "welab", "p_ba": "welba", "p_bi": "welbi"})
table1 = table1.drop(columns=["p_none"])

table2 = pd.read_csv('data/introgression_scans/lpa-sel_predictions_withM/predictions_table.tsv', sep='\t')
# rename columns: p_ab -> selab, p_ba -> selba, p_bi -> selbi and remove p_none
table2 = table2.rename(columns={"p_ab": "selab", "p_ba": "selba", "p_bi": "selbi"})
table2 = table2.drop(columns=["p_none"])

# load introgressed windows
intro = pd.read_csv('data/introgression_scans/bed_files/wel_and_sel_to_lpa_intro.merged.bed', sep='\t', names=["chrom", "start", "end"])

# mhc genes 
ci_genes = pd.read_csv('data/introgression_scans/genes/mhc_ci_genes.bed', sep='\t', names=["chrom", "start", "end", "gene"])
cii_genes = pd.read_csv('data/introgression_scans/genes/mhc_cii_genes.bed', sep='\t', names=["chrom", "start", "end", "gene"])

# all genes
all_genes = pd.read_csv('data/genes.bed', sep='\t', names=["chrom", "start", "end", "gene"])

# filter table to only include the region of interest
table1 = table1[(table1["chrom"] == chrom) & (table1["start"] >= start) & (table1["end"] <= end)]
table2 = table2[(table2["chrom"] == chrom) & (table2["start"] >= start) & (table2["end"] <= end)]
intro = intro[(intro["chrom"] == chrom) & (intro["start"] >= start) & (intro["end"] <= end)]
ci_genes = ci_genes[(ci_genes["chrom"] == chrom) & (ci_genes["start"] >= start) & (ci_genes["end"] <= end)]
cii_genes = cii_genes[(cii_genes["chrom"] == chrom) & (cii_genes["start"] >= start) & (cii_genes["end"] <= end)]
all_genes = all_genes[(all_genes["chrom"] == chrom) & (all_genes["start"] >= start) & (all_genes["end"] <= end)]

# Data
#start1 = np.array(table1["start"], dtype=float)
#end1 = np.array(table1["end"], dtype=float)
#wel_stairs = np.append(start1, end1[-1], axis=None)
#start2 = np.array(table2["start"], dtype=float)
#end2 = np.array(table2["end"], dtype=float)
#sel_stairs = np.append(start2, end2[-1])
#
#welab = np.array(table1["welab"], dtype=float) + np.array(table1["welbi"], dtype=float)
#selab = np.array(table2["selab"], dtype=float) + np.array(table2["selbi"], dtype=float)
#welba = np.array(table1["welba"], dtype=float) + np.array(table1["welbi"], dtype=float)
#selba = np.array(table2["selba"], dtype=float) + np.array(table2["selbi"], dtype=float)

# merge overlapping windows by selecting the maximum value
welab_df = pd.DataFrame(columns=["chrom", "start", "end", "value"])
welab_df["chrom"] = table1["chrom"]
welab_df["start"] = table1["start"]
welab_df["end"] = table1["end"]
welab_df["value"] = table1["welab"] + table1["welbi"]
selab_df = pd.DataFrame(columns=["chrom", "start", "end", "value"])
selab_df["chrom"] = table2["chrom"]
selab_df["start"] = table2["start"]
selab_df["end"] = table2["end"]
selab_df["value"] = table2["selab"] + table2["selbi"]
# merge overlapping windows by selecting the maximum value
welab_df = segment_max_intervals(welab_df)
selab_df = segment_max_intervals(selab_df)

# get data for plotting
start1 = np.array(welab_df["start"], dtype=float)
end1 = np.array(welab_df["end"], dtype=float)
wel_stairs = np.append(start1, end1[-1], axis=None)
start2 = np.array(selab_df["start"], dtype=float)
end2 = np.array(selab_df["end"], dtype=float)
sel_stairs = np.append(start2, end2[-1], axis=None)
welab = np.array(welab_df["value"], dtype=float)
selab = np.array(selab_df["value"], dtype=float)

# Stacked area plot for welab and selab
fig, ax = plt.subplots(figsize=(9, 3))

ax.broken_barh([z for z in zip(intro["start"], intro["end"] - intro["start"])], (0., 1.2),
                color='none', alpha=1, label='introgression', linewidth=1.5, edgecolor='black', linestyle='--')
#ax.stackplot((start1+end1)/2, welab, labels=['p_intro_wel'], baseline='zero', alpha=0.5, color='#003c82')
#ax.stackplot((start2+end2)/2, selab, labels=['p_intro_sel'], baseline='zero', alpha=0.5, color='#1b9773')

ax.stairs(welab, wel_stairs, alpha=0.5, label='p intro from ELw', baseline=0, fill=True, facecolor='#003c82', edgecolor='#003c82')
ax.stairs(selab, sel_stairs, alpha=0.5, label='p intro from ELs', baseline=0, fill=True, facecolor='#1b9773', edgecolor='#1b9773')

ax.broken_barh([z for z in zip(all_genes["start"], all_genes["end"] - all_genes["start"])], (-0.25, -0.5),
                color='lightgrey', alpha=1, edgecolor='lightgrey', linewidth=1, label='all genes')

ax.broken_barh([z for z in zip(ci_genes["start"], ci_genes["end"] - ci_genes["start"])], (-0.25, -0.5),
                color='#e3aa00', alpha=1, edgecolor='#e3aa00', linewidth=1, label='MHC class I')

ax.broken_barh([z for z in zip(cii_genes["start"], cii_genes["end"] - cii_genes["start"])], (-0.25, -0.5),
                color='#de7e1e', alpha=1, edgecolor='#de7e1e', linewidth=1, label='MHC class II')

# move legend outside of the plot
ax.legend(loc='upper right', bbox_to_anchor=(1.2, 1), fontsize=8, frameon=False)
ax.set_yticks([0, .25, .5, .75, 1])
ax.set_xticks(np.arange(29_500_000, 32_500_000, 500_000))
ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{round(x / 1e6, 2)} M'))
ax.set_xlabel("Chromosome B2 position (bp)")
ax.set_ylabel("")
ax.text(-0.02, -0.5, 'Genes', va='center', ha='center', transform=ax.get_yaxis_transform(), rotation=90)
ax.text(-0.09, 0.5, 'Probability\nof introgression', va='center', ha='center', transform=ax.get_yaxis_transform(), rotation=90)

# plt.show()
plt.savefig('plots/introgression_scans/introgression_region.pdf', bbox_inches='tight')
plt.close()