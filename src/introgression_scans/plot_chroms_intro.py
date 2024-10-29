import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import argparse

def parse_args():
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="Plot the introgressed regions in the chromosomes from a bed file of introgressed regions.")
    parser.add_argument("--ibed", type=str, help="Input bed file of introgressed regions")
    parser.add_argument("--ofile", type=str, help="Output file to save the plot")
    args = parser.parse_args()
    return args

def colors():
    """
    Define the colors for the introgressed regions and the non-introgressed regions.
    """
    return {
        'none': '#003c82',
        'intro': '#e3aa00'
    }

def plot_title(ibed):
    """
    Get the title for the plot.
    """
    if "wel_and_sel_to_lpa_intro" in ibed:
        return "Introgression from Eurasian to Iberian lynx"
    elif "lpa_to_wel_intro" in ibed:
        return "Introgression from Iberian to Western Eurasian lynx"
    elif "lpa_to_sel_intro" in ibed:
        return "Introgression from Iberian to Southern Eurasian lynx"

def plot_chroms(ibed, colors, ofile):
    """
    Plot the introgressed regions in the chromosomes from a bed file of introgressed regions.
    """
    bed = pd.read_csv(ibed, sep='\t', names=["chrom", "start", "end"])
    chrs = bed.chrom.unique()
    fig, ax = plt.subplots(figsize=(6*1.5, 8*1.5))
    for n, chrom in enumerate(chrs[::-1]):
        wins = bed[bed.chrom == chrom]
        nin = np.array([wins['end'] - wins['start']]).sum()
        ntot = max(wins["end"])
        pc = f'{round(nin/ntot*100, 2)}'
        ax.text(ntot, n*2 + 0.5, f"  {pc}%", verticalalignment='center', horizontalalignment='left', color='black')
        ax.broken_barh([(0, max(wins["end"]))], (n * 2, 1),
                        color=colors['none'], alpha=1, label='no introgression')
        ax.broken_barh([z for z in zip(wins["start"], wins["end"] - wins["start"])], (n * 2, 1), 
                        color=colors['intro'], alpha=1, label="introgression", linewidth=0.2)
        ax.broken_barh([(0, max(wins["end"]))], (n * 2, 1),
                        color='none', alpha=1, edgecolor='black', linewidth=1)
        ax.set_yticks([n*2 + 0.5 for n in range(len(chrs))])
    ax.set_yticklabels([f"{chrom.split('_')[1][3:5]}" for chrom in chrs[::-1]])
    ax.set_xlim(0, max(bed["end"] * 1.1))
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f'{round(x / 1e6, 2)} M'))
    ax.set_xlabel("Chromosome position (bp)")
    ax.set_ylabel("Chromosome")
    ax.set_title(plot_title(ibed))
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc='lower right')
    plt.savefig(ofile, bbox_inches='tight')

def main():
    args = parse_args()
    plot_chroms(args.ibed, colors(), args.ofile)

if __name__ == "__main__":
    main()
