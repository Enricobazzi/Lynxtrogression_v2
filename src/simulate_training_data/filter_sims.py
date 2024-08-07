# dummy arguments
import gzip
import argparse
import os

argparser = argparse.ArgumentParser('filter msOut and anc files')
argparser.add_argument('--idir', type=str, help='input directory with msOut and anc files')
argparser.add_argument('--odir', type=str, help='output directory for filtered msOut and anc files')
argparser.add_argument('--n_sites', type=int, help='minimum number of segregating sites')
argparser.add_argument('--n_sims', type=int, help='total number of simulations to keep')
args = argparser.parse_args()

idir = args.idir
odir = args.odir
n_sites = args.n_sites
n_sims = args.n_sims

print(f'Filtering msOut and anc files in {idir}')
print(f'Keeping simulations with at least {n_sites} segregating sites')
print(f'Keeping {n_sims} simulations')
print(f'Writing filtered msOut and anc files in {odir}')

def get_msout_lines(msout):
    """
    Extract lines from msOut file
    """
    with gzip.open(msout, 'rt') as f:
        lines = f.readlines()
    return lines

def get_segsites(lines):
    """
    Extract list of segsites per sim from msOut file
    """
    segsites = []
    for line in lines:
        if 'segsites' in line:
            segsites.append(int(line.strip().split()[1]))
    return segsites

def get_first_three(lines):
    """
    Extract ms command from msOut file
    """
    return list(lines[:3])

def get_last_line(lines):
    """
    Extract last line from msOut file
    """
    return lines[-1]


def get_ms_blocks(lines):
    """
    Split msout file in blocks separated by lines starting with '//'
    """
    # remove first three lines (ms command)
    mslines = lines.copy()
    mslines.pop(0)
    mslines.pop(0)
    mslines.pop(0)
    # remove last line (admixture proportion from msmodified)
    mslines.pop()
    # create blocks
    blocks = []
    current_block = []
    for line in mslines:
        if line.startswith('//') and current_block:
            blocks.append(current_block)
            current_block = [line]
        else:
            current_block.append(line)
    if current_block:
        blocks.append(current_block)
    return blocks

def get_anc_blocks(anc, segsites):
    """
    Split anc file in blocks with lines equal to the number of segsites
    """
    with gzip.open(anc, 'rt') as f:
        lines = f.readlines()
    anc_blocks = []
    current_line = 0
    for i, n in enumerate(segsites):
        start_index = current_line
        end_index = current_line + n
        block_lines = lines[start_index:end_index]
        anc_blocks.append(block_lines)
        current_line = end_index
    return anc_blocks

# get the files
msout = f'{idir}/mig.msOut.gz'
anc = f'{idir}/out.anc.gz'

# get lines of msOut file
lines = get_msout_lines(msout)
# parse the lines extracting segsites, ms simulation blocks, first three lines and last line
segsites = get_segsites(lines)
ms_blocks = get_ms_blocks(lines)
ms_first_three = get_first_three(lines)
ms_last = get_last_line(lines)
# parse the anc file in blocks
anc_blocks = get_anc_blocks(anc, segsites)

# keep ms_blocks and anc_blocks with more than n_sites
# until n_sims is reached
ms_blocks_filtered = []
anc_blocks_filtered = []
for s, block, anc_block in zip(segsites, ms_blocks, anc_blocks):
    if s >= n_sites:
        ms_blocks_filtered.append(block)
        anc_blocks_filtered.append(anc_block)
    if len(ms_blocks_filtered) == n_sims:
        break

# write filtered msOut and anc files
with open(f'{odir}/mig.msOut', 'w') as f:
    f.write(''.join(ms_first_three))
    for block in ms_blocks_filtered:
        f.write(''.join(block))
    f.write(ms_last)
with open(f'{odir}/out.anc', 'w') as f:
    for block in anc_blocks_filtered:
        f.write(''.join(block))

# gzip the files
os.system(f'gzip {odir}/mig.msOut')
os.system(f'gzip {odir}/out.anc')
