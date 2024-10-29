# dummy arguments
import gzip
import argparse
import os
import numpy as np

argparser = argparse.ArgumentParser('filter msOut and anc files')
argparser.add_argument('--idir', type=str, help='input directory with msOut and anc files')
argparser.add_argument('--odir', type=str, help='output directory for filtered msOut and anc files')
argparser.add_argument('--n_sites', type=int, help='minimum number of segregating sites')
argparser.add_argument('--n_sims', type=int, help='total number of simulations to keep')
argparser.add_argument('--pop_sizes', type=str, help='comma separated list of population sizes')
argparser.add_argument('--migration', type=str, help='ab, ba, bi or none')
args = argparser.parse_args()

idir = args.idir
odir = args.odir
n_sites = args.n_sites
n_sims = args.n_sims
pop_sizes = args.pop_sizes
pop_sizes = tuple(list(map(int, pop_sizes.split(','))))
migration = args.migration

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

def get_ws_we(anc_block, n_sites, migration, pop_sizes):
    """
    Get the start and end of the window in the anc block. At least 1 site has to be introgressed in the window.
    Migration direction and population sizes are needed to determine introgression happening in the window.
    Try 20 times before giving up.
    """
    # get start positions of windows
    starts = [n for n in range(0, len(anc_block), n_sites)]
    # if the last window is not complete, change it to start before
    if len(anc_block) - starts[-1] < n_sites:
        starts[-1] = len(anc_block) - n_sites
    for _ in range(20):
        # choose a random window
        ws = np.random.choice(starts)
        we = ws + n_sites
        # check if there is introgression in the window
        if migration == 'ab':
            c = [line.strip()[pop_sizes[0]:].count('1') for line in anc_block[ws:we]]
            c = sum(c)
        elif migration == 'ba':
            c = [line.strip()[:pop_sizes[0]].count('1') for line in anc_block[ws:we]]
            c = sum(c)
        elif migration == 'bi':
            c = [line.strip().count('1') for line in anc_block[ws:we]]
            c = sum(c)
        elif migration == 'none':
            c = 1
        if c > 0:
            return ws, we
    return None, None

def filter_anc_block(anc_block, ws, we):
    """
    Filter the anc block to only include the window of interest
    """
    return anc_block[ws:we]

def filter_ms_block(ms_block, ws, we):
    """
    Filter the ms block to only include the window of interest
    """
    # modify the positions line of the ms block
    positions = []
    positions.append('positions:')
    for p in ms_block[2].split('  ')[ws+1:we+1]:
        positions.append(p)
    filtered_ms_block = []
    for n, line in enumerate(ms_block):
        if n == 0:
            filtered_ms_block.append(line)
        elif n == 1:
            filtered_ms_block.append(f'segsites: {n_sites}\n')
        elif n == 2:
            filtered_ms_block.append(f"{'  '.join(positions)}\n")
        else:
            filtered_ms_block.append(f'{line[ws:we]}\n')
    return filtered_ms_block

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
        ws, we = get_ws_we(anc_blocks[0], n_sites, migration, pop_sizes)
        if ws:
            block = filter_ms_block(block, ws, we)
            anc_block = filter_anc_block(anc_block, ws, we)
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
