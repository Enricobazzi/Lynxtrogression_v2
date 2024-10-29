import os
import gzip
import numpy as np

def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description='Convert msOut file to npz file')
    parser.add_argument('--idir', type=str, help='input directory with msOut file')
    parser.add_argument('--ofile', type=str, help='output npz file')
    parser.add_argument('--pop_sizes', type=str, help='comma separated list of population sizes')
    return parser.parse_args()

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

def read_msfile(ms_file):
    """
    Read ms file and return lines
    """
    with gzip.open(ms_file, 'rt') as f:
        lines = f.readlines()
    return lines

def get_matrices(blocks, pop_sizes):
    """
    Get matrices from ms blocks
    """
    gts = {}
    for p in range(sum(pop_sizes)):
        if p < pop_sizes[0]:
            gts[f'pop1_{p}'] = []
        else:
            gts[f'pop2_{p}'] = []
    for block in blocks:
        for p in range(sum(pop_sizes)):
            if p < pop_sizes[0]:
                for gt in block[p+3].strip():
                    gts[f'pop1_{p}'].append(int(gt))
            else:
                for gt in block[p+3].strip():
                    gts[f'pop2_{p}'].append(int(gt))
    sechMatrix = np.stack([np.array(gts[f'pop1_{p}']) for p in range(pop_sizes[0])], axis=0).T
    simMatrix = np.stack([np.array(gts[f'pop2_{p}']) for p in range(pop_sizes[0], sum(pop_sizes))], axis=0).T
    return sechMatrix, simMatrix

def main():
    args = parse_args()
    idir = args.idir
    ms_file = f'{os.path.abspath(idir)}/mig.msOut.gz'
    ofile = args.ofile
    pop_sizes = args.pop_sizes
    pop_sizes = tuple(list(map(int, pop_sizes.split(','))))
    blocks = get_ms_blocks(read_msfile(ms_file))
    sechMatrix, simMatrix = get_matrices(blocks, pop_sizes)
    np.savez_compressed(ofile, sechMatrix = sechMatrix, simMatrix = simMatrix)

if __name__ == '__main__':
    main()