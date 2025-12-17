import gzip
import argparse
import os
import numpy as np

def parse_args():
    argparser = argparse.ArgumentParser('add genotyping error to simulations based on average depth and number of low-depth individuals per population')
    argparser.add_argument('--idir', type=str, help='input directory with msOut and anc files')
    argparser.add_argument('--odir', type=str, help='output directory for msOut and anc files with genotyping error')
    argparser.add_argument('--pop_sizes', type=str, help='comma separated list of population sizes')
    argparser.add_argument('--nbad', type=str, help='comma separated list of number of low-depth individuals per population')
    argparser.add_argument('--avg_depth', type=float, help='average sequencing depth for low-depth individuals')
    return argparser.parse_args()

def get_msout_lines(msout):
    """
    Extract lines from msOut file
    """
    with gzip.open(msout, 'rt') as f:
        lines = f.readlines()
    return lines

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

def get_matrix_from_block(ms_block):
    """
    Create genotype matrix from ms block
    """
    block_lines = ms_block[3:-1].copy()
    matrix = np.zeros((len(block_lines), len(block_lines[0].strip())), dtype=int)
    for i, line in enumerate(block_lines):
        matrix[i, :] = np.array(list(line.strip()), dtype=int)
    return matrix

def get_gtmatrix_from_matrix(matrix):
    """
    Create genotype matrix from haplotype matrix
    """
    gtmatrix = np.zeros((matrix.shape[0]//2, matrix.shape[1]), dtype=int)
    for i in range(matrix.shape[0]//2):
        gtmatrix[i, :] = matrix[2*i, :] + matrix[2*i+1, :]
    return gtmatrix

def draw_depth(avg_depth):
    """
    Draw depth from poisson distribution
    """
    depth = np.random.poisson(avg_depth)
    return depth if depth > 0 else 1

def draw_bases(depth):
    """
    Draw bases for a site given depth
    """
    # probability of sampling ref is 0.5
    p_ref = 0.5
    n_ref = np.random.binomial(depth, p_ref)
    n_alt = depth - n_ref
    return n_ref, n_alt

def call_gt(n_ref, n_alt):
    """
    Call genotype based on number of ref and alt bases
    """
    if n_ref == 0:
        return 2
    elif n_alt == 0:
        return 0
    else:
        return 1
 
def get_newgts_for_ind(ind_gtmatrix, avg_depth):
    """
    Get new genotypes for an individual given average depth
    """
    new_gts = np.zeros(ind_gtmatrix.shape, dtype=int)
    for i, gt in enumerate(ind_gtmatrix):
        if gt == 0 or gt == 2:
            new_gts[i] = gt
        else:
            depth = draw_depth(avg_depth)
            n_ref, n_alt = draw_bases(depth)
            new_gts[i] = call_gt(n_ref, n_alt)
    return new_gts

def build_new_gtmatrix(gtmatrix, nbad, avg_depth, pop_sizes):
    """
    Build new genotype matrix with genotyping errors for nbad individuals per population
    without resampling same individual
    """
    new_gtmatrix = gtmatrix.copy()
    # first population
    ind_indexes = np.random.choice(
        pop_sizes[0] // 2,
        size=nbad[0],
        replace=False
    )
    for i, ind_index in enumerate(ind_indexes):
        new_gts = get_newgts_for_ind(gtmatrix[ind_index, :], avg_depth)
        new_gtmatrix[ind_index, :] = new_gts
        #print(f'Pop1 bad ind {i + 1}: index {ind_index}')
    # second population
    ind_indexes = np.random.choice(
        pop_sizes[1] // 2,
        size=nbad[1],
        replace=False
    )
    for i, ind_index in enumerate(ind_indexes):
        new_gts = get_newgts_for_ind(gtmatrix[ind_index + pop_sizes[0] // 2, :], avg_depth)
        new_gtmatrix[ind_index + pop_sizes[0] // 2, :] = new_gts
        #print(f'Pop2 bad ind {i + 1}: index {ind_index + pop_sizes[0] // 2}')
    return new_gtmatrix

def get_new_matrix(matrix, gtmatrix, new_gtmatrix):
    """
    Get new haplotype matrix based on new genotype matrix
    """
    new_matrix = matrix.copy()
    diffs = np.where(gtmatrix != new_gtmatrix)
    # modify original matrix (haplotypes) based on new genotype matrix
    for row, col in zip(diffs[0], diffs[1]):
        new_gt = new_gtmatrix[row, col]
        # if new gt is 0, change the haplotype that was 1 to 0
        if new_gt == 0:
            hap_indices = [2*row, 2*row+1]
            for hi in hap_indices:
                if new_matrix[hi, col] == 1:
                    new_matrix[hi, col] = 0
                    break
        # if new gt is 2, change the haplotype that was 0 to 1
        elif new_gt == 2:
            hap_indices = [2*row, 2*row+1]
            for hi in hap_indices:
                if new_matrix[hi, col] == 0:
                    new_matrix[hi, col] = 1
                    break
    return new_matrix

def build_new_ms_block(ms_block, new_matrix):
    """
    Build new ms block from new haplotype matrix
    """
    new_ms_block = [ms_block[0], ms_block[1], ms_block[2]]
    for row in range(new_matrix.shape[0]):
        new_line = ''.join(new_matrix[row, :].astype(str)) + '\n'
        new_ms_block.append(new_line)
    new_ms_block.append(ms_block[-1])
    return new_ms_block

def main():
    args = parse_args()
    idir = args.idir
    odir = args.odir
    pop_sizes = args.pop_sizes
    pop_sizes = tuple(list(map(int, pop_sizes.split(','))))
    nbad = args.nbad
    nbad = tuple(list(map(int, nbad.split(','))))
    avg_depth = args.avg_depth
    msout = f'{idir}/mig.msOut.gz'
    # read msOut file
    lines = get_msout_lines(msout)
    ms_blocks = get_ms_blocks(lines)
    ms_first_three = get_first_three(lines)
    ms_last = get_last_line(lines)
    
    # test on first ms block
    # write filtered msOut and anc files
    with open(f'{odir}/mig.msOut', 'w') as f:
        f.write(''.join(ms_first_three))
        for ms_block in ms_blocks:
            matrix = get_matrix_from_block(ms_block)
            gtmatrix = get_gtmatrix_from_matrix(matrix)
            new_gtmatrix = build_new_gtmatrix(gtmatrix, nbad, avg_depth, pop_sizes)
            new_matrix = get_new_matrix(matrix, gtmatrix, new_gtmatrix)
            diffs = np.sum(matrix != new_matrix)
            print(f'Number of genotype changes in this block: {diffs}, out of {matrix.shape[0] * matrix.shape[1]} genotypes ({(diffs / (matrix.shape[0] * matrix.shape[1])) * 100:.2f}%)')
            new_ms_block = build_new_ms_block(ms_block, new_matrix)
            f.write(''.join(new_ms_block))
        f.write(ms_last)
    # gzip the file
    os.system(f'gzip {odir}/mig.msOut')
    # copy anc file without changes
    os.system(f'cp {idir}/out.anc.gz {odir}/out.anc.gz')

if __name__ == '__main__':
    main()
