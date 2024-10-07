import argparse
import numpy as np
from seriate import seriate
from scipy.spatial.distance import pdist, cdist
from scipy.optimize import linear_sum_assignment
import os

def parse_args():
    parser = argparse.ArgumentParser(description = "Write binary fastas from a npz file")
    parser.add_argument("--ifile", help = "The input npz file with the data")
    parser.add_argument("--odir", help = "The output directory to store the fastas")
    parser.add_argument("--pop_sizes", help = "The population sizes")
    parser.add_argument("--shape", help = "The shape of the data that was used to train the model")
    parser.add_argument("--step_size", help = "The step size for the sliding window")
    return parser.parse_args()

def load_npz(ifile):
    ifile = np.load(ifile)
    # NOTE: the pop1 = sech and pop2 = sim is hardcoded from the way the data was generated
    # by src/introgression_scans/phasedVcfToNpz.py
    pop1_x = ifile['sechMatrix'].T
    pop2_x = ifile['simMatrix'].T

    x = np.vstack((pop1_x, pop2_x))

    # destroy the perfect information regarding
    # which allele is the ancestral one
    for k in range(x.shape[1]):
        if np.sum(x[:,k]) > x.shape[0] / 2.:
            x[:,k] = 1 - x[:,k]
        elif np.sum(x[:,k]) == x.shape[0] / 2.:
            if np.random.choice([0, 1]) == 0:
                x[:,k] = 1 - x[:,k]
    return x

def seriate_x(x, metric = 'cosine'):
    Dx = pdist(x, metric = metric)
    Dx[np.where(np.isnan(Dx))] = 0.
    ix = seriate(Dx, timeout = 0)
    return x[ix], ix

args = parse_args()
npz_file = args.ifile
pop_sizes = args.pop_sizes
pop_sizes = tuple(list(map(int, pop_sizes.split(','))))
out_shape = args.shape
out_shape = tuple(list(map(int, out_shape.split(','))))
step_size = int(args.step_size)
sorting = 'seriate_match'
sorting_metric = 'cosine'
pop1_size = pop_sizes[0]
pop2_size = pop_sizes[1]
pop_size = out_shape[1]
odir = args.odir

# load the data
x = load_npz(npz_file)
# break x into blocks of snps with selected window size and step size
w_size = out_shape[-1]
blocks = []
for i in range(0, x.shape[1], step_size):
    if x[:, i:i+w_size].shape[1] == w_size:
        blocks.append(x[:, i:i+w_size])

# create new blocks by sorting the populations
for z, block in enumerate(blocks):
    print(f'writing fasta of window: {z} in {odir}', end = '\r')
    # get the indices of the first population and upsample if needed
    x1_indices = list(range(pop1_size))
    n = pop_size - pop1_size
    replace = True if n > pop1_size else False
    if n > 0:
        x1_indices = x1_indices + list(np.random.choice(range(pop1_size), n, replace = replace))
    x1 = block[x1_indices,:]
    # get the indices of the second population and upsample if needed
    x2_indices = list(range(pop1_size, pop1_size + pop2_size))
    n = pop_size - pop2_size
    replace = True if n > pop2_size else False
    if n > 0:
        x2_indices = x2_indices + list(np.random.choice(range(pop1_size, pop1_size + pop2_size), n, replace = replace))
    x2 = block[x2_indices,:]
    # sort the populations by population 1 (x1)
    x1, ix1 = seriate_x(x1, metric = sorting_metric)
    D = cdist(x1, x2, metric = sorting_metric)
    D[np.where(np.isnan(D))] = 0.
    i, j = linear_sum_assignment(D)
    x2 = x2[j,:]
    # create the new block
    new_block = np.stack((x1, x2), axis=0)
    fasta = f'{os.path.abspath(odir)}/window_{z}.fa'
    for pop, seq in enumerate(new_block):
        for i in range(seq.shape[0]):
            with open(fasta, 'a') as file:
                file.write(f'>pop{pop+1}_{i}\n')
                # change 0 for "A" and 1 for "T"
                ats = ['A' if x == 0 else 'T' for x in seq[i]]
                file.write(''.join(ats) + '\n')

print(f'writing fasta of window: {z} in {odir} done')
