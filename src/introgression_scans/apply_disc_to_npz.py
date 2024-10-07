import numpy as np
import argparse
from seriate import seriate
from scipy.spatial.distance import pdist, cdist
from scipy.optimize import linear_sum_assignment
from torchvision_mod_layers import resnet34
import torch
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(description = "Apply the discriminative embedding to a npz file")
    parser.add_argument("--ifile", help = "The input npz file with the data")
    parser.add_argument("--ofile", help = "The output csv file with the predictions")
    parser.add_argument("--weights", help = "The weights file")
    parser.add_argument("--pop_sizes", help = "The population sizes")
    parser.add_argument("--shape", help = "The shape of the data that was used to train the model")
    parser.add_argument("--step_size", help = "The step size for the sliding window")
    parser.add_argument("--in_channels", help = "The number of populations")
    parser.add_argument("--n_classes", help = "The number of classes in the model")
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

def softmax(logits):
    exp_logits = np.exp(logits - np.max(logits, axis=-1, keepdims=True))
    probabilities = exp_logits / np.sum(exp_logits, axis=-1, keepdims=True)
    return probabilities

# parse the arguments
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
weights = args.weights
ofile = args.ofile
in_channels = args.in_channels
n_classes = args.n_classes

### load and preprocess the data ###

x = load_npz(npz_file)
# break x into blocks of snps with selected window size and step size
w_size = out_shape[-1]
blocks = []
for i in range(0, x.shape[1], step_size):
    if x[:, i:i+w_size].shape[1] == w_size:
        blocks.append(x[:, i:i+w_size])
# create new blocks by sorting the populations
new_blocks = []
for block in blocks:
    print(f'processing block {n}', end = '\r')
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
    new_blocks.append(np.stack((x1, x2), axis=0))

### load the model ###

# torch device
device = torch.device('cpu')
# load model and weights
model = resnet34(in_channels = int(in_channels), num_classes = int(n_classes))
model = model.to(device)
checkpoint = torch.load(weights, map_location = device)
model.load_state_dict(checkpoint)
model.eval()

### apply the model to the data ###

# initialize the arrays to store the probabilities
wins = np.arange(0, len(new_blocks), dtype=int)
p_ab = np.zeros((len(new_blocks),), dtype=float)
p_ba = np.zeros((len(new_blocks),), dtype=float)
p_bi = np.zeros((len(new_blocks),), dtype=float)
p_none = np.zeros((len(new_blocks),), dtype=float)
print('generating predictions...')
# apply the model to the data
for i in range(len(new_blocks)):
    # expand as it's expecting batched data (batch size of 1)
    X = np.expand_dims(np.array(new_blocks[i]), axis=0)
    # get predictions
    with torch.no_grad():
        x = torch.FloatTensor(X).to(device)
        y_pred = model(x)
    y_pred = y_pred.detach().cpu().numpy()
    # softmax the predictions
    y_pred = softmax(y_pred.flatten())
    # get the probabilities
    p_ab[i] = y_pred[0]
    p_ba[i] = y_pred[1]
    p_bi[i] = y_pred[2]
    p_none[i] = y_pred[3]

# create df and save
df = pd.DataFrame({'window': wins, 'p_ab': p_ab, 'p_ba': p_ba, 'p_bi': p_bi, 'p_none': p_none})
df.to_csv(ofile, index = False)
