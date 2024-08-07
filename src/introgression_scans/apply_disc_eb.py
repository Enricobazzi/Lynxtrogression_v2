import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import torch
import h5py
import sys
import pandas as pd
import glob
from torchvision_mod_layers import resnet34
from data_loaders import H5DisDataGenerator
from scipy.special import expit
import argparse

# read arguments
parser = argparse.ArgumentParser()
parser.add_argument('--weights', type = str, help= 'path to the weights of the model')
parser.add_argument('--ifile', type = str, help = 'path to the input hdf5 file')
parser.add_argument('--ofile', type = str, help = 'path to the output csv file with the predictions')
parser.add_argument('--in_channels', type = int, help = 'number of input channels', default=2)
parser.add_argument('--n_classes', type = int, help = 'number of classes', default=4)
args = parser.parse_args()

weights = args.weights
ifile = args.ifile
ofile = args.ofile
in_channels = args.in_channels
n_classes = args.n_classes

# torch device
device = torch.device('cpu')

# load model and weights
model = resnet34(in_channels = int(in_channels), num_classes = int(n_classes))
model = model.to(device)
checkpoint = torch.load(weights, map_location = device)
model.load_state_dict(checkpoint)
model.eval()

# # read data
# ifile = h5py.File(ifile, 'r')
# keys = list(ifile.keys())
#     
# shape = tuple(ifile['shape'])
# l = shape[-1]
#     
# Y = np.zeros((l, int(n_classes)))
# count = np.zeros((l, int(n_classes)))
# 
# # create df to store predictions
# data = []
# 
# for key in keys:
#     try:
#         int(key)
#     except:
#         continue
#         
#     X = np.array(ifile[key]['x_0'])
#     X = np.squeeze(X, axis=1)
#     pos = np.array(ifile[key]['positions'])
#     
#     # let's forward the whole batch through segm
#     with torch.no_grad():
#         x = torch.FloatTensor(X).to(device)
#         y_pred = model(x)
#     y_pred = y_pred.detach().cpu().numpy()
# 
#     for i in range(len(y_pred)):
#         start = pos[i][0]
#         end = pos[i][-1]
#         y_argmax = np.argmax(y_pred[i])
#         # exp_logits = np.exp(y_pred[i])
#         # probs = exp_logits / np.sum(exp_logits)
#         # store predictions in df
#         data.append({'start': start, 'end': end, 'p_ab': y_pred[0], 'p_ba': y_pred[1], 'p_bi': y_pred[2], 'p_none': y_pred[3], 'pred': y_argmax})
# 
# # create pandas df
# df = pd.DataFrame(data)
# # sort by start
# df = df.sort_values(by='start')
# # add window number column
# df['window'] = range(0, len(df))
# # save to csv
# df.to_csv(ofile, index=False)

ifile = h5py.File(args.ifile, 'r')
keys = list(ifile.keys())
    
shape = tuple(ifile['shape'])
l = shape[-1]
    
Y = np.zeros((l, int(args.n_classes)))
count = np.zeros((l, int(args.n_classes)))

for key in keys:
    try:
        int(key)
    except:
        continue
    
    X = np.array(ifile[key]['x_0'])
    X = np.squeeze(X, axis=1)
    pos = np.array(ifile[key]['positions'])
    indices_ = np.array(ifile[key]['pi'])
    
    # let's forward the whole batch through segmentation
    with torch.no_grad():
        x = torch.FloatTensor(X).to(device)
        y_pred = model(x)
    
    y_pred = y_pred.detach().cpu().numpy()
    
    for k in range(indices_.shape[0]):
        ip = indices_[k]
        
        Y[ip,:] += y_pred[k]
        count[ip,:] += 1
        
Y = expit(Y / count)

# save to csv
df = pd.DataFrame(Y)
df.to_csv(ofile, index=False)