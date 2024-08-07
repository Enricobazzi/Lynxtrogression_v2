# -*- coding: utf-8 -*-
import os
import argparse
import logging

import torch
from torch.nn import CrossEntropyLoss, KLDivLoss, NLLLoss, SmoothL1Loss, BCELoss, BCEWithLogitsLoss
import torch

from torch import nn

import torch.optim as optim
from torch.optim.lr_scheduler import ReduceLROnPlateau
import torch.utils.data
import torch.utils.data.distributed
import torch.nn.functional as F
import torch.distributed as dist
import copy

from layers import LexStyleNet
from data_loaders import H5DisDataGenerator
import h5py

import numpy as np
from sklearn.metrics import accuracy_score
import pandas as pd

import glob
from torchvision_mod_layers import resnet18, resnet34, resnet50, resnext50_32x4d
from sklearn.metrics import accuracy_score, confusion_matrix
import pandas as pd

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns

idir = 'data/introgression_scans/simulations'
ifiles = dict()
    
_ = glob.glob(os.path.join(idir, '*.hdf5'))
    
for f in _:
    t = f.split('/')[-1].split('.')[0]
    
    ifiles[t] = [f]
    
print(ifiles)
generator = H5DisDataGenerator(ifiles, batch_size = int(16))
print(generator.classes)
print(sorted(generator.classes))
x, y = generator.get_batch()
