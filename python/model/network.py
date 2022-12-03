import os, glob, argparse
from typing import Dict
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

from datasets.data_utils import *

class contactScorer(nn.Module):
    # device = "cuda" if torch.cuda.is_available() else "cpu"
    """ Network for computing the interface score """
    def __init__(self, config, device):
        super().__init__()

        # Set general params
        self.device = device

        # Set network parameters
        self.min_distance = config['min_distance'] if 'min_distance' in config else 2.5 #distances in angstroms
        self.max_distance = config['max_distance'] if 'max_distance' in config else 13.5
        self.n_rbf = config['encoding']['nRBF'] if 'encoding' in config else 12
        self.nNodes = config['model']['nNodes'] if 'model' in config else 64
        self.nLinLayer = config['model']['nLinLayer'] if 'model' in config else 1
        self.output = config['output'] if 'output' in config else 'bindingScore'
        self.config = config

        # Define RBF encoding tensors
        self.register_buffer('mu',torch.unsqueeze(torch.linspace(self.min_distance,self.max_distance,self.n_rbf,dtype=torch.float),dim=1))
        self.register_buffer('sigma2',torch.tensor([(self.max_distance-self.min_distance)/self.n_rbf],dtype=torch.float))

        # Set background amino acid probabilities
        self.register_buffer('bg_aa_prob',torch.unsqueeze(torch.tensor([prob for aa3,prob in aa3tosurfProb.items()],dtype=torch.float),dim=0))

        # Define network
        self.flatten = nn.Flatten()
        # self.network_stack = nn.Sequential(
        #     nn.Linear(16*self.n_rbf,self.nNodes),
        #     nn.LogSigmoid(),
        #     nn.Linear(self.nNodes, 20)
        # )
        self.layers = nn.ModuleList([
            nn.Linear(16*self.n_rbf,self.nNodes),
            nn.LogSigmoid(),
            nn.Linear(self.nNodes, 20)])

        # layers = [nn.Linear(16*self.n_rbf,self.nNodes),nn.LogSigmoid()]
        # for x in range(1,self.nLinLayer):
        #     layers += [nn.Linear(self.nNodes, self.nNodes),nn.LogSigmoid()]
        # layers += [nn.Linear(self.nNodes, 20)]
        # self.network_stack = nn.Sequential(*layers)

    def rbfEncoding(self, X_distances):
        """ Encode backbone atom distances using a series of gaussian radial basis functions
        Args
        ----
        X_distances : torch.tensor.float
            The distances between N, Ca, C, and O atoms of two residues
            Shape: n_batches x 16 x 1
        Returns
        -------
        X_encoding : torch.tensor.float
            The final tensor after flattening along the 1st dimension
            Shape: n_batches x 16 * self.n_rbf
        """
        return self.flatten(torch.exp(-torch.square((X_distances-self.mu)/self.sigma2)))
        
    def forward(self, X: Dict[str, torch.Tensor]) -> torch.Tensor:
        """ Predict the interface score of the residue pair
        Args
        ----
        X : dictionary{X_distances,X_aa}
        X_distances: torch.tensor.float
            Shape: n_batches x 16 x 1
        X_aa: torch.tensor.int8
            Shape: n_batches x 1
        Returns
        -------
        X_encoding : torch.tensor.float
            The interface score for each residue pair in the batch
            Shape: n_batches x 1
        """
        X_distances = X['distances']
        print('X_distances',X_distances.shape,torch.isnan(X_distances).any(),X_distances)
        n_batches = X_distances.shape[0]
        if self.output == 'bindingScore':
            X_aa = X['aa']
            # print('X_aa',X_aa.shape,X_aa[0])

            # Compute numerator
            X_encoding = self.rbfEncoding(X_distances)
            print('X_encoding',X_encoding.shape,torch.isnan(X_encoding).any(),X_encoding[0])
            for layer in self.layers:
                X_encoding = layer(X_encoding)
                print('X_encoding',X_encoding.shape,torch.isnan(X_encoding).any(),X_encoding[0])
            aaProbDist = X_encoding
            # aaProbDist = self.network_stack(X_encoding) # n_batches X 20
            print('aaProbDist',aaProbDist.shape,torch.isnan(aaProbDist).any(),aaProbDist[0])
            aaProbNum = torch.gather(aaProbDist,dim=1,index=X_aa)
            print('aaProbNum',aaProbNum.shape,torch.isnan(aaProbNum).any(),aaProbNum[0])

            # Compute denominator 
            aaProbDenom = torch.gather(self.bg_aa_prob.expand(n_batches,20),dim=1,index=X_aa)

            score = -torch.log(aaProbNum/aaProbDenom)
            # print(score)
            return score

        elif self.output == 'aminoAcidProbability':
            X_encoding = self.rbfEncoding(X_distances)
            aaProbDist = self.network_stack(X_encoding) # n_batches X 20
            return aaProbDist
            