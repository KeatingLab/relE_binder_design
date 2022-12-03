import glob
import math
from enum import Enum
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import torch
from torch.utils.data import Dataset

class ContactDataset(Dataset):
    def __init__(self, data_file):
        self.structure_data = pd.read_csv(data_file)

    def __len__(self):
        return len(self.structure_data)

    def __getitem__(self, idx):
        ### input values
        # one-hot encoding of target aa
        aaIdx = self.structure_data.loc[idx,'targetResAAIdx']
        aaEncoding = torch.tensor(np.array([aaIdx]),dtype=torch.int64)
        
        # 16 distance values between backbone atoms of two residues
        distances = torch.unsqueeze(torch.tensor(self.structure_data.loc[idx,['N-N','N-Ca','N-C','N-O',
                                                             'Ca-N','Ca-Ca','Ca-C','Ca-O',
                                                             'C-N','C-Ca','C-C','C-O',
                                                             'O-N','O-Ca','O-C','O-O']],
                                                              dtype=torch.float),dim=0)
        
        ### output values
        # the interface contact score
        score = torch.unsqueeze(torch.tensor(self.structure_data.loc[idx,'bindingScore'],dtype=torch.float),dim=0)
        aaDist = torch.unsqueeze(torch.nn.functional.softmax(torch.tensor(self.structure_data.loc[idx,['A','C','D','E','F',
                                                                                        'G','H','I','K','L',
                                                                                        'M','N','P','Q','R',
                                                                                        'S','T','V','W','Y']],
                                                                                        dtype=torch.float)+1,
                                                                                        dim=0),dim=0)
        return ({'aa':aaEncoding,
                'distances':distances},
                {'bindingScore':score,
                 'aminoAcidProbability':aaDist})