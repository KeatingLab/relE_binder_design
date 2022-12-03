from inspect import ArgInfo
import os, sys, glob, argparse, json
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

parser = argparse.ArgumentParser(description='Train three-layer perceptron')
parser.add_argument('--name', type=str, required=True,
                    help='A name specifying the network that is being trained')
parser.add_argument('--trainData', type=str, required=True,
                    help='The path to the training data table')
parser.add_argument('--valData', type=str, required=True,
                    help='The path to the validation data table')
parser.add_argument('--testData', type=str, required=True,
                    help='The path to the test data table')               
parser.add_argument('--batchSize', type=int,
                    help='', default = 512)
parser.add_argument('--epochs', type=int,
                    help='', default = 200)
parser.add_argument('--learningRate', type=float,
                    help='The learning rate that is used to compute the parameter updates during training', default = 1e-3)
parser.add_argument('--pathToRepo', type=str,
                    help='The path to the python directory of the interfaceGenerator repo on the local system')
parser.add_argument('--earlyStopping', type=int,
                    help='If provided, will stop training if the validation loss stops decreasing. The value provided is the "patience", i.e. the number of epochs to continue training while the validation loss is not improving before giving up. The final model weights are always from whichever epoch yielded the best validation error',
                    default=-1)
parser.add_argument('--json', type=str,
                    help='A .json file specifying the network properties',
                    default="")
args = parser.parse_args()

# Library code
sys.path.append(args.pathToRepo)
from datasets.data_utils import *
from datasets.data_loading import *
from model.network import *
from train_model.training import *
from train_model.train_utils import *

if __name__ == "__main__":
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"Using {device} device")

    networkConfig = {}
    if (args.json != ""):
        with open(args.json) as file:
            networkConfig = json.load(file)
    print("network configs: ",networkConfig)
    model = contactScorer(networkConfig,device).to(device)
    print(model)

    shuffle=True

    # Initialize the optimizer
    optimizer = torch.optim.SGD(model.parameters(), lr=args.learningRate, momentum = 0.9)

    # Initialize the early_stopping object
    patience = args.earlyStopping if args.earlyStopping != -1 else args.epochs
    early_stopping = EarlyStopping(patience=10, delta=0, verbose=True)

    # Load the data
    train_data = ContactDataset(args.trainData)
    val_data = ContactDataset(args.valData)
    test_data = ContactDataset(args.testData)

    # Create data loaders
    train_dataloader = DataLoader(train_data, batch_size=args.batchSize, shuffle=shuffle)
    val_dataloader = DataLoader(val_data, batch_size=args.batchSize)
    test_dataloader = DataLoader(test_data, batch_size=args.batchSize)

    if (networkConfig['output']=='bindingScore'):
        # Initialize the loss function
        loss_fn = nn.MSELoss()
    elif (networkConfig['model']=='predictAADist'):
        # Initialize the loss function
        loss_fn = nn.CrossEntropyLoss() #expects unnormalized scores

    # Train the network
    trained_model, loss_df = train_model(args.epochs, train_dataloader, val_dataloader, model, loss_fn, optimizer, early_stopping, args.name)

    # Test the network after training
    test_loss = test(test_dataloader, model, loss_fn)
    print(f"Final test loss: {test_loss:>7f}")

    # Write the final weights out
    model_scripted = torch.jit.script(model)
    model_scripted.save(args.name+'.pt')

    print("Done!")