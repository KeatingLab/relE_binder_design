import os
import glob
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

import torch
from torch import nn
from torch.utils.data import Dataset
from torch.utils.data import DataLoader

from datasets.data_utils import *
from datasets.data_loading import *
from train_model.train_utils import *

def train_epoch(dataloader, model, loss_fn, optimizer):
    model.train()
    size = len(dataloader.dataset)
    num_batches = len(dataloader)
    total_train_loss = 0.0
    for batch, (X, y) in enumerate(dataloader):
        for key in X.keys():
            X[key] = X[key].to(model.device)
        y = y[model.config['output']].to(model.device)

        # Compute prediction and loss
        pred = model(X)
        loss = loss_fn(pred, y) #mean loss over the batch

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if batch % (num_batches//10) == 0:
            loss, current = loss.item(), batch * len(y)
            print(f"Mean train loss: {loss:>7f}  [{current:>5d}/{size:>5d}]")
            total_train_loss += loss 
        else:
            total_train_loss += loss.item()
    
    return total_train_loss / num_batches

def test(dataloader, model, loss_fn):
    model.eval()
    total_test_loss = 0.0
    num_batches = len(dataloader)
    with torch.no_grad():
        for X, y in dataloader:
            for key in X.keys():
                X[key] = X[key].to(model.device)
            y = y[model.config['output']].to(model.device)

            pred = model(X)
            total_test_loss += loss_fn(pred, y).item() #mean loss over the batch
    
    return total_test_loss / num_batches

def train_model(n_epoch, train_dataloader, val_dataloader, model, loss_fn, optimizer, early_stopping, name):

    # Check the loss before any training
    val_loss = test(val_dataloader,model,loss_fn)
    print(f"Validation loss before training: {val_loss:>7f}")
    torch.save(model.state_dict(), "initial_weights.pt")

    mean_train_loss_per_epoch = list()
    mean_val_loss_per_epoch = list()

    for epoch in range(1,n_epoch + 1):
        # train step
        train_loss = train_epoch(train_dataloader, model, loss_fn, optimizer)
        mean_train_loss_per_epoch.append(train_loss)
        print(f"Mean training loss over epoch {epoch}: {mean_train_loss_per_epoch[-1]:>7f}")

        # validate step
        val_loss = test(val_dataloader,model,loss_fn)
        mean_val_loss_per_epoch.append(val_loss)
        print(f"Mean validation loss during epoch {epoch}: {mean_val_loss_per_epoch[-1]:>7f}")

        early_stopping(mean_val_loss_per_epoch[-1],model)

        if early_stopping.early_stop:
            print("Terminating early as stopping condition has been met")
            break

        loss_df = createTrainDataFrame(epoch, mean_train_loss_per_epoch, mean_val_loss_per_epoch)
        loss_df.to_csv('loss.csv',index=False)
        writeTrainingCurve(name, loss_df)

    # store the overfit model, for debugging
    print('test')
    torch.save(model.state_dict(),'final_epoch.pt')

    # load the checkpoint with the lowest validation loss
    model.load_state_dict(torch.load('checkpoint.pt'))

    return model, loss_df
