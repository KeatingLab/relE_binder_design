import os
import glob
import math
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from .data_utils import *

parser = argparse.ArgumentParser(description='Load data tables and construct training data tables with random 80/10/10 train/val/test split')
parser.add_argument('dataPath', type=str,
                    help='The path to the data tables (include wildcards if there are multiple file)')
parser.add_argument('dataName', type=str,
                    help='The unique name of the provided dataset')
args = parser.parse_args()

if __name__ == "__main__":
    aaprob_df = loadDataAndCreateMasterDataFrame(args.dataName, args.dataPath)
    splitData(aaprob_df, args.dataName)

    print("Done!")