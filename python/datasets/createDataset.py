import os
import glob
import math
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

from createDataset import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load residue pair statistics and generate structure scores')
    parser.add_argument('dataPath', type=str,
                        help='The path to the data tables (include wildcards if there are multiple file)')
    parser.add_argument('dataName', type=str,
                        help='The unique name of the provided dataset')
    args = parser.parse_args()

    # 


    print("Done!")