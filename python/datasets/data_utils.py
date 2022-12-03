import glob
import math
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

aa3 = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 
       'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
       'LEU', 'LYS', 'MET', 'PHE', 'PRO', 
       'SER', 'THR', 'TRP', 'TYR', 'VAL']
aa3toaa1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P', 
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}
aa1toaa3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP',
    'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'G': 'GLY',
    'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
    'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER',
    'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}
aa3tosurfProb = {
    'ALA':0.048049,
    'ARG':0.090171,
    'ASN':0.062437,
    'ASP':0.095417,
    'CYS':0.002533,
    'GLN':0.070898,
    'GLU':0.125364,
    'GLY':0.029287,
    'HIS':0.023232,
    'ILE':0.016676,
    'LEU':0.033629,
    'LYS':0.132090,
    'MET':0.009195,
    'PHE':0.012781,
    'PRO':0.070291,
    'SER':0.068929,
    'THR':0.058319,
    'TRP':0.006683,
    'TYR':0.016602,
    'VAL':0.027414
} # P(aa|rSASA > 0.05)

def loadDataAndCreateMasterDataFrame(dataName: str, dataPath: str):
    # Load all data
    print("Loading data...")
    #path = '/home/gridsan/sswanson/keatinglab_shared/swans/learningInterfaceScoringFunction/generateData/220824_interfaceScoreSimple/pixelDB728_uniqueBindingModes_path-targetchains-binderchain_*_trainingData.csv'
    paths = glob.glob(dataPath)

    df_list = list()
    for i,path in enumerate(paths):
        df = pd.read_csv(path)
        df_list.append(df)
    data_df = pd.concat(df_list,ignore_index=True)
    # data_df.to_csv('PixelDB728_uniqueBindingModes_allData.csv',index=False)
    data_df.to_csv(dataName+'_targetAAProb_allData.csv',index=False)
    print(len(data_df)," total datapoints")

    # for each contact, generate 20 possible contacts, each with a distinct target AA
    df_cols = ['originalIdx','targetResAA1','targetResAAIdx',
            'N-N', 'N-Ca', 'N-C', 'N-O', 'Ca-N', 'Ca-Ca', 'Ca-C','Ca-O', 
            'C-N', 'C-Ca', 'C-C', 'C-O', 'O-N', 'O-Ca', 'O-C', 'O-O',
            'nTotalMatches','nativeAA','Score']
    minNMatches = 20

    data_list = list()
    for i,row in data_df.iterrows():
        for j,aa in enumerate(aa3):
            data = list()
            data += [i,aa3toaa1[aa],j]
            data += list(row[6:22])
            data += [row['nTotalMatches']]
            data += [aa == row['targetResName']]
            if (row['nTotalMatches'] < minNMatches):
                data += [0.0]
            else:
                data += [-np.log2(((row[aa3toaa1[aa]]+1)/(row['nTotalMatches']+1))/(aa3tosurfProb[aa]))]
            data_list.append(data)
            
    aaprob_df = pd.DataFrame(data_list,columns=df_cols)
    
    sns.histplot(data=aaprob_df,x='Score',bins=30,binrange=[-15,15])
    plt.savefig(dataName+"_scores_histplot.png",dpi=300)

    return aaprob_df

def splitData(data_df: pd.DataFrame, dataName: str):
    # randomly split data into training and testing sets 
    shuffle_df = data_df.sample(frac=1,random_state=0).reset_index(drop=True)

    nDatapoints = len(shuffle_df)
    nTrain = math.floor(nDatapoints*0.8)

    train_df = shuffle_df.loc[:nTrain]
    val_df = shuffle_df.loc[nTrain+1:math.floor(nTrain+1+((len(shuffle_df)-nTrain)/2))]
    test_df = shuffle_df.loc[math.floor(nTrain+(len(shuffle_df)-nTrain)/2)+1:].reset_index(drop=True)

    train_df.to_csv(dataName+"_targetAAProb_train80.csv",index=False)
    val_df.to_csv(dataName+"_targetAAProb_val10.csv",index=False)
    test_df.to_csv(dataName+"_targetAAProb_test10.csv",index=False)

    return (train_df,val_df,test_df)