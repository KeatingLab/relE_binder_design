import glob
import json
import random

import pandas as pd
import seaborn as sns
import numpy as np


# Functions for processing residue pair statistics, calculating structure scores, and creating feature files

def groupByResPairType(data_df):
    # group by close-in-chain and close-in-space 
    data_df = data_df[data_df['Ca_distance']<12.5].copy(deep=True)
    cic_df = data_df[(data_df['distance_in_chain']<=8)&(data_df['same_chain']==True)]
    cis_df = data_df[(data_df['distance_in_chain']>8)|(data_df['same_chain']==0)]
    
    print(f"{len(cic_df)} close-in-chain residue pairs and {len(cis_df)} residue pairs")
    return cic_df,cis_df

def getStatisticsForReferenceState(data_df,prefix='',bins=[0,6,12.5,2000]):
    cic_df,cis_df = groupByResPairType(data_df)
    cic_grouped = cic_df.groupby('distance_in_chain').median().reset_index()
    distanceinchain_data = dict(zip(cic_grouped['distance_in_chain'],cic_grouped['n_matches']))

    print(distanceinchain_data)

    # define bins
    # bins = [0,6,12.5,2000]
    cis_df['Ca_distance_bin'] = pd.cut(cis_df['Ca_distance'],bins)
    cis_df

    cis_grouped = cis_df.groupby('Ca_distance_bin').median().reset_index()
    distanceinspace_data = dict(zip(cis_grouped['Ca_distance_bin'],cis_grouped['n_matches']))
    print(distanceinspace_data)

    # store as json
    cic_path = prefix + '_distanceinchain_data.json'
    with open(cic_path,'w') as file:
        json.dump(distanceinchain_data,file)

    cis_path = prefix + '_distanceinspace_data.json'
    with open(cis_path,'w') as file:
        json.dump(distanceinspace_data,file)

    distanceinchain_data,distanceinspace_data

def scoreCloseInChain(n_matches,distance_in_chain,data,pseudocount=1):
    return -np.log((n_matches+pseudocount)/(data[distance_in_chain]+pseudocount))

# deprecated: now CIS interactions are classified as SS and BB
def scoreCloseInSpace(n_matches,Ca_distance_interval,data,pseudocount=1):
    return -np.log((n_matches+pseudocount)/(data[Ca_distance_interval]+pseudocount))

def getScoresForAllResiduePairs():
    # get scores
    score_list = []
    type_list = []
    for i,row in df.iterrows():
        if i % 100000 == 0:
            print(i)
        if row['same_chain'] and row['distance_in_chain'] <= 7:
            # close in chain
            score = scoreCloseInChain(row['n_matches'],row['distance_in_chain'],distanceinchain_data)
            score_type = 'cic'
        else:
            # close in space
            score = scoreCloseInSpace(row['n_matches'],row['Ca_distance_bin'],distanceinspace_data)
            score_type = 'cis'
            
        score_list.append(score)
        type_list.append(score_type)
        
    df['score'] = score_list
    df['score_type'] = type_list
    return df

def aggPairScores():
    # new approach for aggregating
    pdb_id_list = []
    ri_residx_list = []
    ri_resnum_list = []
    ri_chainid_list = []
    ri_neighbors_list = []
    score_list = []
    for i,(pdb_id,prot_df) in enumerate(df.groupby('pdb_id')):
        if i % 100 == 0:
            print(i)
        for (resIdx,chainID,resnum),res_pairs_df in prot_df.groupby(['Ri_resIdx','Ri_chainID','Ri_resnum']):
            # compute close in chain score
            cic_df = res_pairs_df[res_pairs_df['score_type']=='cic']
            cic_scores = np.array(cic_df.score)
            cic_chain_dist = np.array(cic_df.distance_in_chain)
            cic_weights = np.exp(-cic_chain_dist/10)
            cic_weights = cic_weights / cic_weights.sum()
            cic_score = np.inner(cic_scores,cic_weights)

            # compute close in space score
            pseudovalue = 0.0
            cis_scores = np.array(list(res_pairs_df[res_pairs_df['score_type']=='cis'].score) + [pseudovalue])
            cis_score = cis_scores.mean()

            score = (cic_score + cis_score) / 2
            
            pdb_id_list.append(pdb_id)
            ri_residx_list.append(resIdx)
            ri_resnum_list.append(resnum)
            ri_chainid_list.append(chainID)
            ri_neighbors_list.append(len(res_pairs_df))
            score_list.append(score)
            
    # # group by Ri (need to duplicate contacts for this to properly work)
    # res_score_df = df.groupby(['pdb_id','Ri_resIdx','Ri_chainID','Ri_resnum']).agg(
    # #     mean_score = pd.NamedAgg(column='score',aggfunc='mean'),
    #     size = pd.NamedAgg(column='score',aggfunc='size')
    # ).reset_index()
    # res_score_df = res_score_df.sort_values(by='Ri_resIdx')
    # res_score_df['score'] = score_list
    # res_score_df

    # # group by Ri (need to duplicate contacts for this to properly work)
    # res_score_df = df.groupby(['pdb_id','Ri_resIdx','Ri_chainID','Ri_resnum']).agg(
    #     mean_score = pd.NamedAgg(column='score',aggfunc='mean'),
    #     size = pd.NamedAgg(column='score',aggfunc='size')
    # ).reset_index()
    # res_score_df

    res_score_df = pd.DataFrame({'pdb_id':pdb_id_list,
                                'Ri_resIdx':ri_residx_list,
                                'Ri_chainID':ri_chainid_list,
                                'Ri_resnum':ri_resnum_list,
                                'Ri_neighbors':ri_neighbors_list,
                                'score':score_list})

    # res_score_df = pd.DataFrame({'pdb_id':pdb_id_list,
    #                              'Ri_resIdx':ri_residx_list,
    #                              'score':score_list})
    res_score_df.to_csv('230301_dataset_res_structscore.csv')


def getResIdxToResInfoDict(df,pdb_id):
    r2dict = {row['Ri_resIdx']:(row['Ri_chainID'],row['Ri_resnum']) for i,row in df[df['pdb_id']==pdb_id].iterrows()}
    r2dict2 = {row['Rj_resIdx']:(row['Rj_chainID'],row['Rj_resnum']) for i,row in df[df['pdb_id']==pdb_id].iterrows()}
    return {**r2dict,**r2dict2}

def writeScoreDFtoCSV(res_score_df,pdb_id,df,suffix=''):
    r2dict = getResIdxToResInfoDict(df,pdb_id)
    low_score_df = res_score_df[(res_score_df['pdb_id']==pdb_id)][['pdb_id','Ri_resIdx','mean_score']].copy(deep=True)
    chain_list = []
    resnum_list = []
    for i,row in low_score_df.iterrows():
        ri_info = r2dict[row['Ri_resIdx']]
        chain_list.append(ri_info[0])
        resnum_list.append(ri_info[1])
    low_score_df['chain_id'] = chain_list
    low_score_df['resnum'] = resnum_list
    low_score_df.to_csv(f"230111_{pdb_id}_score.csv") if suffix == '' else low_score_df.to_csv(f"230111_{pdb_id}_{suffix}_score.csv")
    

# Functions for creating feature files

'''
1. Load existing feature files for those structures in the subsets
2. Add fine-tuning data to each
3. Copy to a new directory

One complication: sometimes, data is missing for a residue in the structure. The most certain workaround is to load the red.pdb using the original parseCoords function, get the residue info, and then perform some kind of merge (where any missing value is 0).
'''

def createResIdxDF(res_info):
    idx_list = []
    chain_list = []
    resnum_list = []
    for i,(chain_id,resnum) in enumerate(res_info):
#         print(i,chain_id,res_num)
        idx_list.append(int(i))
        chain_list.append(chain_id)
        resnum_list.append(str(resnum))
        
    return pd.DataFrame({'Ri_resIdx':idx_list,'Ri_chainID':chain_list,'Ri_resnum':resnum_list})

def createResIdxDFfromScoreDF(df,res_score_df,pdb_id):
    r2dict = getResIdxToResInfoDict(df,pdb_id)
    low_score_df = res_score_df[(res_score_df['pdb_id']==pdb_id)][['pdb_id','Ri_resIdx','score']].copy(deep=True)
    low_score_df = low_score_df.sort_values(by='Ri_resIdx')
    chain_list = []
    resnum_list = []
    for i,row in low_score_df.iterrows():
        ri_info = r2dict[row['Ri_resIdx']]
        chain_list.append(ri_info[0])
        resnum_list.append(str(ri_info[1]))
    low_score_df['Ri_chainID'] = chain_list
    low_score_df['Ri_resnum'] = resnum_list
    return low_score_df
    
def mergeScoreInfoWithResInfo(pdb_id,df,res_score_df):
    try:
        _, _, res_info = parseCoords(f"/data1/groups/keatinglab/fosterb/data/multichain_clean/{pdb_id}/{pdb_id}.red.pdb",
                False,False,False,False)
    except:
        return None
    df1 = createResIdxDF(res_info)
    df2 = createResIdxDFfromScoreDF(df,res_score_df,pdb_id)
    df2 = df2.drop_duplicates(['Ri_chainID','Ri_resnum'],keep='first')
# #     display(df1)
#     print(df1.dtypes)
# #     display(df2)
#     print(df2.dtypes)
    return df1.merge(df2[['pdb_id','Ri_chainID','Ri_resnum','score']],on=['Ri_chainID','Ri_resnum'],how='left')


def createFeatureFiles():
    path_to_data = '/data1/groups/keatinglab/fosterb/data/multichain_features'
    path_to_new_data = '/data1/groups/keatinglab/swans/TERMinator/230218_finetune_multichain_features'
    subsets = {'train_filt_subsample':'/data1/groups/keatinglab/fosterb/data/multichain_features/train_filt_subsample.in',
            'validation_subsample':'/data1/groups/keatinglab/fosterb/data/multichain_features/validation_subsample.in',
            'test_subsample':'/data1/groups/keatinglab/fosterb/data/multichain_features/test_subsample.in'}

    for name,path in subsets.items():
        print(name,path)
        # Get the PDB IDs corresponding to this subset
        with open(path,'r') as file:
            pdb_ids = [x.rstrip() for x in file]
        print(len(pdb_ids)," total structures")
        for pdb_id in pdb_ids:
            print(pdb_id)
            old_feature_path = os.path.join(path_to_data,pdb_id,pdb_id+'.features')
            old_len_path = os.path.join(path_to_data,pdb_id,pdb_id+'.length')
            new_dir_path = os.path.join(path_to_new_data,pdb_id)
            new_feature_path = os.path.join(new_dir_path,pdb_id+'.features')
            new_len_path = os.path.join(new_dir_path,pdb_id+'.length')
            # create new dir
            try:
                os.mkdir(new_dir_path)
            except FileExistsError:
                pass
            if os.path.isfile(new_feature_path) and os.path.isfile(new_len_path):
                print('Feature files exist, skipping...')
                continue
            # properly merge structure score with residue info
            filt_df = mergeScoreInfoWithResInfo(pdb_id,df,res_score_df)
            if filt_df is None:
                print(f"Issue parsing {pdb_id}, skipping...")
                continue
            # copy the length file
            shutil.copy(old_len_path,new_len_path)
            # load the feature file and add structure score info
            with open(old_feature_path,'rb') as file:
                feature_dict = pickle.load(file)
            assert len(filt_df) == len(feature_dict['sequence'])
            feature_dict['sscore'] = np.array(list(filt_df['score']))
            # write to new location
            with open(new_feature_path,'wb') as file:
                pickle.dump(feature_dict,file)
    #         raise ValueError