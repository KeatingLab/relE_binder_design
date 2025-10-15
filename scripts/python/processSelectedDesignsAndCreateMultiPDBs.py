import json, sys, os, pathlib, re
from pathlib import Path
import shutil
import glob
from copy import deepcopy

import pandas as pd
import matplotlib.pyplot as plt

from Bio.Seq import Seq
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import one_to_three 
from Bio.PDB import PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqIO.FastaIO import SimpleFastaParser


sys.path.insert(0,'/data1/groups/keatinglab/swans/repos/peptide_binder_design')
from python.src.utils.pdbutils import multiPDBLoader,readStructuresFromSilentFile


def main():

    ### Params
    path_to_csv="/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/230712_relEbinderdesigns_final6000_aminoaciddnasequences_todia.csv"
    final_df = pd.read_csv(path_to_csv,index_col=0)

    data_path_dict = {
        'denovo_2seed_coordinator':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/1_createSilentFiles/denovo_2seed_coordinator/createsilent',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/1_createSilentFiles/denovo_2seed_coordinator/relax',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/denovo_2seed_coordinator'
        },
        'denovo_2seed_mpnn':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_thesis/rosetta_analysis/0_silentFromPDBs_denovo2seedMPNN',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_thesis/rosetta_analysis/1_relax_denovo2seedMPNN',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/denovo_2seed_mpnn'
        },
        'denovo_3seed_coordinator':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/1_createSilentFiles/denovo_3seed_coordinator/createsilent',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/1_createSilentFiles/denovo_3seed_coordinator/relax',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/denovo_3seed_coordinator'
        },
        'denovo_3seed_mpnn':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_3seeds/rosetta_analysis/denovoMPNN/0_silentFromPDBs_denovoMPNN',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_3seeds/rosetta_analysis/denovoMPNN/1_relax_denovoMPNN',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/denovo_3seed_mpnn'
        },
        'denovopartialdiff_2seed_coordinator':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/rf_diffusion/partialdiffusion_diversifyingbackbones/rosetta_analysis/0_silentFromPDBs_2',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/rf_diffusion/partialdiffusion_diversifyingbackbones/rosetta_analysis/1_relax_2',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/denovopartialdiff_2seed_coordinator'
        },
        'denovopartialdiff_3seed_coordinator':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/rf_diffusion/partialdiffusion_diversifyingbackbones/rosetta_analysis/0_silentFromPDBs_2',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/rf_diffusion/partialdiffusion_diversifyingbackbones/rosetta_analysis/1_relax_2',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/denovopartialdiff_3seed_coordinator'
        },
        'relBextension_2seed_coordinator':{
            '1_fused':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_2seedextensions/rosetta_analysis/0_silentFromPDBs_2',
            '2_relaxed':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/round_4_2seedextensions/rosetta_analysis/1_relax_2',
            '3_af':'/data1/groups/keatinglab/swans/binderDesign_relE/analysis/combined_analysis/alphafold/2_afInitGuess/relBextension_2seed_coordinator'
        }
    }

    # NOTE: renaming the structures as pipeline_rank (the full names are garbled by the PyMol loader)
    name_remap = {row['name']:f"{row['pipeline']}_{int(row['rank'])}"for _,row in final_df.iterrows()}
    relax_name_remap = {f"{row['name']}_0001":f"{row['pipeline']}_{int(row['rank'])}"for _,row in final_df.iterrows()}
    alphafold_relax_name_remap = {f"{row['name']}_0001_af2pred":f"{row['pipeline']}_{int(row['rank'])}"for _,row in final_df.iterrows()}

    full_dir = "topranked100"
    pathlib.Path(full_dir).mkdir(exist_ok=True)
    os.chdir(full_dir)
    for pipeline in data_path_dict:
        subset_list = list(final_df[(final_df['pipeline']==pipeline)&(final_df['rank']<=100)]['name'])
        multipdbname = f"{pipeline}_fused_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['1_fused'],multipdbname,subset_list,name_remap)

        relaxed_subset_list = [f"{name}_0001" for name in subset_list]
        multipdbname = f"{pipeline}_relaxed_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['2_relaxed'],multipdbname,relaxed_subset_list,relax_name_remap)

        af_relaxed_subset_list = [f"{name}_0001_af2pred" for name in subset_list]
        multipdbname = f"{pipeline}_alphafold_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['3_af'],multipdbname,af_relaxed_subset_list,alphafold_relax_name_remap)
    os.chdir("..")

    full_dir = "random100"
    pathlib.Path(full_dir).mkdir(exist_ok=True)
    os.chdir(full_dir)
    for pipeline in data_path_dict:
        subset_list = list(final_df[(final_df['pipeline']==pipeline)].sample(100,random_state=100)['name'])
        multipdbname = f"{pipeline}_fused_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['1_fused'],multipdbname,subset_list,name_remap)

        relaxed_subset_list = [f"{name}_0001" for name in subset_list]
        multipdbname = f"{pipeline}_relaxed_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['2_relaxed'],multipdbname,relaxed_subset_list,relax_name_remap)

        af_relaxed_subset_list = [f"{name}_0001_af2pred" for name in subset_list]
        multipdbname = f"{pipeline}_alphafold_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['3_af'],multipdbname,af_relaxed_subset_list,alphafold_relax_name_remap)
    os.chdir("..")

    full_dir = "allstructures"
    pathlib.Path(full_dir).mkdir(exist_ok=True)
    os.chdir(full_dir)
    for pipeline in data_path_dict:
        subset_list = list(final_df[final_df['pipeline']==pipeline]['name'])
        multipdbname = f"{pipeline}_fused_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['1_fused'],multipdbname,subset_list,name_remap)

        relaxed_subset_list = [f"{name}_0001" for name in subset_list]
        multipdbname = f"{pipeline}_relaxed_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['2_relaxed'],multipdbname,relaxed_subset_list,relax_name_remap)

        af_relaxed_subset_list = [f"{name}_0001_af2pred" for name in subset_list]
        multipdbname = f"{pipeline}_alphafold_{full_dir}"
        print(multipdbname)
        readStructuresFromSilentFile(data_path_dict[pipeline]['3_af'],multipdbname,af_relaxed_subset_list,alphafold_relax_name_remap)
    os.chdir("..")

if __name__ == '__main__':
    main()
    print("Done!")