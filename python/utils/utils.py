import glob
import json
import pickle
import pathlib
import subprocess

import numpy as np
import pandas as pd

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import one_to_three 
from Bio.PDB import Superimposer
from Bio.PDB import Structure,Model
from Bio.PDB import NeighborSearch
import os
import warnings
import shutil

## for loading structures from multi-entry PDB file and threading a new sequence on
def getPDBsFromMultiEntryFile(pdbNames,multientryPDBpath,outDir):
    pdbNames_set = set(pdbNames)
    pdbNames_found = set()
    print(f"Searching for {len(pdbNames_set)} structures in multi-entry PDB file")
    # open multi-entry pdb file
    with open(multientryPDBpath,'r') as file:
        # get lines corresponding to selected PDBs
        filenames = []
        for line in file:
            if line[:6] == 'HEADER':
                name = line[10:].rstrip()
                if name not in pdbNames_set:
                    continue
                if name in pdbNames_found:
                    print(f"Warning: {name} was observed twice, there are duplicate names in the PDB file")
                pdbNames_found.add(name)
                pdb_lines = []
                pdb_line = file.readline()
                while pdb_line[:3] != 'END':
                    if pdb_line[:3] != 'TER':
                        pdb_lines.append(pdb_line)
                    pdb_line = file.readline()
                with open(os.path.join(outDir,name+'.pdb'),'w') as pdb_file:
                    for pdb_line in pdb_lines:
                        pdb_file.write(pdb_line)
                filenames.append(name)

    print(f"Found {len(pdbNames_found)} structures in the file")
    # return the list of paths
    return filenames
    
def loadAndModifyPDBs(filtered_df,multientryPDBpath,targetPDBchains=None,prefix='',
                      binder_chain_id='0',name_col='name',
                      seq_col='',score_col='',renumber=False):
    # extract the PDBs of interest from the multi-entry file
    originalPDBDir = 'extractedPDBs' if prefix == '' else f"{prefix}_extractedPDBs"
    os.makedirs(originalPDBDir,exist_ok=True)
    filenames = getPDBsFromMultiEntryFile(filtered_df[name_col],multientryPDBpath,originalPDBDir)
    
    # load with Bio.PDB parser
    modifiedPDBDir = 'modifiedPDBs' if prefix == '' else f"{prefix}_modifiedPDBs"
    os.makedirs(modifiedPDBDir,exist_ok=True)
    io = PDBIO()
    p = PDBParser(PERMISSIVE=1)
    warnings.filterwarnings('ignore')
    for i,row in filtered_df.iterrows():
        name = row[name_col]
        binder_structure = p.get_structure(name,os.path.join(originalPDBDir,name+'.pdb'))
        
        # set sequence/bfactor
        if seq_col != '':
            for i,residue in enumerate(binder_structure[0][binder_chain_id]):
                residue.resname = one_to_three(row[seq_col][i])
                if score_col == '':
                    continue
                for j,atom in enumerate(residue):
                    atom.set_bfactor(row[score_col])
                
        # combine with target chain 
        if targetPDBchains:
            binder_chain = binder_structure[0][binder_chain_id]
            binder_chain.detach_parent()
            if renumber:
                binder_chain.id = 'A'
            model = Model.Model(0)
            new_binder_structure = Structure.Structure('')
            new_binder_structure.add(model)
            new_binder_structure[0].add(binder_chain)
            if renumber:
                chain_names = ['B','C','D','E']
                residue_idx_counter = len(list(binder_chain.get_residues())) + 1
                for i,chain in enumerate(targetPDBchains):
                    for res in chain.get_residues():
                        res.detach_parent()
                        res.id = (res.id[0],residue_idx_counter,' ')
                        residue_idx_counter+=1
                    chain.detach_parent()
                    chain.id = chain_names[i]
                    new_binder_structure[0].add(chain)
            else:
                for i,chain in enumerate(targetPDBchains):
                    new_binder_structure[0].add(chain)
            binder_structure = new_binder_structure
        
        # write modified file
        io.set_structure(binder_structure)
        io.save(os.path.join(modifiedPDBDir,modifiedPDBDir+'_'+name+'.pdb'))
    warnings.filterwarnings('default')

# For comparing structure predictions to designs
def getBBAtoms(residue):
    bb_atoms = []
    all_atoms = {a.name:a for a in residue.get_atoms()}
    for atom_name in ['N','CA','C','O']:
        if atom_name in all_atoms:
            bb_atoms.append(all_atoms[atom_name])
        else:
            raise ValueError(f"Backbone atom {atom_name} not found in {residue}")
    return bb_atoms

def getBBAtomsFromChain(chain):
    return [atom for residue in chain.get_residues() for atom in getBBAtoms(residue)]
        
def calcRMSD(atoms1,atoms2):
    atoms1 = np.array(atoms1)
    atoms2 = np.array(atoms2)
    return np.sqrt(np.mean(np.square(atoms1-atoms2)))

def alignStructureByChains(fixed_chain,mobile_chain,mobile,sup):
    # find the optimal superposition between the selected chains
    fixed_atoms = getBBAtomsFromChain(fixed_chain)
    mobile_atoms = getBBAtomsFromChain(mobile_chain)
    
    assert len(fixed_atoms) == len(mobile_atoms)
    
    sup.set_atoms(fixed_atoms,mobile_atoms)
    print('Aligned chains with RMSD:',sup.rms)
    
    # apply to the whole mobile structure
    sup.apply([a for a in mobile.get_atoms()])
    return sup.rms

def getRMSD(design_structure,
            design_target_chain_id,design_peptide_chain_id,
            alphafold_structure,
            alphafold_target_chain_id,alphafold_peptide_chain_id,
            sup=None ):
    if sup == None:
        sup = Superimposer()

    # Superimpose by target chain (E in design_structure, C in alphafold_structure)
    design_target_chain = design_structure[0][design_target_chain_id]
    alphafold_target_chain = alphafold_structure[0][alphafold_target_chain_id]
    alignStructureByChains(design_target_chain,alphafold_target_chain,alphafold_structure,sup)
    
    # Extract peptide chains (A in design_structure, B in alphafold_structure)
    design_peptide_chain = design_structure[0][design_peptide_chain_id]
    alphafold_peptide_chain = alphafold_structure[0][alphafold_peptide_chain_id]
    
    # Calculate RMSD over backbone atoms
    return calcRMSD(getBBAtomsFromChain(design_peptide_chain),getBBAtomsFromChain(alphafold_peptide_chain))

def loadAlphaFoldResults(dir_path: str):
    glob_path = os.path.join(dir_path,'output/*/ranking_debug.json')
    paths = glob.glob(glob_path)
    print(f"{len(paths)} paths, ex: {paths[0]}")

    name_list = []
    path_list = []
    itptm_list = []
    plddt_list = []

    for path in paths:
        with open(path,'r') as file:
            name = path.split('/')[-2]
            ranking = json.load(file)
            pdb_path = path[:-len('/ranking_debug.json')]+'/unrelaxed_'+ranking['order'][0]+'.pdb'
            pickle_path = path[:-len('/ranking_debug.json')]+'/result_'+ranking['order'][0]+'.pkl'
            with open(pickle_path,'rb') as file:
                dump = pickle.load(file)
            plddt = dump['plddt'][0:len(dump['seqs'][0])].mean()
            
            name_list.append(name)
            path_list.append(pdb_path)
            itptm_list.append(ranking['iptm+ptm'][ranking['order'][0]])
            plddt_list.append(plddt)
            
    bind_ptm_df = pd.DataFrame({'name':name_list,
                        'path':path_list,
                        'iptm+ptm':itptm_list,
                        'plddt':plddt_list})
    return bind_ptm_df

def findSeedsCloseToSelection(structure_path,chain_id_list,binder_scoring_iterable_dataset,distance_cutoff=5):
    # load structure to search
    parser = PDBParser(QUIET=False)
    structure = parser.get_structure("",structure_path)
    atoms_to_search_against = [a for chain_id in chain_id_list for a in structure[0][chain_id].get_atoms()]
    ns = NeighborSearch(atoms_to_search_against)

    with open("atom_overlaps.csv",'w') as file:
        file.write(f"name,n_atom_overlaps"+"\n")

        # iterate over seeds in dataloader
        # dataset = BinderScoringIterableDataset(args.binder_dataset,args.target_pdb,27500,binder_subset,complex_only=complex_only)
        dataset_iter = iter(binder_scoring_iterable_dataset)
        assert binder_scoring_iterable_dataset.load_binder == True and binder_scoring_iterable_dataset.load_complex == False 
        for binder_data in dataset_iter:
            for i,seed in enumerate(binder_data):
                name = seed[0]['pdb']
                coords = seed[0]['coords'] # coords have dim of L x 4 x 3
                ca_coords = coords[:,0] # L x 3
                nearby_list = []
                for atom_coord in ca_coords:
                    # search each Ca carbon for proximal atoms
                    nearby_list += ns.search(atom_coord,distance_cutoff)
                nearby_set = set(nearby_list)
                # print(f"{name},{len(nearby_set)}")
                file.write(f"{name},{len(nearby_set)}"+"\n")


def hamming_distance(s1,s2):
    assert len(s1) == len(s2)
    return sum(a != b for a, b in zip(s1, s2))

def selectTopNDistinctSequences(sequences,N,hamming_distance_cutoff=1,verbose=False):
    ''' Select N sequences with >= hamming distance

    This function assumes that sequences are sorted in order of precedence
    '''
    selected_seqs = []
    for current_seq in sequences:
        unique = True
        for selected_seq in selected_seqs:
            if hamming_distance(current_seq,selected_seq) < hamming_distance_cutoff:
                # sequence is redundant to one already selected
                unique = False
                continue
        if unique:
            selected_seqs.append(current_seq)
        if len(selected_seqs) >= N:
            break

    if verbose:
        print(f"Found {len(selected_seqs)} sequences with hamming distance >= {hamming_distance_cutoff}")
    return selected_seqs
