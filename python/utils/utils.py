import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import one_to_three 
from Bio.PDB import Superimposer
import os
import warnings
import shutil

import glob
import json
import pickle

## for loading structures from multi-entry PDB file and threading the sequence on

def getPDBsFromMultiEntryFile(pdbNames,multientryPDBpath,outDir):
    pdbNames_set = set(pdbNames)
    # open multi-entry pdb file
    with open(multientryPDBpath,'r') as file:
        # get lines corresponding to selected PDBs
        filenames = []
        for line in file:
            if line[:6] == 'HEADER':
                name = line[10:].rstrip()
                if name not in pdbNames_set:
                    continue
                print(name)
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
        
    # return the list of paths
    return filenames
    
def loadAndModifyPDBs(filtered_df,multientryPDBpath,targetPDBchain,prefix='',pdb_chain_id='A',
                        seq_col='sequence',name_col='name'):
    # extract the PDBs of interest from the multi-entry file
    originalPDBDir = 'extractedPDBs' if prefix == '' else f"{prefix}_extractedPDBs"
    os.makedirs(originalPDBDir,exist_ok=True)
    print(len(filtered_df[name_col]))
    filenames = getPDBsFromMultiEntryFile(filtered_df[name_col],multientryPDBpath,originalPDBDir)
    print(len(filenames),filenames[0:4])
    
    # load with Bio.PDB parser
    modifiedPDBDir = 'modifiedPDBs' if prefix == '' else f"{prefix}_modifiedPDBs"
    os.makedirs(modifiedPDBDir,exist_ok=True)
    io = PDBIO()
    p = PDBParser(PERMISSIVE=0)
    warnings.filterwarnings('ignore')
    for i,row in filtered_df.iterrows():
        name = row[name_col]
        s = p.get_structure(name,os.path.join(originalPDBDir,name+'.pdb'))
        
        # set sequence/bfactor
        for i,residue in enumerate(s[0][pdb_chain_id]):
            residue.resname = one_to_three(row[seq_col][i])
            for j,atom in enumerate(residue):
                atom.set_bfactor(row['score'])
                
        # combine with target chain
        s[0].add(targetPDBchain)
        
        # write modified file
        io.set_structure(s)
        io.save(os.path.join(modifiedPDBDir,name+'.pdb'))
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
            sup=None):
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
