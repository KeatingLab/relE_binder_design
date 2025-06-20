import argparse
import os
import regex as re
import glob
import json
import pickle
import pathlib
import subprocess
import shutil
import warnings
import sys
import time
from copy import deepcopy
import string
import random

import pandas as pd
from Bio.File import as_handle
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import one_to_three 
from Bio.PDB import PDBIO
from Bio.PDB.Chain import Chain
from Bio.PDB.PDBExceptions import PDBConstructionWarning

import numpy as np
import pandas as pd

# Methods using Bio.PDB library to alter PDB structures

def modifyStructure(structure, 
                    chain_ids: set, 
                    renumber_chain_ids: dict = None,
                    remap_chain_ids: dict = None):
    # filter out chains
    for chain in structure[0].get_chains():
        if chain.id not in chain_ids:
            structure[0].detach_child(chain.id)
    
    # renumber THEN remap chains
    if renumber_chain_ids is not None:
        for chain in structure[0].get_chains():
            renumberChain(chain,renumber_chain_ids[chain.id])
    if remap_chain_ids is not None:
        remapChains(structure,remap_chain_ids)

def remapChains(structure,remap_chain_ids: dict):
    """Switch the chain IDs in a structure

    Arguments:
    - structure: the Bio.PDB.Structure to be altered
    - remap_chain_ids: dictionary containing chain ids to be swapped (e.g. remap_chain_ids['A'] = 'B')
    """
    # verify that the new chain IDs are unique
    assert len(remap_chain_ids) == len(set(remap_chain_ids.values()))
    
    # detach chains before modifying to avoid chain id clashes
    remapped_chains = list()
    for chain in structure[0].get_chains():
        if chain.id not in remap_chain_ids:
            raise ValueError(f"{chain.id} is not in the structure")
        structure[0].detach_child(chain.id)
        chain.id = remap_chain_ids[chain.id]
        remapped_chains.append(chain)
    for chain in remapped_chains:
        structure[0].add(chain)

def renumberChain(chain,start_resnum: int):
    """Set the residue number of a chain in ascending order starting from `start_resnum`
    """
    chain_res = list(chain.get_residues())
    for i,res in enumerate(chain_res):
        chain.detach_child(res.id)
        res.id = (res.id[0],start_resnum+i,' ')
    for res in chain_res:
        chain.add(res)

def replaceChainSequence(chain, aa_seq: str):
    assert len(aa_seq) == len(chain)
    for i,residue in enumerate(chain.get_residues()):
        residue.resname = one_to_three(aa_seq[i])

def setChainBFactor(chain, values: list):
    assert len(values) == len(chain)
    for i,residue in enumerate(chain.get_residues()):
        residue['CA'].set_bfactor(values[i])

# Methods for writing/reading to a silent file

def writeStructuresToSilentFile(structures: list, 
                                silent_binary_name: str, 
                                silentfrompdbs: str = '/data1/groups/keatinglab/swans/repos/silent_tools/silentfrompdbs'):
    # Write temporary PDB files
    warnings.filterwarnings('ignore')
    io = PDBIO()
    PDBDir = silent_binary_name + '_temp'
    shutil.rmtree(PDBDir,ignore_errors=True)
    pathlib.Path(PDBDir).mkdir(parents=True, exist_ok=True)
    for structure in structures:
        io.set_structure(structure)
        io.save(os.path.join(PDBDir,f"{structure.get_id()}.pdb"))

    # External call to silentfrompdbs
    command = f"{silentfrompdbs} {PDBDir}/*.pdb > {silent_binary_name}_binder_designs.silent"
    CompletedProcess = subprocess.run(command,shell=True,env=os.environ.copy())
    if CompletedProcess.returncode != 0:
        raise ValueError('Subprocess returned an error when trying to create a silent file')

    # delete PDBs and work on next batch
    shutil.rmtree(silent_binary_name + '_temp')
    warnings.filterwarnings('default')


def readStructuresFromSilentFile(silent_dir: str, 
                                 multipdbname: str,
                                 name_subset: list,
                                 name_remap: dict):
    
    # Before writing: check if the file already exists and contains all of the files
    if os.path.exists(f"{multipdbname}.pdb"):
        # File exists, check if it contains all of the names
        with multiPDBLoader(f"{multipdbname}.pdb") as loader:
            if name_remap:
                new_name_subset = [name_remap[name] for name in name_subset]
                if loader.areStructuresInFile(new_name_subset):
                    print(f"All structures already exist in {multipdbname}.pdb, skip")
                    return
                else:
                    print(f"{multipdbname}.pdb is missing some files, delete and remake it")
                    os.remove(f"{multipdbname}.pdb")
            elif loader.areStructuresInFile(name_subset):
                print(f"All structures already exist in {multipdbname}.pdb, skip")
                return
            else:
                print(f"{multipdbname}.pdb is missing some files, delete and remake it")
                os.remove(f"{multipdbname}.pdb")
                # raise ValueError("File exists, but does not contain all the structures. If you want to replace, delete the file")
    else:
        print(f"{multipdbname}.pdb does not exist yet...")

    # Get all silent files in the repo and concatenate into a single file
    silent_file_paths = glob.glob(silent_dir+'/*.silent')
    assert len(silent_file_paths) > 0
    print(f"Found {len(silent_file_paths)} silent files in the directory")

    tmp_dir_name = "temp_dir"
    pathlib.Path(tmp_dir_name).mkdir(parents=False, exist_ok=True)

    # External call to concatenate the files
    catsilentfile = f"{tmp_dir_name}/{multipdbname}_cat.silent"
    command = f"cat {' '.join(silent_file_paths)} > {catsilentfile}"
    CompletedProcess = subprocess.run(command,shell=True,env=os.environ.copy())
    if CompletedProcess.returncode != 0:
        raise ValueError('Subprocess returned an error when trying to concatenate a silent file')
    print(f"Concatenated all files into {catsilentfile}")

    with multiPDBWriter(f"{multipdbname}.pdb") as loader:
        # Extract the files in batches (to avoid clogging the directory with too many individual PDBs)
        batch_size = 100
        for batch_index_start in range(0,len(name_subset),batch_size):
            print(len(name_subset))
            print(batch_index_start,batch_index_start+batch_size)
            name_subset_batch = name_subset[batch_index_start:batch_index_start+batch_size]
            assert len(name_subset_batch) > 0
            print(f"Adding batch of {len(name_subset_batch)} pdb structures")
            
            # Extract PDBs from the silent file
            command = f"silentextractspecific {catsilentfile} {' '.join(name_subset_batch)}"
            CompletedProcess = subprocess.run(command,shell=True,env=os.environ.copy())
            if CompletedProcess.returncode != 0:
                raise ValueError('Subprocess returned an error when trying to extract structures from a silent file')
            
            extracted_pdb_paths = [f"{x}.pdb" for x in name_subset_batch]
            for path in extracted_pdb_paths:
                shutil.move(path,f"{tmp_dir_name}/{path}")
            
            # Add to the multiPDB
            extracted_pdb_paths = [f"{tmp_dir_name}/{x}.pdb" for x in name_subset_batch]

            if name_remap:
                new_name_list = [name_remap[name] for name in name_subset_batch]
                loader.addStructuresFromPDBs(extracted_pdb_paths,new_name_list)
            else:
                loader.addStructuresFromPDBs(extracted_pdb_paths)

            # Remove the PDB files
            for path in extracted_pdb_paths:
                os.remove(path)
    
    os.remove(catsilentfile)
    os.remove(str(f"{catsilentfile}.idx"))
    os.rmdir(tmp_dir_name)

# Extend the PDBParser class to parse lines from a file to enable reading of multiPDB files
class PDBParserEXT(PDBParser):

    def __init__(self,
                PERMISSIVE=True,
                get_header=False,
                structure_builder=None,
                QUIET=True,
                is_pqr=False):
        """Create a modified PDBParser object.

        Documentation copied from https://github.com/biopython/biopython/blob/master/Bio/PDB/PDBParser.py

        The PDB parser call a number of standard methods in an aggregated
        StructureBuilder object. Normally this object is instantiated by the
        PDBParser object itself, but if the user provides his/her own
        StructureBuilder object, the latter is used instead.

        Arguments:
         - PERMISSIVE - Evaluated as a Boolean. If false, exceptions in
           constructing the SMCRA data structure are fatal. If true (DEFAULT),
           the exceptions are caught, but some residues or atoms will be missing.
           THESE EXCEPTIONS ARE DUE TO PROBLEMS IN THE PDB FILE!.
         - get_header - unused argument kept for historical compatibility.
         - structure_builder - an optional user implemented StructureBuilder class.
         - QUIET - Evaluated as a Boolean. If true, warnings issued in constructing
           the SMCRA data will be suppressed. If false (DEFAULT), they will be shown.
           These warnings might be indicative of problems in the PDB file!
         - is_pqr - Evaluated as a Boolean. Specifies the type of file to be parsed.
           If false (DEFAULT) a .pdb file format is assumed. Set it to true if you
           want to parse a .pqr file instead.

        """
        super().__init__(PERMISSIVE,get_header,structure_builder,QUIET,is_pqr)


    def get_structure(self, id, file = None, lines = None):
        """Return the structure.

        Arguments:
         - id - string, the id that will be used for the structure
         - file - name of the PDB file OR an open filehandle
         - lines - a list of line strings from a PDB file
         
        """
        if file is None and lines is None:
            raise ValueError("Must provide a file or lines")
        with warnings.catch_warnings():
            if self.QUIET:
                warnings.filterwarnings("ignore", category=PDBConstructionWarning)

            self.header = None
            self.trailer = None
            # Make a StructureBuilder instance (pass id of structure as parameter)
            self.structure_builder.init_structure(id)

            if file is not None:
                with as_handle(file) as handle:
                    lines = handle.readlines()
                    if not lines:
                        raise ValueError("Empty file.")
                    self._parse(lines)
            else:
                self._parse(lines)

            self.structure_builder.set_header(self.header)
            # Return the Structure instance
            structure = self.structure_builder.get_structure()

        return structure

class multiPDBLoader:
    """Load structures from multiPDB file

    MultiPDB files are convenient because they group distinct PDBs into a single file and can be directly read by PyMol
    """

    def __init__(self, 
                 file_path: str):
        self.path = file_path
        self.fileh = open(self.path,'r')
        self.parser = PDBParserEXT()
        self.idx2pos = dict()
        self.idx2name = dict()
        self.name2pos = dict()
        self.maxpos = -1
        self.maxidx = -1

        print("Scanning multi-PDB file...")
        start = time.time()
        self.__getFilePositions()
        stop = time.time()
        print(f"Scanned the file in {(stop-start):.3f} s and found {self.maxidx + 1} structures")

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.fileh.close() 

    def getNStructures(self):
        return self.maxidx + 1
    
    def areStructuresInFile(self, structure_names: list):
        for name in structure_names:
            if not self.isStructureInFile(name):
                return False
        return True

    def isStructureInFile(self, structure_name: str):
        if structure_name in self.name2pos:
            return True
        return False

    def loadPDBByName(self, name: str):
        if not self.isStructureInFile(name):
            raise ValueError(f"name: {name} not in multiPDB file")
        return self.__loadPDB(self.name2pos[name],name)

    def loadPDBByIdx(self, idx: int):
        if idx not in self.idx2pos:
            raise ValueError(f"Structure with idx: {idx} not in multiPDB file")
        return self.__loadPDB(self.idx2pos[idx],self.idx2name[idx])
    
    def copyToNewMultiPDB(self, name_subset: set, new_filename: str, name_remapping: dict = None):
        if len(name_subset) > 100:
            print("Warning: large multiPDB files will slow down PyMol")
        with multiPDBWriter(new_filename) as file:
            for name in name_subset:
                self.fileh.seek(self.name2pos[name])
                lines = self.__readLines()
                if name_remapping is not None:
                    name = name_remapping[name]
                file.addStructureLines(name,lines)

    def sampleToNewMultiPDB(self, nToSample: int, new_filename: str, name_remapping: dict = None):
        if nToSample > 100:
            print("Warning: large multiPDB files will slow down PyMol")
        if nToSample > self.getNStructures():
            raise ValueError(f"Cannot sample > {self.getNStructures()} structures")
        # Get a list of random indices 
        all_indices = list(range(self.getNStructures()))
        random.shuffle(all_indices)
        sel_indices = all_indices[:nToSample]
        with multiPDBWriter(new_filename) as file:
            for idx in sel_indices:
                self.fileh.seek(self.idx2pos[idx])
                lines = self.__readLines()
                name = self.idx2name[idx]
                if name_remapping is not None:
                    name = name_remapping[name]
                file.addStructureLines(name,lines)

    def writeToRosettaSilent(self,
                             names: list,
                             name2seqdict: dict,
                             silent_binary_name: str,
                             target_structure = None,
                             remap_names: dict = None,
                             silentfrompdbs: str = '/data1/groups/keatinglab/swans/repos/silent_tools/silentfrompdbs',
                             binder_chain_id: str = '0'):
        
        print(f"Load {len(names)} structures and add to silent file")

        # Load each unique backbone structure 
        name_subset = set(names)
        name2structure = {x:self.loadPDBByName(x) for x in name_subset}

        # We create a new structure for each backbone + sequence pair
        name_list = []
        all_seq_list = []
        structure_list = []
        for name,seq_list in name2seqdict.items():
            for i,seq in enumerate(seq_list):
                # copy before modifying object
                structure = deepcopy(name2structure[name])
                
                # set sequence to designed seq
                replaceChainSequence(structure[0][binder_chain_id],seq)

                # NOTE: if remapping names, don't add an extra sequence index
                new_name = f"{structure.get_id()}_{i}" if remap_names is None else f"{remap_names[structure.get_id()]}_{i}"
                
                name_list.append(new_name)
                all_seq_list.append(seq)

                structure.id = new_name

                ## Reformat for compatibility with alphafold initial guess
                # renumber binder chain
                modifyStructure(structure,{binder_chain_id},{binder_chain_id:1},{binder_chain_id:'A'})
                
                # add target chain, renumbering according to length of binder chain
                if target_structure is not None:
                    target_structure_copy = deepcopy(target_structure)
                    resnum_start = len(structure[0]['A'])+1
                    for i,target_chain in enumerate(target_structure_copy.get_chains()):
                        target_structure_copy[0].detach_child(target_chain.id)
                        
                        # renumber and alter chain id
                        renumberChain(target_chain,resnum_start)
                        target_chain.id = string.ascii_uppercase[i+1]

                        structure[0].add(target_chain)
                        resnum_start += len(target_chain)

                structure_list.append(structure)

        if silentfrompdbs != "":
            writeStructuresToSilentFile(structure_list,silent_binary_name,silentfrompdbs)
        else:
            print("No path was provided to silentFromPDBs, will skip running")
        return pd.DataFrame({'name':name_list,'sequence':all_seq_list})
    

    def writeExactToRosettaSilent(self,
                                backbone_names: list,
                                seqs: list,
                                design_names: list,
                                silent_binary_name: str,
                                target_structure = None,
                                silentfrompdbs: str = '/data1/groups/keatinglab/swans/repos/silent_tools/silentfrompdbs',
                                binder_chain_id: str = '0'):
        
        assert len(backbone_names) == len(seqs) == len(design_names)
        assert len(design_names) == len(set(design_names))
        
        # Load each unique backbone structure 
        bb_name_subset = set(backbone_names)
        print(f"Load {len(bb_name_subset)} structures and add to silent file")
        bbname2structure = {x:self.loadPDBByName(x) for x in bb_name_subset}

        # We create a new structure for each backbone + sequence pair
        structure_list = []
        for backbone_name,seq,design_name in zip(backbone_names,seqs,design_names):
            # copy before modifying object
            structure = deepcopy(bbname2structure[backbone_name])
            
            # set sequence to designed seq
            replaceChainSequence(structure[0][binder_chain_id],seq)

            structure.id = design_name

            ## Reformat for compatibility with alphafold initial guess
            # renumber binder chain
            modifyStructure(structure,{binder_chain_id},{binder_chain_id:1},{binder_chain_id:'A'})
            
            # add target chain, renumbering according to length of binder chain
            if target_structure is not None:
                target_structure_copy = deepcopy(target_structure)
                resnum_start = len(structure[0]['A'])+1
                for i,target_chain in enumerate(target_structure_copy.get_chains()):
                    target_structure_copy[0].detach_child(target_chain.id)
                    
                    # renumber and alter chain id
                    renumberChain(target_chain,resnum_start)
                    target_chain.id = string.ascii_uppercase[i+1]

                    structure[0].add(target_chain)
                    resnum_start += len(target_chain)

            structure_list.append(structure)

        if silentfrompdbs != "":
            writeStructuresToSilentFile(structure_list,silent_binary_name,silentfrompdbs)
        else:
            print("No path was provided to silentFromPDBs, will skip running")
        return pd.DataFrame({'name':design_names,'sequence':seqs})

    def __getFilePositions(self):
        # Read through the whole file one time to record all of the contents
        self.fileh.seek(0)
        line = self.fileh.readline()
        if line == '' or line == 'FILEEND\n':
            print("File is empty")
            self.maxpos = -1
            self.maxidx = -1
            return
        idx = 0
        while line != '':
            header = self.__readHeader(line)
            if header != False:
                start_pos = self.fileh.tell()
                self.idx2pos[idx] = start_pos
                self.idx2name[idx] = header
                self.name2pos[header] = start_pos
                idx += 1
            line = self.fileh.readline()
        self.maxpos = start_pos
        self.maxidx = idx

    def __loadPDB(self, pos, name):
        if pos < 0 or pos > self.maxpos:
            raise ValueError(f"pos: {pos} not in multiPDB file")
        self.fileh.seek(pos)
        lines = self.__readLines()
        assert len(lines) > 0
        return self.parser.get_structure(name,None,lines)
    
    def __readHeader(self,line):
        # slightly faster check to eliminate non-matching lines
        if line[:6] != "HEADER":
            return False
        header_split = line.split()
        if header_split[0] == 'HEADER' and len(header_split) == 2:
            return header_split[1].rstrip()
        else:
            return False
    
    def __readLines(self):
        lines = []
        line = ''
        while re.match("^END",line) is None:
            line = self.fileh.readline()
            lines.append(line)
        return lines

class multiPDBWriter:
    """For writing multiPDB files
    """

    def __init__(self,file_path):
        self.path = file_path
        self.fileh = open(self.path,'w')
        self.writer = PDBIO()

    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.fileh.write("FILEEND\n")
        self.fileh.close()

    def addStructuresFromPDBs(self, pdb_paths: list, names: list = None, selChainID: str = ''):
        assert len(pdb_paths) > 0
        if names is not None:
            assert len(pdb_paths) == len(names)
        parser = PDBParserEXT()
        for i,pdb_path in enumerate(pdb_paths):
            if i % max(1,int(len(pdb_paths)/10)) == 0:
                print(f"Adding structure {i}/{len(pdb_paths)}...")
            if names is not None:
                name = names[i]
                structure = parser.get_structure(name,pdb_path)
            else:
                name = pathlib.Path(pdb_path).stem
                structure = parser.get_structure(name,pdb_path)
            if selChainID != '':
                chains = list(structure.get_chains())
                for chain in chains:
                    if chain.get_id() != selChainID:
                        structure[0].detach_child(chain.get_id())
                if len(structure) == 0:
                    raise ValueError(f"No chains matching {selChainID} in {structure.get_id()}")
            self.addStructure(structure,name)

    def addStructure(self, structure, name):
        self.writer.set_structure(structure)
        self.__writeHeader(name)
        self.writer.save(self.fileh)

    def addStructureLines(self, name, lines):
        self.__writeHeader(name)
        for line in lines:
            self.fileh.write(line)

    def __writeHeader(self,name):
        self.fileh.write(f"HEADER    {name}"+"\n")

# class silentFileWriter:
#     """For writing rosetta silent files
    
    
#     """


