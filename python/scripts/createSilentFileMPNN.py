import sys,os,pathlib,math
import argparse

import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.SeqIO.FastaIO import SimpleFastaParser

sys.path.insert(0,"/data1/groups/keatinglab/swans/repos/peptide_binder_design")
from python.src.utils.pdbutils import multiPDBLoader

def getSeqsFromFasta(path):
    name = ''
    seq_list = []
    with open(path,"r") as file:
        for header,sequence in SimpleFastaParser(file):
            header_split = header.split(',')
            if len(header_split) == 8:
                # native seq line
                name = header_split[0]
            elif len(header_split) == 5:
                # designed seq line
                seq_list.append(sequence)
    assert name != '' and seq_list != []
    return name,seq_list

def removePrefix(name,prefix):
    return name[len(prefix):]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load MPNN sequences and peptide structures and create rosetta silent files')
    parser.add_argument('--mpnn_dir', type=str, help='Directory where MPNN was run')
    parser.add_argument('--multientry_pdb', type=str, help='')
    parser.add_argument('--silent_name', type=str, help='')
    parser.add_argument('--target_pdb', type=str, default='')
    parser.add_argument('--remove_name_prefix', type=str, default='')
    parser.add_argument('--n_batches', type=int, help='')
    parser.add_argument('--batch_id', type=int, help='')
    parser.add_argument('--silentfrompdbs', type=str, help='', default='/home/gridsan/sswanson/local_code_mirror/silent_tools/silentfrompdbs')
    args = parser.parse_args()
    print(sys.argv)

    # Load all fasta files from directory and select files that are pertinent to this batch
    dir_path = os.path.join(args.mpnn_dir,"outputs/seqs/")
    fasta_list = os.listdir(dir_path)
    assert len(fasta_list) > 0
    fasta_list.sort() #for consistency
    print(f"Found {len(fasta_list)} fasta files in the MPNN directory, ex: {fasta_list[0]}")
    batch_size = math.ceil(len(fasta_list) / (args.n_batches))
    # batch_size = len(fasta_list) // (args.n_batches)
    fasta_batch_list = fasta_list[args.batch_id*batch_size:args.batch_id*batch_size+batch_size]
    fasta_name_list = [removePrefix(pathlib.Path(x).stem,args.remove_name_prefix) for x in fasta_batch_list]
    fasta_path_list = [os.path.join(dir_path,x) for x in fasta_batch_list]

    # load target pdb
    target_structure = None
    if (args.target_pdb != ""):
        p = PDBParser(PERMISSIVE=1)
        target_structure = p.get_structure('',args.target_pdb)

    name2seqdict = dict()
    for fasta_path in fasta_path_list:
        name,seq_list = getSeqsFromFasta(fasta_path)
        name = removePrefix(name,args.remove_name_prefix)
        assert name not in name2seqdict
        name2seqdict[name] = seq_list

    silent_name = args.silent_name + "_" + str(args.batch_id)
    if target_structure is not None:
        with multiPDBLoader(args.multientry_pdb) as loader:
            name2seq_df = loader.writeToRosettaSilent(fasta_name_list,
                                        name2seqdict,
                                        silent_name,
                                        target_structure,
                                        silentfrompdbs=args.silentfrompdbs)
    else:
        with multiPDBLoader(args.multientry_pdb) as loader:
            name2seq_df = loader.writeToRosettaSilent(fasta_name_list,
                                        name2seqdict,
                                        silent_name,
                                        silentfrompdbs=args.silentfrompdbs)
            
    name2seq_df.to_csv(f"name2seq_{args.batch_id}.csv")

    print("Done!")