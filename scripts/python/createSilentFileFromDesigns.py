import sys,math
import argparse

import pandas as pd
from Bio.PDB.PDBParser import PDBParser

sys.path.insert(0,"/data1/groups/keatinglab/swans/repos/peptide_binder_design")
from python.src.utils.pdbutils import multiPDBLoader

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load sequence and peptide structures and create rosetta silent files')
    parser.add_argument('--selected_designs', type=str, help='a CSV file containing a name and sequence column')
    parser.add_argument('--multientry_pdb', type=str, help='')
    parser.add_argument('--silent_name', type=str, help='')
    parser.add_argument('--target_pdb', type=str, help='')
    parser.add_argument('--n_batches', type=int, help='')
    parser.add_argument('--batch_id', type=int, help='')
    parser.add_argument('--binder_chain_id', type=str, help='')
    parser.add_argument('--name_col', type=str, help='')
    parser.add_argument('--seq_col', type=str, help='')
    parser.add_argument('--design_name_col',type=str,help='If there is already a name for each backbone + sequence, then use that')
    parser.add_argument('--silentfrompdbs', type=str, help='', default='/data1/groups/keatinglab/swans/repos/silent_tools/silentfrompdbs')
    args = parser.parse_args()
    print(sys.argv)

    # select backbones from the DF
    sel_df = pd.read_csv(args.selected_designs)
    batch_size = math.ceil(len(sel_df) / (args.n_batches))
    batch_df = sel_df[args.batch_id*batch_size:args.batch_id*batch_size+batch_size]

    # load target pdb
    p = PDBParser(PERMISSIVE=1)
    target_structure = p.get_structure('',args.target_pdb)

    silent_name = args.silent_name + "_" + str(args.batch_id)
    name2seqlist = {name:[seq for seq in group_df[args.seq_col]] for name,group_df in batch_df.groupby(args.name_col)}

    with multiPDBLoader(args.multientry_pdb) as loader:
        name2seq_df = loader.writeExactToRosettaSilent(batch_df[args.name_col],
                                                  batch_df[args.seq_col],
                                                  batch_df[args.design_name_col],
                                                  silent_name,
                                                  target_structure,
                                                  binder_chain_id=args.binder_chain_id,
                                                  silentfrompdbs=args.silentfrompdbs)
        
    name2seq_df.to_csv(f"name2seq_{args.batch_id}.csv")

    print("Done!")