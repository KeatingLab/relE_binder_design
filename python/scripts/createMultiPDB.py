import sys,os,glob
import argparse

import pandas as pd
from Bio.PDB.PDBParser import PDBParser

sys.path.insert(0,"/data1/groups/keatinglab/swans/repos/peptide_binder_design")
from python.src.utils.pdbutils import multiPDBWriter

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load PDBs and add to multientry PDB file')
    parser.add_argument('--pdb_dir',type=str,required=True)
    parser.add_argument('--multientry_name', type=str, required=True)
    parser.add_argument('--chain_id',type=str,default='')
    # parser.add_argument('--target_pdb', type=str, help='')
    # parser.add_argument('--n_batches', type=int, help='')
    # parser.add_argument('--batch_id', type=int, help='')
    # parser.add_argument('--silentfrompdbs', type=str, help='', default='/home/gridsan/sswanson/local_code_mirror/silent_tools/silentfrompdbs')
    args = parser.parse_args()
    print(sys.argv)

    pdb_list = glob.glob(os.path.join(args.pdb_dir,"*.pdb"))
    assert len(pdb_list) > 0
    pdb_list.sort() #for consistency
    print(f"Found {len(pdb_list)} pdb files in the directory, ex: {pdb_list[0]}")

    fullpath_list = [os.path.join(args.pdb_dir,x) for x in pdb_list]
    
    with multiPDBWriter(args.multientry_name) as loader:
        loader.addStructuresFromPDBs(fullpath_list,None,args.chain_id)

    print("Done!")