import sys
import argparse
import subprocess

import pandas as pd
from Bio.PDB.PDBParser import PDBParser

sys.path.insert(0,"/data1/groups/keatinglab/swans/repos/peptide_binder_design")
from python.src.utils.pdbutils import multiPDBLoader

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Load silent file and extract peptide sequence')
    parser.add_argument('--silent_file',type=str,help='')
    parser.add_argument('--target_sequence',type=str,help='this sequence is ignored',default='')
    # parser.add_argument('--n_batches', type=int, help='')
    parser.add_argument('--batch_id', type=int, help='')
    parser.add_argument('--silentsequence', type=str, help='', default='/data1/groups/keatinglab/swans/repos/silent_tools/silentsequence')
    args = parser.parse_args()
    print(sys.argv)

    # External call to silentfrompdbs
    command = f"{args.silentsequence} {args.silent_file}"
    CompletedProcess = subprocess.run(command,
                                      shell=True,
                                      capture_output=True,
                                      text=True)
    if CompletedProcess.returncode != 0:
        raise ValueError('Subprocess returned an error when trying to create a silent file')
    
    name_list = []
    sequence_list = []
    for line in str(CompletedProcess.stdout).split('\n'):
        if line == "":
            continue
        line_split = line.split(' ')
        assert len(line_split) > 2
        chain_sequences,name = line_split[:-1],line_split[-1]
        name_list.append(name)
        if len(chain_sequences) == 1:
            sequence_list.append(chain_sequences[0])
        else:
            non_target_sequences = [sequence for sequence in chain_sequences if sequence != args.target_sequence]
            assert len(non_target_sequences) == 1
            sequence_list.append(non_target_sequences[0])

    name2seq_df = pd.DataFrame({'name':name_list,'sequence':sequence_list})
    name2seq_df.to_csv(f"name2seq_{args.batch_id}.csv")

    print("Done!")