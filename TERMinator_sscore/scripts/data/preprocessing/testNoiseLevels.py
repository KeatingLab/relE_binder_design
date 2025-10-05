import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import argparse
import copy
from tqdm import tqdm

DIR = os.path.dirname(os.path.abspath(__file__))
assert DIR[0] == "/", "DIR should be an abspath"

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(DIR))), 'terminator', 'data'))
from data import load_file, get_flex, generate_noise, calc_dihedrals_bond_lengths

if __name__ == '__main__':
    # TODO: create argument for parser to take in folder w proper dir structure for pdbs
    print("starting main")
    parser = argparse.ArgumentParser('Run NRGTEN to calculate dynamic signatures.')
    parser.add_argument("--in_folder", help=" folder to find TERM file", required=True)
    parser.add_argument("--pdb_folder", help="folder to find raw pdb files", required=True)
    parser.add_argument("--out_folder", help="folder to hold results", required=True)
    parser.add_argument('--flex_folder', help='folder to find flex files', required=True)
    parser.add_argument('--pdb_list', help="list of pdbs", required=True)
    parser.add_argument('--verbose', help="whether to save updated structures and ramachandran plots for each pdb file", default=False)
    args = parser.parse_args()

    """
        in_folder : str
        folder to find TERM file.
    pdb_folder : str
        folder to find raw pdb files.
    flex_folder : str
        folder to find flex files.
    pdb_id : str
        PDB ID to load.
    flex_type : str
        methodology to calculate flex data
    noise_level : float
        std of noise to add
    bond_length_noise_level : float
        std of noise to add to bond lengths 
    min_protein_len : int
        minimum cutoff for loading TERM file.
    max_protein_len : int
        maximumum cutoff for loading protein to fit on GPU.
    num_ensembles : int
        number of conformational ensembles for each protein.
"""
    if not os.path.isdir(args.out_folder):
        os.mkdir(args.out_folder)
    with open(args.pdb_list, 'r') as f:
        pdb_list = f.readlines()
    pos_variation = np.ones((20, 20))
    pdb_list = [pdb_list[0]]
    for pdb_id in tqdm(pdb_list):
        for i_n, n in enumerate([0.1]):
            for i_b, b in enumerate([0]):
                pdb_id = pdb_id.strip()
                noise_level = 0.5
                bond_noise_level = 0
                res = load_file(args.in_folder, args.pdb_folder, args.flex_folder, pdb_id, flex_type=None, noise_level=0.0, 
                                                                        bond_length_noise_level=0, min_protein_len=30, max_protein_len=0, num_ensembles=1)
                if res is not None:
                    pdb_data, total_term_length, seq_len, flex = res
                orig_dihedrals, orig_bond_lengths = calc_dihedrals_bond_lengths(X=pdb_data['coords'], chain_lens=pdb_data['chain_lens'])
                all_new_dihedrals, all_new_bond_lengths = [], []
                for i in range(10):
                    noise = generate_noise(flex_type='random_torsion_batch', noise_level=n, size=pdb_data['coords'].shape, X=pdb_data['coords'], bond_length_noise_level=b, chain_lens=pdb_data['chain_lens'], expected_dihedrals=orig_dihedrals)
                    noise = noise.reshape((noise.shape[0]*noise.shape[1], noise.shape[2]))
                    pos_variation[i_n, i_b]= np.mean(abs(noise))
                    print(np.mean(abs(noise)))
                    if args.verbose:
                        filename = os.path.join(args.pdb_folder, pdb_id, pdb_id + ".red.pdb")
                        newlines = [f"MODEL       {str(i)}\n"]
                        i_atom = 0
                        with open(filename, 'r') as fp:
                            for line in fp:
                                data = line.strip()
                                if data[:3] == 'TER' or data[:6] == 'REMARK':
                                    if data[:3] == 'TER':
                                        continue
                                    if i == 0:
                                        newlines.append(line)
                                    continue
                                if data[:3] == 'END' and data[3:6] != 'MDL':
                                    newlines.append(line)
                                    continue
                                elif data[:5] == 'MODEL':
                                    cur_model = data[11:14].strip()
                                    newlines.append(line)
                                    continue
                                elif data[:6] == 'ENDMDL':
                                    newlines.append(line)
                                    continue
                                x = data[30:38].strip()
                                y = data[38:46].strip()
                                z = data[46:54].strip()
                                x = " " + str(round(float(x) + noise[i_atom, 0], len(x.split('.')[1]))).ljust(len(x), '0').rjust(len(data[30:37]))
                                y = " " + str(round(float(y) + noise[i_atom, 1], len(y.split('.')[1]))).ljust(len(y), '0').rjust(len(data[38:45]))
                                z = " " + str(round(float(z) + noise[i_atom, 2], len(z.split('.')[1]))).ljust(len(z), '0').rjust(len(data[46:53]))
                                data = data[:30] + x + y + z + data[54:] + '\n'
                                i_atom += 1
                                newlines.append(data)
                        out_filename = os.path.join(args.out_folder, pdb_id + ".pdb")
                        with open(out_filename, 'a') as f:
                            for line in newlines:
                                    f.write(line)
                            f.write("ENDMDL\n")
                        noise = noise.reshape((int(noise.shape[0]/4), 4, 3))
                        pdb_data['coords'] += noise
                        new_dihedrals, new_bond_lengths = calc_dihedrals_bond_lengths(X=pdb_data['coords'], chain_lens=pdb_data['chain_lens'], expected_bond_lengths=orig_bond_lengths, expected_dihedrals=orig_dihedrals)
                        if i == 0:
                            all_new_dihedrals, all_new_bond_lengths = new_dihedrals, new_bond_lengths
                        else:
                            all_new_dihedrals = np.concatenate((all_new_dihedrals, new_dihedrals))
                            all_new_bond_lengths = np.concatenate((all_new_bond_lengths, new_bond_lengths))
                    
                    ## Plot results
                    plt.scatter(orig_dihedrals[0], orig_dihedrals[1], c=len(orig_dihedrals[0])*[1], s=np.pi*3, alpha=0.5)
                    plt.title('Original Ramachandran plot')
                    plt.xlabel('phi')
                    plt.ylabel('psi')
                    axes = plt.gca()
                    y_range = axes.get_ylim()
                    x_range = axes.get_xlim()
                    plt.savefig(os.path.join(args.out_folder, pdb_id + "_orig_ramachandran.png"))
                    plt.close()

                    plt.scatter(all_new_dihedrals[0], all_new_dihedrals[1], c=len(all_new_dihedrals[0])*[1], s=np.pi*3, alpha=0.5)
                    plt.title('New Ramachandran plot')
                    plt.xlabel('phi')
                    plt.ylabel('psi')
                    axes = plt.gca()
                    axes.set_ylim(y_range)
                    axes.set_xlim(x_range)
                    plt.savefig(os.path.join(args.out_folder, pdb_id + "_new_ramachandran.png"))
                    plt.close()

                    _, bins, _ = plt.hist(orig_bond_lengths[0], alpha=0.5, label='orig', density=True)
                    plt.hist(all_new_bond_lengths[0], bins, alpha=0.5, label='new', density=True)
                    plt.legend(loc='upper right')
                    plt.title('N-CA Distances')
                    plt.show()
                    plt.savefig(os.path.join(args.out_folder, pdb_id + "_NCA.png"))
                    plt.close()

                    _, bins, _ = plt.hist(orig_bond_lengths[1], alpha=0.5, label='orig', density=True)
                    plt.hist(all_new_bond_lengths[1], bins, alpha=0.5, label='new', density=True)
                    plt.legend(loc='upper right')
                    plt.title('CA-C Distances')
                    plt.show()
                    plt.savefig(os.path.join(args.out_folder, pdb_id + "_CAC.png"))
                    plt.close()

                    _, bins, _ = plt.hist(orig_bond_lengths[2], alpha=0.5, label='orig', density=True)
                    plt.hist(all_new_bond_lengths[2], bins, alpha=0.5, label='new', density=True)
                    plt.legend(loc='upper right')
                    plt.title('C-N Distances')
                    plt.show()
                    plt.savefig(os.path.join(args.out_folder, pdb_id + "_CN.png"))
                    plt.close()

    #np.savetxt("full_test_noise_results.csv", pos_variation, delimiter=",")


