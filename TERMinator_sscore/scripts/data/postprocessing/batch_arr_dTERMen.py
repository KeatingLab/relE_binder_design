"""Submit a batch array job on SLURM to run multiple batched dTERMen jobs.

Usage:
    .. code-block::

        python batch_arr_dTERMen.py \\
            --output_dir <dir_containing_etabs_folder> \\
            --pdb_root <pdb_root> \\
            --dtermen_data <dtermen_data_root> \\
            [--batch_size <batch_size>]

See :code:`python batch_arr_dTERMen.py --help` for more info.
"""

import argparse
import glob
import os
import subprocess
import random
import sys

# for autosummary import purposes
sys.path.insert(0, os.path.dirname(__file__))
from search_utils import find_pdb_path

DIR = os.path.dirname(os.path.abspath(__file__))
assert DIR[0] == "/", "DIR should be an abspath"

sys.path.insert(0, os.path.join(os.path.dirname(DIR), 'preprocessing'))
from parseCoords import parseCoords

def write_data(index, pdbs, batch_size):
    """Writes pdb information to new file and returns file name and num_batches

    Args
    ----
    index : int
        Index of dTERMen job batch
    pdbs : list
        List of pdb names to include in batch
    batch_size : int
        Number of dTERMen runs to run per node

    Returns
    -------
    batch_arr_list : str
        Path to list of pdb files to run
    num_batches : int
        Info on number of cpus needed for each node
    """
    batch_arr_list = os.path.join(args.output_dir, f"{basename}_batch_arr_{i_job_batch}.list")

    with open(batch_arr_list, 'w') as fp:
        random.shuffle(pdbs)
        for pdb in pdbs:
            fp.write(pdb + "\n")
    fp.close()

    num_batches = (len(pdbs) // batch_size) + 1
    return batch_arr_list, num_batches


if __name__ == '__main__':
    # TODO: create argument for parser to take in folder w proper dir structure for pdbs
    parser = argparse.ArgumentParser('Run dTERMen for testing.')
    parser.add_argument('--output_dir', help='Output directory', required=True)
    parser.add_argument("--pdb_root", help="The root for all raw data PDB databases", required=True)
    parser.add_argument('--dtermen_data', help="Root directory for dTERMen runs", required=True)
    parser.add_argument('--batch_size', help='number of dTERMen runs to run per node', default=5, type=int)
    parser.add_argument('--job_batch_size', help='number of dTERMen jobs to run at once', default=100, type=int)
    parser.add_argument('--dtermen_type', help='type of dtermen etab scoring to perform', default='enerTable')
    args = parser.parse_args()
    os.chdir(DIR)

    pdbs = []
    output_path = os.path.join(args.output_dir, 'etabs')
    basename = os.path.basename(args.output_dir)
    batch_arr_lists = []
    num_batch_list = []
    i_job_batch = 0
    i_pdb = 0
    script_lists = []
    script_list = []
    pdb_lists = []
    print(f"Starting to look for files in {output_path}.")
    for filename in glob.glob(os.path.join(output_path, '*.etab')):
        pdb_id = os.path.basename(filename)[:-5]
        output_path_test = f"{output_path}/design_{pdb_id}.out"
        if os.path.isfile(output_path_test):
            with open(output_path_test, 'r') as f:
                lines = f.readlines()
            if len(lines) > 0 and lines[-1].find("recovery") > -1:
                continue
        print(pdb_id)
        if args.dtermen_type == 'design':
            pdb_path = find_pdb_path(pdb_id, args.pdb_root)
            os.system(f"cp {pdb_path} {output_path}/{pdb_id}.pdb")
            os.system((f"sed -e \"s|ID|{pdb_id}|g\" "
                    f"-e \"s|OUTPUTDIR|{output_path}|g\" "
                    f"-e \"s|POSTDIR|{DIR}|g\" "
                    f"< run_dTERMen.sh "
                    f" >{output_path}/run_{pdb_id}.sh"))
        elif args.dtermen_type == 'enerTable':
            pdb_path = find_pdb_path(pdb_id, args.dtermen_data)
            _, seq = parseCoords(pdb_path, False, False, False, False)
            os.system(f"cp {pdb_path} {output_path}/{pdb_id}.pdb")
            os.system((f"sed -e \"s|ETAB|{filename}|g\" "
                     f"-e \"s|ID|{pdb_id}|g\" "
                    f"-e \"s|PDBDIR|{output_path}|g\" "
                    f"-e \"s|SEQUENCE|{seq}|g\" "
                    f"< run_designsequence.sh "
                    f" >{output_path}/run_{pdb_id}.sh"))
        if i_pdb > args.job_batch_size:
            print("finishing old batch and starting new batch")
            batch_arr_list, num_batches = write_data(i_job_batch, pdbs, args.batch_size)
            batch_arr_lists.append(batch_arr_list)
            num_batch_list.append(num_batches)
            script_lists.append(script_list)
            pdb_lists.append(pdbs)
            script_list = []
            pdbs = []
            i_job_batch += 1
            i_pdb = 0
        i_pdb += 1
        script_list.append(f"{output_path}/run_{pdb_id}.sh")
        pdbs.append(pdb_id)
    if pdbs:
        batch_arr_list, num_batches = write_data(i_job_batch, pdbs, args.batch_size)
        batch_arr_lists.append(batch_arr_list)
        num_batch_list.append(num_batches)
    # bid = -1
    print(batch_arr_lists, num_batch_list)
    # for i_job_batch, (batch_arr_list, script_list, pdb_list) in enumerate(zip(batch_arr_lists, script_lists, pdb_lists)):
    #     print(f"Running batch {i_job_batch}.")
    #     process_list = []
    #     for script, pdb_id in zip(script_list, pdb_list):
    #         print(script)
    #         if args.dtermen_type == 'design':
    #             process = subprocess.Popen(
    #                 (f"sbatch --parsable --array=0-{num_batch_list[i_job_batch]} --mem={args.batch_size * 2000} --mincpu={args.batch_size} "
    #                 f"{os.path.join(DIR, 'batch_arr_dTERMen.sh')} {output_path} {batch_arr_list} {args.batch_size}")).read()
    #         elif args.dtermen_type == 'enerTable':
    #             with open(os.path.join(output_path, 'design_'+pdb_id+'.out'), 'wb') as out, open(os.path.join(output_path, 'design_'+pdb_id+'.err'), 'wb') as err:
    #                 process = subprocess.Popen(["/bin/bash", script], stdout=out, stderr=err)
    #         process_list.append(process)
    #     for process in process_list:
    #         process.wait()
        # os.system(f"sbatch --dependency=afterany:{bid} sum_res.sh {args.output_dir} {args.dtermen_data}")
