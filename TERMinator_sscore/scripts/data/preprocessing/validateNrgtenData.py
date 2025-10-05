"""Compare NRGTEN results from full-atom and backbone-only structures

Usage:
    .. code-block::

        python validateNrgtenData.py \\
            --full_root \\
            --full_list
            --backbone_root \\
            --backbone_list \\
            --out_dir \\
            --num_processes \\

"""

import argparse
import os
import multiprocessing as mp
from tqdm import tqdm
import pickle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from parseCoords import parseCoords


def calc_corr(x, y):
    """Calculates correlation between x and y
    Args
    ----
    x : list
        X data
    y : list
        Y data

    Returns
    -------
    corr : float
        Correlation between x and y
    """
    corr_matrix = np.corrcoef(x, y)
    corr = corr_matrix[0,1]
    corr = corr**2
    return corr

def plot_results(full, backbone, title, out_file, comp_type):
    
    """"Plots results with r^2 shown
    Args
        ----
        full : list
            Data to plot on x-axis
        backbone: list
            Data to plot on y-axis
        title : str
            Title for plot
        out_file : str
            Location to save plots
        comp_type : str
            Type of comparison to plot
    """

    def scatter_hist(x, y, ax, ax_histx, ax_histy, title_x, title_y):
        """Plots points with histograms on axes

        Args
        ----
        x : list
            Data to plot on x-axis
        y : list
            Data to plot on y-axis
        ax : plt.ax
            Axis to plot data on
        ax_histx : plt.ax
            Histogram for x-axis
        ax_histy : plt.ax
            Histogram for y-axis
        title_x : str
            Title for x-axis
        title_y : str
            Title for y-axis
        Returns 
        -------
        """
        ax_histx.tick_params(axis="x", labelbottom=title_x)
        ax_histy.tick_params(axis="y", labelleft=title_y)

        ax.scatter(x, y, s=1)
        ax.set_xlim([0, np.max(x)])
        ax.set_ylim([0, np.max(y)])

        binwidth = 0.25
        x_lim = (int(np.max(x)/binwidth) + 1) * binwidth
        y_lim = (int(np.max(y)/binwidth) + 1) * binwidth
        x_bins = np.arange(0, x_lim + binwidth, binwidth)
        y_bins = np.arange(0, y_lim + binwidth, binwidth)
        ax_histx.hist(x, bins=x_bins)
        ax_histy.hist(y, bins=y_bins, orientation='horizontal')

        R_sq = calc_corr(x, y)
        ax.text(0.95, 0.01, f"R^2: {R_sq}.",
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='black', fontsize=15)

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.005


    rect_scatter = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]

    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_axes(rect_scatter)
    ax_histx = fig.add_axes(rect_histx, sharex=ax)
    ax_histy = fig.add_axes(rect_histy, sharey=ax)
    ax.set_title(title)

    if comp_type == "structure":
        x_title = "Full atom"
        y_title = "Backbone only"
    elif comp_type == "model":
        x_title = "Dynamic signatures"
        y_title = "Conformational ensembles"
    elif comp_type.find("b_factor") > -1:
        x_title = "B-factors"
        y_title = "Dynamic signatures" if comp_type.find("sigs") > -1 else "Conformational ensembles"
    scatter_hist(full, backbone, ax, ax_histx, ax_histy, x_title, y_title)
    plt.savefig(out_file)

def compare_results(full_root, backbone_root, pdb_root, pdb, out_dir, analysis_level, verbose, analysis_type):
    """Calculates nrgten features for given pdb file

    Args
    ----
    full_root : str
        Root for NRGTEN results for full-atom structures
    backbone_root : str
        Root for NRGTEN results for backbone-only structures
    pdb_root : str
        Location of pdb structures with B-factor information
    pdb
        ID of current pdb file to load
    out_dir : str
        Location of directory to save outputs of analysis
    analysis_level : str
            Whether to analyze std at a residue, atomic, or dimension level
    verbose : bool
        Whether or not to save intermediate analysis steps
    analysis_type : str
        What type of data to analyze
    Returns
    -------
    full_dynam_sigs : list
        List of full-atom dynamic signatures values as floats
    full_conf_ensemble : list
        List of full-atom conf ensemble positional variation values as floats
    backbone_dynam_sigs : list
        List of backbone-only dynamic signatures values as floats
    full_conf_ensemble : list
        List of backbone-only conf ensemble positional variation values as floats
    """

    def parse_sig_file(sig_file):
        """Parses dynamic signatures file

        Args
        ----
        sig_file : str
            Path to dynamc signatures file

        Returns
        -------
        sigs : list
            List of dynamic signatures values as floats
        """
        with open(sig_file, 'r') as f:
            sigs = f.readlines()
        f.close()
        sigs = [float(sig) for sig in sigs]
        return sigs

    def parse_ensemble_file(ensemble_file, trim=False, analysis_level="residue"):
        """Parses conf ensemble file

        Args
        ----
        ensemble_file : str
            Path to conf ensemble file
        trim : bool
            Whether to trim ensemble file to only backbone atoms
        analysis_level : str
            Whether to analyze std at a residue, atomic, or dimension level

        Returns
        -------
        pos_variation : list
            List of positional variation values as floats
        """
        if trim:
            trimed_file = ensemble_file[:-4] + "_trim.pdb"
            print(f"trimming {trimed_file}")
            with open(trimed_file, 'w') as of:
                with open(ensemble_file, 'r') as f:
                    lines = f.readlines()
                for line in lines:
                    if line[:4] == "ATOM" and line[13:16].strip() not in ['C', 'CA', 'N', 'O']:
                        continue
                    of.write(line)
            of.close()
            f.close()
            ensemble_file = trimed_file
            print("trimming done")
        try:
            coords, _ = parseCoords(ensemble_file, ensemble=True, save=False, verbose=False)
            if len(coords) == 1:
                chain = next(iter(coords.keys()))
                coords_tensor = coords[chain]
            else:
                chains = sorted(coords.keys())
                coords_tensor = np.vstack([coords[c] for c in chains])
            pos_variation = np.std(coords_tensor, axis=-1)
            if analysis_level == "residue":
                pos_variation = np.mean(pos_variation[:,1,:], axis=1)
                pos_variation = pos_variation.flatten()
                # pos_variation = np.mean(pos_variation, axis=(1, 2))
            elif analysis_level == "atom":
                pos_variation = np.mean(pos_variation, axis=2)
                pos_variation = pos_variation.flatten()
            else:
                pos_variation = pos_variation.flatten()
            return pos_variation
        except Exception as e:
            print(f"{e} with {ensemble_file}")
            return []

    def parse_b_factor_file(b_factor_file, analysis_level="residue"):
        """Parses b_factor file

        Args
        ----
        b_factor_file : str
            Path to b_factor file
        analysis_level : str
            Whether to analyze std at a residue or atomic level

        Returns
        -------
        b_factors : list
            List of B-factors for each backbone atom
        """
        try:
            b_factors = parseCoords(b_factor_file, ensemble=False, save=False, verbose=False, record_b_factor=True)
            if len(b_factors) == 1:
                chain = next(iter(b_factors.keys()))
                b_factors_tensor = b_factors[chain]
            else:
                chains = sorted(b_factors.keys())
                b_factors_tensor = np.vstack([b_factors[c] for c in chains])
            if analysis_level == "residue":
                b_factors_tensor = b_factors_tensor[:,1]
            return b_factors_tensor.flatten(), len(list(b_factors.keys()))
        except Exception as e:
            print(f"Error {e} with {b_factor_file}.")
            return []

    if analysis_type.find("signatures") > -1:
        full_file_signatures = os.path.join(full_root, pdb[:-1] + "_dynamic_signatures.txt")
        backbone_file_signatures = os.path.join(backbone_root, pdb[:-1] + "_dynamic_signatures.txt")
        # backbone_file_signatures = os.path.join("/data1/groups/keatinglab/fosterb/data/multichain_nrgten_large_motion_789", pdb[:-1] + "_dynamic_signatures.txt") # JFM
    if analysis_type.find("ensembles") > -1:
        full_file_ensembles = os.path.join(full_root, pdb[:-1] + "_conf_ensemble.pdb")
        backbone_file_ensembles = os.path.join(backbone_root, pdb[:-1] + "_conformational_ensembles.pdb")
    if analysis_type.find("b_factors") > -1:
        b_factors_file = os.path.join(pdb_root, pdb[:-1], pdb[:-1] + ".red.pdb")
    try:
        array_lens = []
        chains = -1
        if analysis_type.find("signatures") > -1:
            full_sigs = parse_sig_file(full_file_signatures)
            backbone_sigs= parse_sig_file(backbone_file_signatures)
            array_lens += [len(full_sigs), len(backbone_sigs)]
        if analysis_type.find("ensembles") > -1:
            full_ensembles = parse_ensemble_file(full_file_ensembles, analysis_level=analysis_level)
            backbone_ensembles = parse_ensemble_file(backbone_file_ensembles, analysis_level=analysis_level)
            array_lens += [len(full_ensembles), len(backbone_ensembles)]
        if analysis_type.find("b_factors") > -1:
            b_factors, chains = parse_b_factor_file(b_factors_file, analysis_level=analysis_level)
            array_lens += [len(b_factors)]
        if max(array_lens) != min(array_lens):
            [print(f"Array length: {array_len}") for array_len in array_lens] 
            raise ValueError
        array_len = array_lens[0]
        if verbose:
            if not os.path.isdir(os.path.join(out_dir, pdb[:-1])):
                os.mkdir(os.path.join(out_dir, pdb[:-1]))


            # Plot results
            sigs_analyzed = False
            ensembles_analyzed = False
            to_return = []
            if analysis_type.find("signatures") > -1:
                dynamic_outfile = os.path.join(out_dir, pdb[:-1], pdb[:-1] + "_dynamic_sigs.png")
                plot_results(full_sigs, backbone_sigs, 'Dynamic Signatures', dynamic_outfile, comp_type="structure")
                to_return += full_sigs, backbone_sigs
                sigs_analyzed = True
            if analysis_type.find("ensembles") > -1:
                ensemble_outfile = os.path.join(out_dir, pdb[:-1], pdb[:-1] + "_conf_ensembles.png")
                plot_results(full_ensembles, backbone_ensembles, 'Conformation Ensembles', ensemble_outfile, comp_type="structure")
                to_return += full_ensembles, backbone_ensembles
                ensembles_analyzed = True
            if args.analysis_type.find("b_factors") > -1:
                if sigs_analyzed:
                    b_factor_sigs_full_outfile = os.path.join(out_dir, pdb[:-1], pdb[:-1] + "sigs_bfactors_full.png")
                    plot_results(b_factors, full_sigs, 'B-factors vs Sigs Full', b_factor_sigs_full_outfile, comp_type="b_factors_sigs")
                    b_factor_sigs_backbone_outfile = os.path.join(out_dir,  pdb[:-1], pdb[:-1] + "sigs_bfactors_backbone.png")
                    plot_results(b_factors, backbone_sigs, 'B-factors vs Sigs Backbone', b_factor_sigs_backbone_outfile, comp_type="b_factors_sigs")
                if ensembles_analyzed:
                    b_factor_ensembles_full_outfile = os.path.join(out_dir,  pdb[:-1], pdb[:-1] + "ensembles_bfactors_full.png")
                    plot_results(b_factors, full_ensembles, 'B-factors vs Ensembles Full', b_factor_ensembles_full_outfile, comp_type="b_factors_ensembles")
                    b_factor_ensembles_backbone_outfile = os.path.join(out_dir,  pdb[:-1], pdb[:-1] + "ensembles_bfactors_backbone.png")
                    plot_results(b_factors, backbone_ensembles, 'B-factors vs Ensembles Backbone', b_factor_ensembles_backbone_outfile, comp_type="b_factors_ensembles")
                to_return += b_factors,
            if sigs_analyzed and ensembles_analyzed:
                full_outfile = os.path.join(out_dir, pdb[:-1], pdb[:-1] + "_full_features.png")
                backbone_outfile = os.path.join(out_dir, pdb[:-1], pdb[:-1] + "_backbone_features.png")
                plot_results(full_sigs, full_ensembles, 'Full Features', full_outfile, comp_type="model")
                plot_results(backbone_sigs, backbone_ensembles, 'Backbone Features', backbone_outfile, comp_type="model")     
            return [to_return, pdb, array_len, chains]

    except Exception as e:
        print(f"Error {e} with {pdb[:-1]}.")
        return [[[], [], [], [], []], pdb, -1, -1]



if __name__ == '__main__':
    parser = argparse.ArgumentParser('Compare NRGTEN results from full-atom and backbone-only structures.')
    parser.add_argument("--full_root", help="Toot for NRGTEN results for full-atom structures", required=True)
    parser.add_argument("--full_list", help="List of full-atom structures", required=True)
    parser.add_argument("--backbone_root", help="Root for NRGTEN results for backbone-only structures", required=True)
    parser.add_argument("--backbone_list", help="List of backbone-only structures", required=True)
    parser.add_argument('--out_dir', help="Location of directory to save outputs of analysis", required=True)
    parser.add_argument("--pdb_dir", help="Location of pdb structures with B-factor information", default=None)
    parser.add_argument('--num_processes', help="Number of processes for parallelism", default=32)
    parser.add_argument('--analysis_level', help="Level to preform conformational ensemble analysis", default="residue")
    parser.add_argument('--verbose', help="Whether or not to save intermediate analysis steps", default=True)
    parser.add_argument('--update', help="Whether or not to recalculate", default=True)
    parser.add_argument('--analysis_type', help="What type of data to analyze", default="signatures_ensembles")
    args = parser.parse_args()
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)


    file_paths = []
    if args.analysis_type.find("signatures") > -1:
        full_dynamic_outfile = os.path.join(args.out_dir, "full_dynamic_sigs.list")
        backbone_dynamic_outfile = os.path.join(args.out_dir, "backbone_dynamic_sigs.list")
        file_paths += [full_dynamic_outfile, backbone_dynamic_outfile]
    if args.analysis_type.find("ensembles") > -1:
        full_ensemble_outfile = os.path.join(args.out_dir, "full_ensembles.list")
        backbone_ensemble_outfile = os.path.join(args.out_dir, "backbone_ensembles.list")
        file_paths += [full_ensemble_outfile, backbone_ensemble_outfile]
    if args.analysis_type.find("b_factors") > -1:
        b_factors_outfile = os.path.join(args.out_dir, "b_factors.list")
        file_paths += [b_factors_outfile]

    if args.update or not (sum([os.path.exists(file_path) for file_path in file_paths]) == len(file_paths)):
        print("starting validation")
        shared_pdbs = []
        with open(args.full_list, 'r') as f:
            shared_pdbs = set(f.readlines())
        f.close()
        with open(args.backbone_list, 'r') as f:
            lines = f.readlines()
            lines = [os.path.basename(line)[:4] + "\n" for line in lines]
            shared_pdbs = list(set(lines) & shared_pdbs)
        f.close()

        # with open(args.backbone_list, 'r') as f:
        #     shared_pdbs = list(set(f.readlines()) & shared_pdbs)
        # f.close()

        full_dynam_sigs_list = []
        backbone_dynam_sigs_list = []
        full_conf_ensembles_list = []
        backbone_conf_ensembles_list = []
        b_factors_list = []
        chains_list = []
        pdb_list = []
        pdb_lens_list = []
        with mp.Pool(int(args.num_processes)) as pool:
            progress = tqdm(total=len(shared_pdbs))

            def update_progress(res):
                del res
                progress.update(1)

            res_list = [
                pool.apply_async(compare_results, (args.full_root, args.backbone_root, args.pdb_dir, pdb, args.out_dir, args.analysis_level, args.verbose, args.analysis_type),
                                    callback=update_progress) for pdb in shared_pdbs
            ]
            pool.close()
            pool.join()
            progress.close()

            for res in res_list:
                data, pdb, pdb_len, chains = res.get()
                if data is not None:
                    sig_offset = 0
                    ensemble_offset = 0
                    if args.analysis_type.find("signatures") > -1:
                        full_dynam_sigs_list.append(data[0])
                        backbone_dynam_sigs_list.append(data[1])
                        sig_offset = 2
                    if args.analysis_type.find("ensembles") > -1:
                        full_conf_ensembles_list.append(data[sig_offset])
                        backbone_conf_ensembles_list.append(data[1 + sig_offset])
                        ensemble_offset = 2
                    if args.analysis_type.find("b_factors") > -1:
                        b_factors_list.append(data[ensemble_offset + sig_offset])
                        chains_list.append(chains)
                    pdb_list.append(pdb)
                    pdb_lens_list.append(pdb_len)
    

        # Save results
        if args.analysis_type.find("signatures") > -1:
            with open(full_dynamic_outfile, 'wb') as fp:
                pickle.dump({'data': full_dynam_sigs_list, 'pdbs': pdb_list, 'pdb_lens': pdb_lens_list}, fp)
            with open(backbone_dynamic_outfile, 'wb') as fp:
                pickle.dump({'data': backbone_dynam_sigs_list, 'pdbs': pdb_list, 'pdb_lens': pdb_lens_list}, fp)
        if args.analysis_type.find("ensembles") > -1:
            with open(full_ensemble_outfile, 'wb') as fp:
                pickle.dump({'data': full_conf_ensembles_list, 'pdbs': pdb_list, 'pdb_lens': pdb_lens_list}, fp)
            with open(backbone_ensemble_outfile, 'wb') as fp:
                pickle.dump({'data': backbone_conf_ensembles_list, 'pdbs': pdb_list, 'pdb_lens': pdb_lens_list}, fp)
        if args.analysis_type.find("ensembles") > -1:
            with open(b_factors_outfile, 'wb') as fp:
                pickle.dump({'data': b_factors_list, 'chains': chains_list, 'pdbs': pdb_list, 'pdb_lens': pdb_lens_list}, fp)
    else:
        if args.analysis_type.find("signatures") > -1:
            with open(full_dynamic_outfile, 'rb') as fp:
                data = pickle.load(fp)
                full_dynam_sigs = data['data']
            with open(backbone_dynamic_outfile, 'rb') as fp:
                data = pickle.load(fp)
                backbone_dynam_sigs = data['data']
        if args.analysis_type.find("signatures") > -1:
            with open(full_ensemble_outfile, 'rb') as fp:
                data = pickle.load(fp)
                full_conf_ensembles = data['data']
            with open(backbone_ensemble_outfile, 'rb') as fp:
                data = pickle.load(fp)
                backbone_conf_ensembles = data['data']
        if args.analysis_type.find("b_factors") > -1:
            with open(b_factors_outfile, 'rb') as fp:
                data = pickle.load(fp)
                b_factors = data['data']
                chain_list = data['chains']
        pdb_list = data['pdbs']
        pdb_lens_list = data['pdb_lens']

    # Save results in csv format and plot results
    print("analysis step 1")
    df = pd.DataFrame(columns=['PDB_ID','Length', 'Num_Chains', 'Dynamic_Signature_R2','Conf_Ensemble_R2','B_factors_vs_Sigs_1','B_factors_vs_Sigs_2','B_factors_vs_Ensembles_1','B_factors_vs_Ensembles_2', 'Sigs_vs_Ensembles_1', 'Sigs_vs_Ensembles_2'])
    df['PDB_ID'] = pdb_list
    df['Length'] = pdb_lens_list
    if len(full_dynam_sigs_list) > 0:
        full_dynam_sigs = np.concatenate(full_dynam_sigs_list)
    if len(backbone_dynam_sigs_list) > 0:
        backbone_dynam_sigs = np.concatenate(backbone_dynam_sigs_list)
    if len(full_conf_ensembles_list) > 0:
        full_conf_ensembles = np.concatenate(full_conf_ensembles_list)
    if len(backbone_conf_ensembles_list) > 0:
        backbone_conf_ensembles = np.concatenate(backbone_conf_ensembles_list)
    if len(b_factors_list) > 0:
        b_factors = np.concatenate(b_factors_list)
    print("analysis step 2")
    sigs_analyzed = False
    ensembles_analyzed = False
    if args.analysis_type.find("signatures") > -1:
        df['Dynamic_Signature_R'] = [calc_corr(x, y) for (x, y) in zip(full_dynam_sigs_list, backbone_dynam_sigs_list)]
        dynamic_outfile = os.path.join(args.out_dir, "dynamic_sigs.png")
        plot_results(full_dynam_sigs, backbone_dynam_sigs, 'Dynamic Signatures', dynamic_outfile, comp_type="structure")
        sigs_analyzed = True
    print("analysis step 3")
    if args.analysis_type.find("ensembles") > -1:
        df['Conf_Ensemble_R'] = [calc_corr(x, y) for (x, y) in zip(full_conf_ensembles_list, backbone_conf_ensembles_list)]
        ensemble_outfile = os.path.join(args.out_dir, "conf_ensembles.png")
        plot_results(full_conf_ensembles, backbone_conf_ensembles, 'Conformation Ensembles', ensemble_outfile, comp_type="structure")
        ensembles_analyzed = True
    print("analysis step 4")
    if args.analysis_type.find("b_factors") > -1:
        df['Num_Chains'] = chains_list
        if sigs_analyzed:
            df['B_factors_vs_Sigs_1'] = [calc_corr(x, y) for (x, y) in zip(b_factors_list, full_dynam_sigs_list)]
            df['B_factors_vs_Sigs_2'] = [calc_corr(x, y) for (x, y) in zip(b_factors_list, backbone_dynam_sigs_list)]
            b_factor_sigs_full_outfile = os.path.join(args.out_dir, "sigs_bfactors_full.png")
            plot_results(b_factors, full_dynam_sigs, 'B-factors vs Sigs Full', b_factor_sigs_full_outfile, comp_type="b_factors_sigs")
            b_factor_sigs_backbone_outfile = os.path.join(args.out_dir, "sigs_bfactors_backbone.png")
            plot_results(b_factors, backbone_dynam_sigs, 'B-factors vs Sigs Backbone', b_factor_sigs_backbone_outfile, comp_type="b_factors_sigs")
        if ensembles_analyzed:
            df['B_factors_vs_Ensembles_1'] = [calc_corr(x, y) for (x, y) in zip(b_factors_list, full_conf_ensembles_list)]
            df['B_factors_vs_Ensembles_2'] = [calc_corr(x, y) for (x, y) in zip(b_factors_list, backbone_conf_ensembles_list)]
            b_factor_ensembles_full_outfile = os.path.join(args.out_dir, "ensembles_bfactors_full.png")
            plot_results(b_factors, full_conf_ensembles, 'B-factors vs Ensembles Full', b_factor_ensembles_full_outfile, comp_type="b_factors_ensembles")
            b_factor_ensembles_backbone_outfile = os.path.join(args.out_dir, "ensembles_bfactors_backbone.png")
            plot_results(b_factors, backbone_conf_ensembles, 'B-factors vs Ensembles Backbone', b_factor_ensembles_backbone_outfile, comp_type="b_factors_ensembles")
    print("analysis step 5")
    print(f"full_dynam_sigs_list len: {len(full_dynam_sigs_list)}; full_conf_ensembles_list len: {len(full_conf_ensembles_list)}.")
    print(f"backbone_dynam_sigs_list len: {len(backbone_dynam_sigs_list)}; backbone_conf_ensembles_list len: {len(backbone_conf_ensembles_list)}.")
    print(f"full_dynam_sigs_list[0] len: {len(full_dynam_sigs_list[0])}; full_conf_ensembles_list len: {len(full_conf_ensembles_list[0])}.")
    print(f"backbone_dynam_sigs_list[0] len: {len(backbone_dynam_sigs_list[0])}; backbone_conf_ensembles_list len: {len(backbone_conf_ensembles_list[0])}.")
    if sigs_analyzed and ensembles_analyzed:
        df['Sigs_vs_Ensembles_1'] = [calc_corr(x, y) for (x, y) in zip(full_dynam_sigs_list, full_conf_ensembles_list)]
        df['Sigs_vs_Ensembles_2'] = [calc_corr(x, y) for (x, y) in zip(backbone_dynam_sigs_list, backbone_conf_ensembles_list)]
        full_outfile = os.path.join(args.out_dir, "full_features.png")
        backbone_outfile = os.path.join(args.out_dir, "backbone_features.png")
        plot_results(full_dynam_sigs, full_conf_ensembles, 'Full Features', full_outfile, comp_type="model")
        plot_results(backbone_dynam_sigs, backbone_conf_ensembles, 'Backbone Features', backbone_outfile, comp_type="model")
    df.to_csv(os.path.join(args.out_dir, 'full_stats.csv'))
    

    

