from pymol import cmd
import glob

def loadAlignedFrames(pathToDir,cleanPrevious=True):
    r"""Method for loading frame-aligned interacting residues and grouping

    Parameters
    ----------
    var1: str
        Path to the directory containing the frame-aligned PDBs

    """
    pass

    # clean up from previous call
    if cleanPrevious:
        cmd.delete('*_aligned')

    pdb_paths = glob.glob(pathToDir+"*.pdb")
    pdb_paths.sort()
    # cmd.set('grid_max',len(pdb_paths))

    for i,path in enumerate(pdb_paths):
        name = path.split('/')[-1]
        aaName = name.split('_')[0]
        group_name = aaName+'_*'
        print(name,aaName)
        cmd.load(path)
        cmd.group(aaName+'_aligned',group_name)

    cmd.set('grid_mode',1)
    return
