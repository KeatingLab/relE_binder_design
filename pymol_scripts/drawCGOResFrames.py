from pymol import cmd, cgo, CmdException

'''

Method for loading the residue frames and drawing them as CGO arrows

'''

def drawCGOResFrames(pathToFile,cleanPrevious=True):
    r"""Loads the PDB, the residue frame file, and then draws the basis vectors

    After running 'buildGreedyClusters' seeds will be written out to a directory
    titled "seed_cluster_members". By providing this as `pathToSeedDir`, the
    seeds will be loaded into the session.

    Parameters
    ----------
    var1: str
        Path to the PDB file of the structure
    var2: str
        Path to the CGO arrows file

    """
    pass

    # clean up from previous call
    if cleanPrevious:
        cmd.delete('points')
        cmd.delete('arrows')

    hradius = 0.2
    radius = 0.1

    # each line is formatted as CHAINID+RESNUM\tCa_pos\tu_pos\tb_pos\tn_pos
    with open(pathToFile,'r') as file:
        res_lines = list(filter(lambda line: line != '', [line.rstrip('\n') for line in file]))

    print(len(res_lines))
    for line in res_lines:
        data = line.split('\t')
        if (len(data) != 5):
            raise ValueError('Wrong number of fields in '+pathToFile+": "+str(len(data)))

        #get residue properties
        chainID = str(data[0][0])
        resNum = int(data[0][1:])
        resName = chainID+str(resNum)

        #get atom coordinates
        ca_coord= '['+','.join(data[1].split())+']'
        u_coord = '['+','.join(data[2].split())+']'
        b_coord = '['+','.join(data[3].split())+']'
        n_coord = '['+','.join(data[4].split())+']'

        #make psuedo atoms
        caName = resName+"_Ca_point"
        uName = resName+"_u_point"
        bName = resName+"_b_point"
        nName = resName+"_n_point"

        cmd.pseudoatom(object=caName,label='',pos=ca_coord)
        cmd.pseudoatom(object=uName,label='u',pos=u_coord)
        cmd.pseudoatom(object=bName,label='b',pos=b_coord)
        cmd.pseudoatom(object=nName,label='n',pos=n_coord)

        #connect with CGO arrows
        uVecName = resName+'_u_vec'
        bVecName = resName+'_b_vec'
        nVecName = resName+'_n_vec'
        cgo_arrow(atom1=caName, atom2=uName, radius=radius, gap=0.0, hlength=0.25, hradius=hradius, color='white', name=uVecName)
        cgo_arrow(atom1=caName, atom2=bName, radius=radius, gap=0.0, hlength=0.25, hradius=hradius, color='white', name=bVecName)
        cgo_arrow(atom1=caName, atom2=nName, radius=radius, gap=0.0, hlength=0.25, hradius=hradius, color='white', name=nVecName)

    cmd.hide('nonbonded',"*_Ca_point")
    cmd.group('points','*_point')
    cmd.group('arrows','*_vec')
    print('done')
    return
