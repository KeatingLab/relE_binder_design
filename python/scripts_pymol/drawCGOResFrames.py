from pymol import cmd, cgo, CmdException

'''

Method for loading the residue frames and drawing them as CGO arrows

'''

def drawCGOResFrames(pathToFile,cleanPrevious=True):
    r"""Loads the residue frame file and then draws the basis vectors

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

    # each line is formatted as res_name\tCa_pos\tx_pos\ty_pos\tz_pos
    with open(pathToFile,'r') as file:
        res_lines = list(filter(lambda line: line != '', [line.rstrip('\n') for line in file]))

    print(len(res_lines))
    for line in res_lines:
        data = line.split('\t')
        if (len(data) != 5):
            raise ValueError('Wrong number of fields in '+pathToFile+": "+str(len(data)))

        #get residue properties
        # chainID = str(data[0][0])
        # resNum = int(data[0][1:])
        # resName = chainID+str(resNum)
        resName = data[0]

        #get atom coordinates
        ca_coord= '['+','.join(data[1].split())+']'
        x_coord = '['+','.join(data[2].split())+']'
        y_coord = '['+','.join(data[3].split())+']'
        z_coord = '['+','.join(data[4].split())+']'

        #make psuedo atoms
        caName = resName+"_Ca_point"
        xName = resName+"_x_point"
        yName = resName+"_y_point"
        zName = resName+"_z_point"

        cmd.pseudoatom(object=caName,label='',pos=ca_coord)
        cmd.pseudoatom(object=xName,label='x',pos=x_coord)
        cmd.pseudoatom(object=yName,label='y',pos=y_coord)
        cmd.pseudoatom(object=zName,label='z',pos=z_coord)

        #connect with CGO arrows
        xVecName = resName+'_x_vec'
        yVecName = resName+'_y_vec'
        zVecName = resName+'_z_vec'
        cgo_arrow(atom1=caName, atom2=xName, radius=radius, gap=0.0, hlength=0.25, hradius=hradius, color='white', name=xVecName)
        cgo_arrow(atom1=caName, atom2=yName, radius=radius, gap=0.0, hlength=0.25, hradius=hradius, color='white', name=yVecName)
        cgo_arrow(atom1=caName, atom2=zName, radius=radius, gap=0.0, hlength=0.25, hradius=hradius, color='white', name=zVecName)

    cmd.hide('nonbonded',"*_Ca_point")
    cmd.group('points','*_point')
    cmd.group('arrows','*_vec')
    print('done')
    return
