from pymol import cmd

def colorByBfactor(filepath: str,obj_name:str,value_col='score'):
    value_list = parseOutputCSV(filepath,value_col)
    setBFactors(obj_name,value_list)


def parseOutputCSV(filepath: str, value_col: str, chainid_col: str = 'chain_id', resnum_col: str = 'resnum'):
    value_list = [] # [(chainid,resnum),value]
    with open(filepath) as file:
        col_num = -1
        for i,line in enumerate(file):
            split = line.rstrip().split(',')
            if i == 0:
                col_num = len(split)
                name2colIdx = {name:i for i,name in enumerate(split)}
            else:
                assert len(split) == col_num
                # print([(split[name2colIdx[chainid_col]],split[name2colIdx[resnum_col]]),split[name2colIdx[value_col]]])
                value_list.append(((split[name2colIdx[chainid_col]],split[name2colIdx[resnum_col]]),float(split[name2colIdx[value_col]])))
    return value_list

def setBFactors(obj_name:str, value_list: list):
    # alter the b factors
    for (chainid,resnum),value in value_list:
        sel = f"(chain {chainid} and resid {resnum}) and {obj_name}"
        val = f"b={value}"
        cmd.alter(sel,val)
    sel = f"({' or '.join(['chain '+chainid+' and resid '+resnum for (chainid,resnum),_ in value_list])}) and " + obj_name
    print('sel',sel)
    # reset the color values
    cmd.spectrum('b','blue_red',sel,minimum=-2,maximum=2)
