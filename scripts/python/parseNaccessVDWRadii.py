'''
Loads the vdw.radii file from Naccess, gathers the radii, and writes them out
in format that is easy to copy.
'''

import regex as re
import pandas as pd

path = '/Users/sebastianswanson/Keating/utilities/repos/naccess/Naccess/vdw.radii'
mapName = 'radii'

# load file
data = list()
with open(path,'r') as file:
    for line in file:
        line  = line.rstrip()
        if len(line) > 0 and line[0] != "#":
            data.append(line.rstrip())

# grab data from lines
atomType = ""
resName = ""
atomName = ""
vdwRadius = 0.0

vdwRadii = list() # [atomType,resName,atomName,vdwRadius]

resPat = '^(RESIDUE)\s(HETATM|ATOM|NUCL|WATER)\s(.{3})'
atomPat = '^(ATOM)\s(.{4})\s(\d.\d{2})'

for line in data:
    resMatch = re.search(resPat,line)
    atomMatch = re.search(atomPat,line)
    if (resMatch != None):
        atomType = resMatch[2]
        resName = resMatch[3]
    elif (atomMatch != None):
        atomName = atomMatch[2].strip()
        vdwRadius = atomMatch[3]
        # print(atomType,resName,atomName,vdwRadius)
        vdwRadii.append([atomType,resName,atomName,vdwRadius])
    else:
        print("Error: line not recognized:",line)

df = pd.DataFrame(vdwRadii,columns=['atomType','resName','atomName','vdwRadius'])
df.to_csv('vdwRadii.csv')

# write in a format that can be copied to cpp file
with open('vdwRadiiForCPP.txt','w') as file:
    for i,row in df.iterrows():
        file.write(mapName+'["'+'"]["'.join([row['resName'],row['atomName']])+'"] = '+row['vdwRadius']+';\n')
