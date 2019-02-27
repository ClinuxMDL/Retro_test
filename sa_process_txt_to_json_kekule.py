from __future__ import print_function
import sys,json
import os
import re
import traceback
import sys
import pickle, gzip
from rdkit.Chem import AllChem
from rdkit import Chem
import copy

def cano(smiles):  # canonicalize smiles by MolToSmiles function
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), kekuleSmiles=True) if (smiles != '') else ''

def unique_str(smiles):
    return ".".join(list(set(smiles.split("."))))


# in_list=[3379, 24175, 29916, 29993, 39520, 45418, 46657, 49009, 107311,153061,171366,185510, 192990]

in1_file=sys.argv[1]
out1_file=sys.argv[2]
out2_file=sys.argv[3]

iscanon=sys.argv[4]


f1=open(in1_file)

needsmi=f1.readlines()
curdict={}
for ix, content in enumerate(needsmi):
    ixstr='%06d' % ix
    curdict[ixstr]=str(content).replace(" ","").replace("<EOS>", "").strip()


if int(iscanon)>0:
    uncan_dict={}
    for k in sorted(curdict.keys()):
        tmpsmi= unique_str(curdict[k])


        try:
            m1 = Chem.MolFromSmiles(tmpsmi)
            Chem.Kekulize(m1)
            tmpsmi = Chem.MolToSmiles(m1, kekuleSmiles=True)
            curdict[k]=tmpsmi
        except:

            uncan_dict[k]=tmpsmi
            print(k,tmpsmi)
            # traceback.print_exc()
            continue
    print("Number of error smiles: %d" % (len(uncan_dict.keys())))
    with open(out2_file, "w") as f:
        json.dump(uncan_dict, f, ensure_ascii=False, sort_keys=True)


print(len(curdict.keys()))
with open(out1_file,"w") as f:
    json.dump(curdict, f, ensure_ascii=False, sort_keys=True)






