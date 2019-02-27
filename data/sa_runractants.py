import json,sys
from rdkit.Chem import AllChem
from rdkit import Chem

# import matplotlib.pyplot as plt




def cano(smiles):  # canonicalize smiles by MolToSmiles function
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if (smiles != '') else ''

infile1="react_template/count_reacts.txt"
infile2="test_all_r1.txt"
step=5000
start_ix= int(0)*step
end_ix= start_ix+step if start_ix+step < 238280 else -1
out1file= "tempJson_{}_{}.json".format(start_ix,end_ix)

rxn_templates=[]
with open(infile1) as f:
    for line in f.readlines():
        count, retro_smarts= line.strip().split(" ")
        if int(count)>=3:
            rxn_templates.append(retro_smarts)


f=open(infile2)
alllines=f.readlines()[1:2]
f.close()

templateDict={}

try:

    for line in alllines:
        line = line.replace("\"", "")
        if " |" in line:
            tmp_fields = line.strip().split(" |")
            line = tmp_fields[0]
        fields = line.strip().split(",")
        id = fields[0]
        rxn_str = fields[1].replace(">", ">>")

        curlist=[]
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_str, useSmiles=True)
        except:
            print(id, rxn_str)
            continue
        AllChem.RemoveMappingNumbersFromReactions(rxn)  # 去原子的号码
        rxn_products = AllChem.ReactionToSmiles(rxn).split(">")[2]
        # rxn_str = fields[2]
        try:
            m1=[Chem.MolFromSmiles(x) for x in rxn_products.split(".")]

            for tpid, rxn in enumerate(rxn_templates):
                # print("template_id: {}, {}".format(i, rxn))
                try:
                    rxn1=AllChem.ReactionFromSmarts(rxn)
                    # print("template_id: {}, {}".format(i, rxn))

                    reactants =rxn1.RunReactants(m1)
                    for i in range(len(reactants)):
                        list1=[Chem.MolToSmiles(reactants[i][j],1) for j in range(len(reactants[i]))]
                        # list1 +=[str(tpid)]
                        curlist.append(cano(".".join(list1)))

                    # print(len(reactants))
                    # for m in reactants:
                    #     print(Chem.MolToSmiles(reactants[m], 1))
                    # tmplist=[(Chem.MolToSmiles(x, 1)) for x in rxn1.RunReactants(m1)[0]]
                    #
                    # templateDict[id].append(cano(".".join(tmplist)))
                    # print("==="*30)
                    # break
                except:
                    # print("This is wrong templates!")
                    continue
        except:

            continue
        if len(curlist)>0:
            templateDict[id]=list(set(curlist))

except KeyboardInterrupt:
    with open(out1file, "w") as f:
        # json.dump(templateDict, )
        json.dump(templateDict, f, ensure_ascii=False, sort_keys=True)




with open(out1file, "w") as f:
    # json.dump(templateDict, )
    json.dump(templateDict, f, ensure_ascii=False, sort_keys=True)