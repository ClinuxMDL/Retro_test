# coding: utf-8
'''
python s1_token_process.py inputfile.txt outputdir mode
'''
from __future__ import print_function
import os
import re
import traceback
import sys
import pickle, gzip
from rdkit.Chem import AllChem
from rdkit import Chem
import copy
import numpy as np


# import parser.Smipar as Smipar

def cano(smiles):  # canonicalize smiles by MolToSmiles function
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if (smiles != '') else ''


#
inputfile = sys.argv[1]
output_path = sys.argv[2]
mode= sys.argv[3]
onepath = output_path

if not os.path.exists(onepath):
    os.makedirs(onepath)

fp_train = open(inputfile)
fp_train_sources = open(os.path.join(onepath, mode+"_sources"), "w")
fp_train_targets = open(os.path.join(onepath, mode+"_targets"), "w")

length_list = []
vocab = {}

fp_train.readline()


def process_trainset():
    for line in fp_train:
        line=line.replace("\"","")
        if " |" in line:
            tmp_fields = line.strip().split(" |")
            line = tmp_fields[0]
        fields = line.strip().split(",")
        id = fields[0]
        rxn_str = fields[1]
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_str, useSmiles=True)
        except:
            print(id,rxn_str)
            # print(id)
            traceback.print_exc()
            continue
        AllChem.RemoveMappingNumbersFromReactions(rxn)  # 去原子的号码
        output_smiles = AllChem.ReactionToSmiles(rxn)
        length_list.append((output_smiles.rfind('>'), len(output_smiles) - output_smiles.rfind('>') - 1))

        reactant_list = []
        agent_list = []
        product_list = []
        split_rsmi = output_smiles.split('>')
        try:
            reactants = cano(split_rsmi[0])
        except:
            print(id, split_rsmi[0])
            traceback.print_exc()
            continue
        try:
            agents = cano(split_rsmi[1])
        except:
            print(id, split_rsmi[1])
            traceback.print_exc()
            continue
        try:
            products = cano(split_rsmi[2])
        except:
            print(id, split_rsmi[2])
            traceback.print_exc()
            continue
        token_regex = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|= |  # |-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"

        def split_smi_to_token(smistr):
            tmp_list=[]
            for smi in smistr.split("."):
                tmp_list += re.split(token_regex, smi)
                tmp_list += '.'
            tmp_list.pop()
            return tmp_list

        reactant_list=split_smi_to_token(reactants)
        agent_list=split_smi_to_token(agents)
        product_list=split_smi_to_token(products)


        def randomize_smiles(smiles, can=False):
            """Perform a randomization of a SMILES string
            must be RDKit sanitizable"""
            m = Chem.MolFromSmiles(smiles)
            ans = list(range(m.GetNumAtoms()))
            np.random.shuffle(ans)
            nm = Chem.RenumberAtoms(m, ans)
            return Chem.MolToSmiles(nm, canonical=can, isomericSmiles=True)

        reactant_list1 = split_smi_to_token(randomize_smiles(reactants))
        agent_list1 = split_smi_to_token(randomize_smiles(agents))
        product_list1 = split_smi_to_token(randomize_smiles(products))



        #
        # for reactant in reactants:
        #     # print(re.split(token_regex, reactant))
        #     # reactant_list += Smipar.parser_list(reactant)
        #     reactant_list += re.split(token_regex, reactant)
        #     reactant_list += '.'
        # for agent in agents:
        #     # agent_list += Smipar.parser_list(agent)
        #     agent_list += re.split(token_regex, agent)
        #     agent_list += '.'
        # for product in products:
        #     # product_list += Smipar.parser_list(product)
        #     product_list += re.split(token_regex, product)
        #     product_list += '.'
        #
        # reactant_list.pop()  # to pop last '.'
        # agent_list.pop()
        # product_list.pop()

        # reactant_list += '>'
        # reactant_list += agent_list

        agent_list += '>'
        agent_list += product_list

        agent_list1 += '>'
        agent_list1 += product_list1

        #argumentation





        product_list = agent_list
        product_list1 = agent_list1

        def strip_token(tmplist):
            for token_i in range(len(tmplist)):
                tmplist[token_i] = tmplist[token_i].strip()

            tmplist1= [i for i in tmplist if i != ""]
            return tmplist1

        out_product_list=strip_token(product_list)
        out_product_list1=strip_token(product_list1)

        out_reactant_list= strip_token(reactant_list)
        out_reactant_list1= strip_token(reactant_list1)


        # for token_i in range(len(product_list)):
        #     product_list[token_i] = product_list[token_i].strip()
        #     product_list1[token_i] = product_list1[token_i].strip()
        #
        # for token_i in range(len(reactant_list)):
        #     reactant_list[token_i] = reactant_list[token_i].strip()
        #     reactant_list1[token_i] = reactant_list1[token_i].strip()



        # out_product_list = []
        # out_reactant_list = []
        # out_product_list1 = []
        # out_reactant_list1 = []
        # out_product_list=strip_token()

        # for token in product_list:
        #     strip_token = token.strip()
        #     if strip_token != "":
        #         out_product_list.append(strip_token)
        # for token in reactant_list:
        #     strip_token = token.strip()
        #     if strip_token != "":
        #         out_reactant_list.append(strip_token)
        print(" ".join(out_product_list), file=fp_train_sources)
        print(" ".join(out_product_list1), file=fp_train_sources)
        print(" ".join(out_reactant_list), file=fp_train_targets)
        print(" ".join(out_reactant_list1), file=fp_train_targets)

    return 1


def process_testset():
    for line in fp_train:
        line = line.replace("\"", "")
        if " |" in line:
            tmp_fields = line.strip().split(" |")
            line = tmp_fields[0]
        fields = line.strip().split(",")
        id = fields[0]
        rxn_str = fields[1].replace(">", ">>")
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_str, useSmiles=True)
        except:
            print(id, rxn_str)
            # print(id)
            traceback.print_exc()
            continue
        AllChem.RemoveMappingNumbersFromReactions(rxn)  # 去原子的号码
        output_smiles = AllChem.ReactionToSmiles(rxn)

        agent_list = []
        product_list = []
        split_rsmi = output_smiles.split('>')

        try:
            agents = cano(split_rsmi[0]).split('.')
        except:
            print(id, split_rsmi[0])
            # print(id)
            traceback.print_exc()
            continue
        try:
            products = cano(split_rsmi[2]).split('.')
        except:
            print(id,split_rsmi[2])
            # print(id)
            traceback.print_exc()
            continue
        token_regex = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|= |  # |-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"

        for agent in agents:
            # agent_list += Smipar.parser_list(agent)
            agent_list += re.split(token_regex, agent)
            agent_list += '.'
        for product in products:
            # product_list += Smipar.parser_list(product)
            product_list += re.split(token_regex, product)
            product_list += '.'

        # reactant_list.pop()  # to pop last '.'
        agent_list.pop()
        product_list.pop()

        # reactant_list += '>'
        # reactant_list += agent_list

        agent_list += '>'
        agent_list += product_list

        product_list = agent_list

        tmp_product_list = copy.copy(product_list)
        # tmp_reactant_list = copy.copy(reactant_list)
        for token_i in range(len(product_list)):
            product_list[token_i] = product_list[token_i].strip()

        out_product_list = []
        # out_reactant_list=[]
        for token in product_list:
            strip_token = token.strip()
            if strip_token != "":
                out_product_list.append(strip_token)

        print(" ".join(out_product_list), file=fp_train_sources)
    return 1


if mode =="train":
    process_trainset()
else:
    process_testset()