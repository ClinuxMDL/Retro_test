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

token_regex = "(\[[^\]]+]|Br?|Cl?|N|O|S|P|F|I|b|c|n|o|s|p|\(|\)|\.|= |  # |-|\+|\\\\|\/|:|~|@|\?|>|\*|\$|\%[0-9]{2}|[0-9])"

def canoOrder(smiles):  # canonicalize smiles by MolToSmiles function
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), canonical=False) if (smiles != '') else ''


def reorder_smiles(smistr):
    tmp_list = smistr.split(".")
    tmp_list1 = sorted(tmp_list, key=lambda i: len(i), reverse=True)

    return ".".join(tmp_list1)


def strip_token(tmplist):
    for token_i in range(len(tmplist)):
        tmplist[token_i] = tmplist[token_i].strip()

    tmplist1 = [i for i in tmplist if i != ""]
    return tmplist1


def split_smi_to_token(smistr):
    tmp_list = []
    for smi in smistr.split("."):
        tmp_list += re.split(token_regex, smi)
        tmp_list += '.'
    tmp_list.pop()
    return tmp_list

#
inputfile = sys.argv[1]
output_path = sys.argv[2]
mode= sys.argv[3]

# inputfile = "train_sub1_wt.txt"
# output_path = "data_test1"
# mode= "train"


onepath = output_path

if not os.path.exists(onepath):
    os.makedirs(onepath)

fp_train = open(inputfile)
fp_train_sources = open(os.path.join(onepath, mode+"_sources"), "w")
fp_train_targets = open(os.path.join(onepath, mode+"_targets"), "w")


fp_train.readline()
def process_trainset():
    for line in fp_train:
        line=line.replace("\"","")
        if " |" in line:
            tmp_fields = line.strip().split(" |")
            line = tmp_fields[0]
        fields = line.strip().split(",")
        id = fields[0]
        rxn_type = fields[1]
        rxn_str= fields[2]
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_str, useSmiles=True)
        except:
            print(id,rxn_str)
            # traceback.print_exc()
            continue
        AllChem.RemoveMappingNumbersFromReactions(rxn)  # 去原子的号码
        output_smiles = AllChem.ReactionToSmiles(rxn)
        # print(output_smiles)
        split_rsmi = output_smiles.split('>')

        try:
            reactants = canoOrder(split_rsmi[0])
        except:
            print(id, split_rsmi[0])
            # traceback.print_exc()
            continue
        try:
            agents = canoOrder(split_rsmi[1])
        except:
            print(id, split_rsmi[1])
            # traceback.print_exc()
            continue
        try:
            products = canoOrder(split_rsmi[2])
        except:
            print(id, split_rsmi[2])
            # traceback.print_exc()
            continue



        # reactants=reorder_smiles(reactants)
        # agents=reorder_smiles(agents)
        # products=reorder_smiles(products)
        # print(reactants+">"+agents+">"+products)



        reactant_list=split_smi_to_token(reactants)
        agent_list=split_smi_to_token(agents)
        product_list=split_smi_to_token(products)

        product_list.insert(0, rxn_type)
        product_list += '>'
        product_list += agent_list


        out_product_list=strip_token(product_list)
        # out_product_list1=strip_token(product_list1)

        out_reactant_list= strip_token(reactant_list)
        # out_reactant_list1= strip_token(reactant_list1)

        # print("".join(out_product_list))

        print(" ".join(out_product_list), file=fp_train_sources)
        print(" ".join(out_reactant_list), file=fp_train_targets)

    return 1


def process_testset():
    for line in fp_train:
        line = line.replace("\"", "")
        if " |" in line:
            tmp_fields = line.strip().split(" |")
            line = tmp_fields[0]
        fields = line.strip().split(",")
        id = fields[0]
        rxn_type=fields[1]
        rxn_str = fields[2].replace(">", ">>")
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_str, useSmiles=True)
        except:
            print(id, rxn_str)
            # print(id)
            # traceback.print_exc()
            continue
        AllChem.RemoveMappingNumbersFromReactions(rxn)  # 去原子的号码
        output_smiles = AllChem.ReactionToSmiles(rxn)

        split_rsmi = output_smiles.split('>')

        try:
            agents = canoOrder(split_rsmi[0])
        except:
            print(id, split_rsmi[0])
            # print(id)
            # traceback.print_exc()
            continue
        try:
            products = canoOrder(split_rsmi[2])
        except:
            print(id,split_rsmi[2])
            # print(id)
            # traceback.print_exc()
            continue

        # reactant_list = split_smi_to_token(reactants)
        agent_list = split_smi_to_token(agents)
        product_list = split_smi_to_token(products)


        # reactant_list += '>'
        # reactant_list += agent_list
        product_list.insert(0,rxn_type)
        product_list += '>'
        product_list += agent_list

        # product_list = agent_list
        out_product_list = strip_token(product_list)

        print(" ".join(out_product_list), file=fp_train_sources)
    return 1


if mode =="train":
    process_trainset()
else:
    process_testset()