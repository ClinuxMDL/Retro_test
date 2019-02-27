
import json,sys, os, argparse
from rdkit import Chem
# import copy

scriptname="s1_token_process.py"
work_data_dir="data"
data_dir="data_token_all"
train_dir="A2_token_demo"
problem="my_reaction_token"
pythondir="/lustre1/lhlai_pkuhpc/wangsw/software/anaconda3/bin"


def cano(smiles):  # canonicalize smiles by MolToSmiles function
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles)) if (smiles != '') else ''

def process_txt_to_tokenfile(_tmpfile, iternum):
    cmdstr = "cd %s; %s/python %s %s %s test_no%s" % (work_data_dir, pythondir, scriptname, _tmpfile, data_dir, str(iternum))
    os.system(cmdstr)
    return "test_no%s_sources" % (iternum)

def decoder_txt_to_smi(_tmpfile, iternum):
    tokenfile=process_txt_to_tokenfile(_tmpfile, iternum)
    outfile= "output_cano%s.txt" %(str(iternum))
    cmdstr="export PATH=$PATH:%s; bash data_decoder.sh %s %s %s/%s %s %s %d" % (pythondir, train_dir, problem, work_data_dir, data_dir, tokenfile, outfile, iternum*2+4)
    os.system(cmdstr)
    return "%s/train/%s" % (train_dir, outfile)



def filter_smidict_to_file(_txtdict, _curdecdict, _tmpfile):

    uncan_dict={}
    for key in sorted(_curdecdict.keys()):
        tmpsmi=_curdecdict[key]
        try:
            tmpsmi = cano(tmpsmi)
            Decdict[key] = tmpsmi
        except:
            uncan_dict[key] = tmpsmi
            print(key, tmpsmi)
            continue

    # generate tmp txt with error smiles
    print("Number of error smiles: %d" %(len(uncan_dict.keys())))
    if len(uncan_dict.keys())>0:
        with open(os.path.join(work_data_dir, _tmpfile), "w") as f:
            f.write("id,reagents>production\n")
            for key in sorted(uncan_dict.keys()):
                f.write(key + "," + _txtdict[key] + "\n")
    else:
        print("No wrong smiles strings.")
        sys.exit(0)

    return uncan_dict

def add_id_to_smidict(_tmpfile, _smifile):
    f1=open(_tmpfile)
    f2=open(_smifile)

    lines1=f1.readlines()[1:]
    lines2=f2.readlines()
    f1.close()
    f2.close()

    _curdecdict={}
    _txtdict={}

    for line1, line2 in zip(lines1,lines2):
        key,smi=line1.split(",", 1)
        _curdecdict[key]= line2.replace(" ","").replace("<EOS>", "").strip()
        _txtdict[key]= smi.strip()

    return _curdecdict, _txtdict



parser = argparse.ArgumentParser()
parser.add_argument('--decoder_file', type=str,default="sub1_t2t_token/train/a.txt", help='')
parser.add_argument('--raw_txt_file', type=str,default="../data/test_sub.txt", help='')
# parser.add_argument('--process_file', type=str,default="error_token1.json", help='')

parser.add_argument('--out1_file', type=str,default="output_token1_r1.json", help='')
parser.add_argument('--out2_file', type=str,default="error_token1_r1.json", help='')
parser.add_argument('--count', type=int,default=2, help='')


opt = parser.parse_args()

smifile=opt.decoder_file
# file2=opt.process_file

txtfile=opt.raw_txt_file

out1_file=opt.out1_file
out2_file=opt.out2_file
# out2_file=sys.argv[4]
count=opt.count
#
# #record raw txt into dict
# f1 = open(txtfile)
# TXTdict = {}
# for line in f1.readlines():
#     key, smi = line.split(",")
#     TXTdict[key] = smi.strip()
# f1.close()
#
#
# # record decoder smiles into dict
# f2=open(file1)
# Decdict={}
# for ix, content in enumerate(f2.readlines()):
#     ixstr='%06d' % ix
#     Decdict[ixstr]=str(content).replace(" ","").replace("<EOS>", "").strip()
# f2.close()

Decdict, TXTdict=add_id_to_smidict(txtfile,smifile)


import copy


curdict=copy.copy(Decdict)

iternum=2
for iter in range(1,iternum):
    uncandict1=filter_smidict_to_file(TXTdict, curdict, "test_no%d.txt" %(iter))
    print(uncandict1)
    outsmifile= decoder_txt_to_smi("test_no%d.txt" %(iter), iter)
    curdict, _ = add_id_to_smidict("%s/test_no%d.txt" %(work_data_dir, iter), outsmifile)


uncandict2=filter_smidict_to_file(TXTdict, curdict, "test_no%d.txt" %(count))
print(uncandict1==uncandict2)

with open(out1_file,"w") as f:
    json.dump(Decdict, f, ensure_ascii=False, sort_keys=True)

# print("Number of error smiles: %d" %(len(uncandict2.keys())))
with open(out2_file,"w") as f:
    json.dump(uncandict2, f, ensure_ascii=False, sort_keys=True)


#
# # process wrong smiles
# rawdict=json.load(open(file1))
# # curdict=json.load(open(file2))
#
#
# tmpfile="test_no"+count+".txt"
#
# with open("data/"+tmpfile, "w") as f:
#     f.write("id,reagents>production\n")
#     for key in sorted(curdict.keys):
#         f.write(key+","+ txtdict[key] +"\n")
#
# scriptname="s1_token_process.py" if mode =="token" else "s1_char_process.py"
# data_dir= "data_token_all" if mode=="token" else "data_char_all"
#
#
#
# cmdstr1="cd data; python %s %s %s test_no%s" % (scriptname, tmpfile, data_dir, count)
# os.system(cmdstr1)
#
# cmdstr2= "bash data_decoder_tiny.sh A1_subword_demo_tiny my_reaction_token data/%s test_no%s_sources output_token_no%s.txt" % (data_dir, count, count)
# os.system(cmdstr2)
#
# decoderfile="A1_subword_demo_tiny/train/output_token_no%s.txt" %(count)
#
# if os.path.exists(decoderfile) :
#     with open(decoderfile) as f:
#         for line in f.readlines():
#             str(content).replace(" ","").replace("<EOS>", "").strip()
# else:
#     print("no ouptut decoder file")
#     sys.exit(0)
#
#
#
# corr_dict={}
# uncan_dict={}
# for k in sorted(curdict.keys()):
#     tmpsmi= curdict[k]
#     try:
#         tmpsmi=cano(tmpsmi)
#         corr_dict[k]=tmpsmi
#     except:
#         uncan_dict[k]=tmpsmi
#         print(k,tmpsmi)
#         # traceback.print_exc()
#         continue
#
# # modify raw output
#
# for k in sorted(corr_dict.keys()):
#     rawdict[k]=corr_dict[k]
#
# with open(out1_file,"w") as f:
#     json.dump(rawdict, f, ensure_ascii=False, sort_keys=True)
#
# with open(out2_file,"w") as f:
#     json.dump(uncan_dict, f, ensure_ascii=False, sort_keys=True)
# print(len(uncan_dict.keys()))
# print("="*50)