import blosum as bl
from collections import defaultdict
from Bio import SeqIO
import pprint
import numpy as np
from Bio import AlignIO
import os
from  Bio.Seq import Seq

filelist = []
for root, dirs, files in os.walk('/mnt/work/nooroka/alignments300'):
    for file1 in files:
        filelist.append(os.path.join(root,file1))
mat = bl.BLOSUM(62)    
valueXY = mat["XY"]
alp = "ACDEFGHIKLMNPQRSTVWY"
gap = "-"
def get_key(d, value):
    for k, v in d.items():
        if v == value:
            return k

blosumDict = dict(mat)
for m in range(len(alp)):
    blosumDict[str(alp[m])+str(gap)] = -10
    blosumDict[str(gap)+str(alp[m])] = -10
blosumDict[str(gap)+str(gap)] = -1
print(blosumDict)
enddict = defaultdict(list)
for name in filelist:
    record_dict = {}
    handle = open(name)
    for record in AlignIO.read(handle, "fasta") :
        print(str(record.id)+" "+str(record.seq))
        id1 = str(record.id)
        seq1 = str(record.seq)
        record_dict[str(id1)].append(str(seq1))
    handle.close()

    print(record_dict)
    name2 = name.split("/")
    os.chdir("/mnt/work/nooroka/blosum/")
    if not os.path.exists("/mnt/work/nooroka/blosum/{}".format(name2[5])):
        os.mkdir("{}".format(name2[5]))
    os.chdir("{}".format(name2[5]))
    text = list(record_dict.values())
    maxlen = 300
    ids = []
   # string = "ALTLSPYYK"

    for i in range(len(text)):
        ids.append(get_key(record_dict,text[i]))
        for k  in range(len(text[i])):
            list1 = []
            for j in range(len(alp)):
                a = str(text[i][k])+str(alp[j])
                list1.append(blosumDict[a])
            enddict[text[i][k]] = list1
            #print('text[i][k] '+str(text[i][k])+" enddict[text[i][k]]  "+str(enddict[text[i][k]]))
    w = open("{}".format(name2[-1]),"w")
    d = 0
    print("ids "+str(ids))
    print("text "+str(text))
    for r in range(len(text)):
        d+=1
        w.write(">"+str(ids[r])+"\n")
        for m in range(len(text[r])):
            stri = ''
            for g in range(len(enddict[text[r][m]])):
                stri = stri+" "+str(enddict[text[r][m]][g])
            w.write(stri+" "+str(text[r][m])+"\n")
        w.write("\n")
        if d%4 == 0:
            w.write("\n")
            w.write("\n")
            w.write("\n")
            w.write("\n")
    w.close()
    print(enddict)
