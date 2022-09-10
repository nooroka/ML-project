#ошибка: ids не такие же, как listsp3
#проверка на THASP
#в record_dict THASP двое
#ошибочно одна последовательность присваивается  AZOSB: изначально их 3, но стало 4
from Bio import SeqIO
import pprint
import numpy as np
from Bio import AlignIO
import os
from  Bio.Seq import Seq
from collections import defaultdict
filelist = []
for root, dirs, files in os.walk('/mnt/work/nooroka/alignments300'):
    for file1 in files:
        filelist.append(os.path.join(root,file1))
def get_key(d, value,valuedict):
    for k, v in d.items():
        for s in range(len(v)):
            if v[s] == value:
                val = valuedict[value]
                print("val "+str(val)+"value "+str(value))
                if len(val) ==1:
                    val2 = val[0]
                    return val2
                else:
                    val2 = val
                    return val2


           
#print("filelist " +str(filelist))
for name in filelist:
    #name = "/mnt/work/nooroka/alignments300/Proteobacteria-60-balanced/RAxML/PB_PF12327_1.fasta.fasta"
    #record_dict = {}
    #name = "PB_PF00958_1.fasta.fasta"
    print(name)
    record_dict = defaultdict(list)
    handle = open(name)
    valuedict = defaultdict(list)
    for record in AlignIO.read(handle, "fasta") :
        record_dict[record.id].append(str(record.seq))
        valuedict[record.seq].append(str(record.id))
    handle.close()
    print("record_dict "+str(record_dict))
    print("reckeys "+str(len(record_dict.keys())))
    name2 = name.split("/")

    os.chdir("/mnt/work/nooroka/one-hotnew")
    if (not os.path.exists("/mnt/work/nooroka/one-hotnew/{}/{}".format(name2[5],name2[6]))) and (not os.path.exists("/mnt/work/nooroka/one-hotnew/{}".format(name2[5]))):
          os.mkdir("/mnt/work/nooroka/one-hotnew/{}".format(name2[5]))
          os.mkdir("/mnt/work/nooroka/one-hotnew/{}/{}".format(name2[5],name2[6]))
    elif (not os.path.exists("/mnt/work/nooroka/one-hotnew/{}/{}".format(name2[5],name2[6]))) and  (os.path.exists("/mnt/work/nooroka/one-hotnew/{}".format(name2[5]))):
          os.mkdir("/mnt/work/nooroka/one-hotnew/{}/{}".format(name2[5],name2[6]))
    os.chdir("/mnt/work/nooroka/one-hotnew/{}/{}".format(name2[5],name2[6]))
    # Join all the sentences together and extract the unique characters
    # from the combined sentences
    #chars = set(''.join(list(record_dict.values())))
    chars = {"G","L","Y","S","E","Q","D","N","F","A","K","R","H","C","V","P","W","I","M","T","-"}
    # Creating a dictionary that maps integers to the characters
    char2int = {'-':0, 'A':1,'C':2,'D':3,'E':4,'F':5, 'G':6, 'H':7, 'I':8,'K':9,'L':10,'M':11, 'N':12, 'P':13, 'Q':14, 'R':15,'S':16, 'T':17, 'V':18, 'W':19, 'Y':20,'X':21}  
    textold = list(record_dict.values()) 
    batch_size = 0
    for d in range(len(textold)):
        for k in range(len(textold[d])):
            batch_size+=1
    text = [item for sublist in textold for item in sublist]
    maxlen = 300
    # Creating lists that will hold our input and target sequences
    input_seq = []
    ids = []

    for i in range(len(text)):
        # Remove last character for input sequence
        input_seq.append(text[i])
        print(get_key(record_dict,text[i],valuedict))
        ids.append(get_key(record_dict,text[i],valuedict))
        #print("Input sequence:".ljust(18), f"'{input_seq[i]}'")
    #print("len ids "+str(len(ids)))
    #print("input seq "+str(input_seq))       
    print("ids "+str(ids))
    se = 0
    ids2 = []
    idcorr = dict()
    se = 0
    for h in range(len(ids)):
        print(ids[h])
        if len(ids[h]) == 5:
            ids2.append(ids[h])
        elif len(ids[h])!=5:
            if (str(ids[h]) in idcorr.keys() and h>= 0):
                se = idcorr.get(str(ids[h]))
                ids2.append(str(ids[h][se]))
                se+=1
                if se <= len(ids[h])-1:
                    idcorr[str(ids[h])] = se
                else:
                    se = 0
                    idcorr[str(ids[h])] = se
            if str(ids[h]) not in idcorr.keys():
                print("se "+str(se))
                se = 0
                ids2.append(str(ids[h][se]))
                se +=1
                idcorr[str(ids[h])] = se

    print("ids2 "+str(ids2))
    for i in range(len(text)):
        input_seq[i] = [char2int[character] for character in input_seq[i]]
        #print("Encodded input sequence:".ljust(25), input_seq[i])
        print()

    dict_size = len(char2int)
    seq_len = maxlen - 1
#print("dict "+str(dict_size))
#print("seq "+str(seq_len))
#print("batch "+str(batch_size))
    def one_hot_encode(sequence, dict_size, seq_len, batch_size):
        # Creating a multi-dimensional array of zeros with the desired output shape
        features = np.zeros((batch_size, seq_len, dict_size), dtype=np.float32)
        print(features.shape)
        # Replacing the 0 at the relevant character index with a 1 to represent that character
        for i in range(batch_size):
            for u in range(seq_len):
                features[i, u, sequence[i][u]] = 1
        return features
    input_seq2 = one_hot_encode(input_seq, dict_size, seq_len, batch_size)
    #print("Input shape: {} --> (Batch Size, Sequence Length, One-Hot Encoding Size)".format(input_seq2.shape))
    #print("input_seq2 "+str(input_seq2))
    d = 0
    d2 = 0
    '''
    w2 = open("/data/nooroka/ML/test3_2.txt".format(name2[-1]),"w")
    for i in range(len(input_seq2)):
        d2+=1
        w2.write(">"+str(ids[i])+"\n")
        for xs2 in input_seq2[i]:
            w2.write(" ".join(map(str, xs2))+" "+str(input_seq[i])+"\n")
        w2.write("\n")
        if d2%4 == 0:
            w2.write("\n")
            w2.write("\n")
            w2.write("\n")
            w2.write("\n")
    w2.close()
    '''
    w = open("/mnt/work/nooroka/one-hotnew/{}/{}/{}.txt".format(name2[5],name2[6],name2[-1]),"w")
    dict1 = defaultdict(list)
    dict2 = defaultdict(list)
    listsp = []
    name44 = name2[-1].split(".")
    op2 = open("/mnt/work/nooroka/ids150822/{}/{}/{}.tre.txt".format(name2[5],name2[6],name44[0]),"r")
    #op2 = open("/data/nooroka/ML/PB_PF00958_1.tre.txt")
    for line in op2:
        line = line.strip()
        listsp.append(line)
    op2.close()
    listsp2 = []
    listsp3 = []
    listsp4 = []

    for k in listsp:
        if k!='':
            listsp2.append(k)
            listsp4.append(k)
        elif k=='':
            listsp3.append(listsp2)
            listsp2 = []


    for i in range(len(input_seq2)):
        #d+=1
        dict1[ids2[i]].append(input_seq2[i])
    list2 = []
    sc  = 0
    for t in range(len(listsp3)):
        for g in range(len(listsp3[t])):
            a = listsp3[t][g]
            sc+=1
            print("dict1el " +str(listsp3[t][g])+" "+str(dict1[listsp3[t][g]]))
            b = dict1[listsp3[t][g]][0]
            del dict1[listsp3[t][g]][0]
            list2.append((a,b))
        

    for h in range(len(list2)):
        d+=1
        w.write(">"+str(list2[h][0])+"\n")
        for xs in list2[h][1]:
            w.write("".join(str(xs))[1:-1]+" "+"\n")
        w.write("\n")
        if d%4 == 0:
            w.write("\n")
            w.write("\n")
            w.write("\n")
            w.write("\n")
    w.close()

