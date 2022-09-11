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
print("filelist " +str(filelist))
for name in filelist:
    #name = "/mnt/work/nooroka/alignments300/Proteobacteria-60-balanced/RAxML/PB_PF12327_1.fasta.fasta"
    #record_dict = {}
#name = "PB_PF18073_1.fasta.fasta"
    print(name)
    record_dict = defaultdict(list)
    handle = open(name)
    valuedict = defaultdict(list)
    text = []
    for record in AlignIO.read(handle, "fasta") :
        record_dict[record.id].append(str(record.seq))
        valuedict[record.seq].append(str(record.id))
        text.append((record.id,record.seq))
    handle.close()
    print("record_dict "+str(record_dict))
    print("reckeys "+str(len(record_dict.keys())))
    name2 = name.split("/")
    os.chdir("/mnt/work/nooroka/one-hotnew2npy")
    if (not os.path.exists("/mnt/work/nooroka/one-hotnew2npy/{}/{}".format(name2[5],name2[6]))) and (not os.path.exists("/mnt/work/nooroka/one-hotnew2npy/{}".format(name2[5]))):
          os.mkdir("/mnt/work/nooroka/one-hotnew2npy/{}".format(name2[5]))
          os.mkdir("/mnt/work/nooroka/one-hotnew2npy/{}/{}".format(name2[5],name2[6]))
    elif (not os.path.exists("/mnt/work/nooroka/one-hotnew2npy/{}/{}".format(name2[5],name2[6]))) and  (os.path.exists("/mnt/work/nooroka/one-hotnew2npy/{}".format(name2[5]))):
          os.mkdir("/mnt/work/nooroka/one-hotnew2npy/{}/{}".format(name2[5],name2[6]))
    os.chdir("/mnt/work/nooroka/one-hotnew2npy/{}/{}".format(name2[5],name2[6]))
    
    # Join all the sentences together and extract the unique characters
    # from the combined sentences
    #chars = set(''.join(list(record_dict.values())))
    chars = {"G","L","Y","S","E","Q","D","N","F","A","K","R","H","C","V","P","W","I","M","T","-"}
    # Creating a dictionary that maps integers to the characters
    char2int = {'-':0, 'A':1,'C':2,'D':3,'E':4,'F':5, 'G':6, 'H':7, 'I':8,'K':9,'L':10,'M':11, 'N':12, 'P':13, 'Q':14, 'R':15,'S':16, 'T':17, 'V':18, 'W':19, 'Y':20,'X':21}  
    #for key in record_dict:
    #    for i in range(len(record_dict[key])):
    #       text.append((key,record_dict[key][i]))
    maxlen = 301
    # Creating lists that will hold our input and target sequences
    input_seq = []
    ids = []
    print("text "+str(text))
    text2 = []
    for i in range(len(text)):
        text2.append((text[i][0],[char2int[character] for character in text[i][1]]))
    print("text2 "+str(text2))
    dict_size = len(char2int)
    seq_len = maxlen - 1
    batch_size = len(text)
    def one_hot_encode(sequence, dict_size, seq_len, batch_size):
        # Creating a multi-dimensional array of zeros with the desired output shape
        features = np.zeros((batch_size, seq_len, dict_size), dtype=np.float32)
        print(features.shape)
        # Replacing the 0 at the relevant character index with a 1 to represent that character
        for i in range(batch_size):
            for u in range(seq_len):
                features[i, u, sequence[i][u]] = 1
        return features
    text3 = []
    input_seq = []
    for k in range(len(text2)):
        input_seq.append(text2[k][1])

    input_seq2 = one_hot_encode(input_seq, dict_size, seq_len, batch_size)
    print(input_seq2)
    for j in range(len(text2)):
        text3.append((text2[j][0],input_seq2[j]))
    print("text3 "+str(text3))
    def func_chunk(lst, n):
        for x in range(0, len(lst), n):
            e_c = lst[x : n + x]

            if len(e_c) < n:
                e_c = e_c + [None for y in range(n - len(e_c))]
            yield e_c
    list3 = list(func_chunk(text3,4))
    list4 = []
    list4n = []

    kl = 0
    list5 = []
    strn =''
    for e in range(len(text3)):
        list5.append(text3[e][1])
        strn+=str(text3[e][0])+"_"
        kl+=1
        if kl%4 == 0:
            arr = np.array(list5)
            print("arr "+str(arr))
            print(arr.shape)
            print(arr.dtype)
            np.save("/mnt/work/nooroka/one-hotnew2npy/{}/{}/{}_{}.npy".format(name2[5],name2[6],name2[-1],strn),arr)
            list5  = []
            strn = ""
         
'''
d = 0
w= open("test6.txt","w")
for h in range(len(text3)):
        d+=1
        w.write(">"+str(text3[h][0])+"\n")
        for xs in text3[h][1]:
            w.write("".join(str(xs))[1:-1]+" "+"\n")
        w.write("\n")
        if d%4 == 0:
            w.write("\n")
            w.write("\n")
            w.write("\n")
            w.write("\n")
w.close()
'''
