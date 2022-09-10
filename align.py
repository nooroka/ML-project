import os
from Bio import SeqIO

def f(lst, n):
     return [lst[i:i + n] for i in range(0, len(lst), n)]
filelist = []
for root, dirs, files in os.walk('/mnt/work/nooroka/ids150822/'):
    for file1 in files:
        filelist.append(os.path.join(root,file1))
for name in filelist:
    op = open(name,"r")
    listid = []
    for line in op:
        if len(line) > 1:
               line = line.strip()
               listid.append(line)
    op.close()
    listend = f(listid, 4)
    print(str(name)+" "+str(listend))
    name2 = name.split("/")
    print(name2)
    os.chdir("/mnt/work/nooroka/alignments")  
    if (not os.path.exists("/mnt/work/nooroka/alignments/{}/{}".format(name2[5],name2[6]))) and (not os.path.exists("/mnt/work/nooroka/alignments/{}".format(name2[5]))):
          os.mkdir("/mnt/work/nooroka/alignments/{}".format(name2[5]))
          os.mkdir("/mnt/work/nooroka/alignments/{}/{}".format(name2[5],name2[6]))
    elif (not os.path.exists("/mnt/work/nooroka/alignments/{}/{}".format(name2[5],name2[6]))) and  (os.path.exists("/mnt/work/nooroka/alignments/{}".format(name2[5]))):
          os.mkdir("/mnt/work/nooroka/alignments/{}/{}".format(name2[5],name2[6]))
    os.chdir("/mnt/work/nooroka/alignments/{}/{}".format(name2[5],name2[6]))
    name3 = name2[7].split(".")
    w = open("{}.fasta".format(name3[0]),"w")
    w.write(str(name)+"\n")
    handle = open("/mnt/work/phylobench/2020-Pfam33/{}/Alignments/{}.afa".format(name2[5],name3[0]),"r")
    dict1 = {}
    for record in SeqIO.parse(handle, "fasta") :
        dict1[record.id] = record.seq
    print(dict1)
    for i in range(len(listend)):
        for d in range(len(listend[i])):
            if listend[i][d] in dict1.keys():
                w.write(">"+str(listend[i][d])+"\n")
                w.write(str(dict1[listend[i][d]])+"\n")
        w.write("\n")
    w.close()
    os.chdir("/mnt/work/nooroka/ids150822")

