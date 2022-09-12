#почему-то MA_PF00587_1.fasta не добавился
#не котируются файлы с один индексами подряд
from Bio import AlignIO
import random
from collections import Counter
from Bio import SeqIO
import os
filelist = []
for root, dirs, files in os.walk('/mnt/work/nooroka/alignments'):
    for file1 in files:
        filelist.append(os.path.join(root,file1))
#print(filelist)
#for name in filelist:
name = "/data/nooroka/ML/FB_PF18787_1.fasta"
print(name)
handle = open(name)
align = AlignIO.read(handle, "fasta")   
name4 = name.split("/")
print("name4 "+str(name4))
print(name4[-1])
'''
    os.chdir("/mnt/work/nooroka/alignments300")
    if (not os.path.exists("/mnt/work/nooroka/alignments300/{}/{}".format(name4[5],name4[6]))) and (not os.path.exists("/mnt/work/nooroka/alignments300/{}".format(name4[5]))):
          os.mkdir("/mnt/work/nooroka/alignments300/{}".format(name4[5]))
          os.mkdir("/mnt/work/nooroka/alignments300/{}/{}".format(name4[5],name4[6]))
    elif (not os.path.exists("/mnt/work/nooroka/alignments300/{}/{}".format(name4[5],name4[6]))) and  (os.path.exists("/mnt/work/nooroka/alignments/{}".format(name4[5]))):
          os.mkdir("/mnt/work/nooroka/alignments300/{}/{}".format(name4[5],name4[6]))
    os.chdir("/mnt/work/nooroka/alignments300/{}/{}".format(name4[5],name4[6]))
'''
    #колонка из одних гэпов не учитывается 
    #w = open("/mnt/work/nooroka/alignments300/{}/{}/{}.fasta".format(name4[5],name4[6],name4[-1]),"w")
w = open("/data/nooroka/ML/test4.txt","w")
d = 0
for record in align:
        #print(record)
        rec = len(record)
        d+=1
print(align)
listalign = []
w.write(str(name)+"\n")
for m in range(0,d,4): 
        if rec < 300:
            print("OK")
            k = 0
            al = align[m:m+4,0:0]
            #k+=len(align)
            while k <= 300:
                b = align[m:m+4,0:rec]
                al += b
                k+=rec
                print("len rec "+str(rec))
            al2 = al[:,0:300]
            AlignIO.write(al2,w,"fasta")
            w.write('\n')

        if rec > 300:
            a = align[m:m+4,0:300]
            AlignIO.write(a,w,"fasta")
            w.write('\n')

            '''
        num2 = random.randint(1, rec)
        k = 0
        while k<300:
            num = random.randint(1, rec)#сохранять порядок колонок!?
            col = align[m:m+4,num-1:num]
            print(col)
            listcol = [] #список с буквами из колонки
            for i in col:
                listcol.append(i[0])
                
            if len(set(listcol)) == 1 and listcol[0] == "-":
                continue
            if k == 0:
                a = col
                print("a "+str(a))
                k+=1
            else:
                a = col+a
                print("a "+str(a))
                k+=1
            
       # w.write(str(k)+"\n")
        if int(k) == 300:
            listalign.append(a)
            AlignIO.write(a, w, "fasta")
                #print(a)
        
        w.write('\n')
        '''
handle.close()
w.close()
    #os.chdir("/mnt/work/nooroka/alignments300")
