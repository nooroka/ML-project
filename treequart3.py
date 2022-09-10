from ete3 import Tree
from collections import defaultdict,Counter
import numpy as np
import sys
import itertools
import random
import os
from Bio import Phylo
tree2 = Tree('/mnt/work/phylobench/2020-Pfam33/Chordata-60-family/ncbi-astral.tre')
#tree2 = Tree("ncbi-astral.tre") 
#tree2 = Tree(sys.argv[1])
#print(tree2)
list22 = []
for node in tree2.traverse("levelorder"):
    if not node.is_leaf():
        for child2 in node.get_children():
            leaf_names2 = child2.get_leaf_names()
            if len(leaf_names2) > 1:
                list22.append(leaf_names2)
list12 = []
for node in tree2.traverse():
    if node.is_leaf():
        list12.append(node.name)
list312 = ["a"]
list322 = ["a"]
branchesex = defaultdict(list)
count = 0
for i in range(len(list22)):
    if len(list322) > 0 and len(list312) > 0:
        if sorted(list22[i])!=sorted(list322[-1]) and sorted(list22[i])!=sorted(list312[-1]):
            count+=1
            list312.append(sorted(list22[i]))
            list322.append(sorted(list(set(list12)-set(list22[i]))))
            branchesex[count].extend([sorted(list22[i]),sorted(list(set(list12)-set(list22[i])))])
#tree = Tree(sys.argv[2])
#ex = str(sys.argv[2])

#ex ="MA_PF01490_1.tre"
#tree = Tree("MA_PF01490_1.tre") 
#print(branchesex)

ex ="/mnt/work/phylobench/2020-Pfam33/Chordata-60-family/Trees/TNT/CH_PF00900_1.tre"
tree = Tree("/mnt/work/phylobench/2020-Pfam33/Chordata-60-family/Trees/TNT/CH_PF00900_1.tre")
print(ex)
list2 = []
for node in tree.traverse("levelorder"):
    if not node.is_leaf():
        for child in node.get_children():
            leaf_names = child.get_leaf_names()
            if len(leaf_names) > 1: #добавляем нетривиальные ветви 
               # print("leaf_names "+str(leaf_names))
                list2.append(leaf_names)
list1 = []
for node in tree.traverse():
    if node.is_leaf():
        list1.append(node.name)
sysa = ex.split("/")
#sysa = sys.argv[2].split("/")
list31 = ["a"]
list32 = ["a"]
'''
os.chdir("/mnt/work/nooroka/ids150822/")
if not os.path.exists("/mnt/work/nooroka/ids150822/{}".format(sysa[5])):
    os.mkdir("{}".format(sysa[5]))
os.chdir("{}".format(sysa[5]))
if not os.path.exists("/mnt/work/nooroka/ids150822/{}/{}".format(sysa[5],sysa[7])):
    os.mkdir("{}".format(sysa[7]))
w = open("/mnt/work/nooroka/ids150822/{}/{}/{}.txt".format(sysa[5],sysa[7],sysa[-1]),"w")
#w.write(str(list1)+"\n")
'''
w = open('/data/nooroka/ML/test2.txt','w')
branches = []
branchesdouble = []
for i in range(len(list2)):
    if len(list32) > 0 and len(list31) > 0:
        if sorted(list2[i])!=sorted(list32[-1]) and sorted(list2[i])!=sorted(list31[-1]):
           if len(list2[i]) > 1 and len(list(set(list1)-set(list2[i]))) > 1:
                if [sorted(list2[i]),sorted(list(set(list1)-set(list2[i])))] in branchesex.values() or [sorted(list(set(list1)-set(list2[i]))),sorted(list2[i])] in branchesex.values():
                    print("set1 " +str(sorted(list2[i]))+"\t"+"set2 "+str(sorted(list(set(list1)-set(list2[i]))))+ " + ")
                else:
                    print("set1 " +str(sorted(list2[i]))+"\t"+"set2 "+str(sorted(list(set(list1)-set(list2[i]))))+ " - ")
                list31.append(sorted(list2[i]))
                list32.append(sorted(list(set(list1)-set(list2[i]))))
                branches.append([sorted(list2[i]),sorted(list(set(list1)-set(list2[i])))])
                branchesdouble.append([sorted(list(set(list1)-set(list2[i]))),sorted(list2[i])])

dict2= {}
#нужно взять один лист рандомно из каждого листа child2 
#что делать с leafost5?

dictdfs = {}
listall = []
for i in range(len(branches)):
    br = branches[i]
    #w.write(str(br)+"\n")
    root = tree.get_tree_root()
    for node in tree.traverse("preorder"):
        if not node.is_leaf():
            leaf_names5 = node.get_leaf_names()
            leafost5 = sorted(list(set(list1)-set(leaf_names5))) #вторая часть разбиения
            #dictdfs[tuple(leaf_names5)] = "black"
            #dictdfs[tuple(leafost5)] = "black"
            listall.append(leafost5)
            if ((sorted(leaf_names5) == br[0] and sorted(leafost5) == br[1])) or  (sorted(leaf_names5) == br[1] and sorted(leafost5) == br[0]):#без условия смотрит и корень
              if tuple(sorted(leaf_names5)) not in dictdfs.keys():
                #print("NONO")
               # print("dictdfs " +str(dictdfs))

                children = node.get_children()
                if tuple(sorted(children[0].get_leaf_names())) not in dictdfs.keys() and tuple(sorted(children[1].get_leaf_names())) not in dictdfs.keys():
                    print("children[0] "+str(sorted(children[0].get_leaf_names())))#другая часть разбиения мб представлена в виде двух нод и тд
                    print("children[1] "+str(sorted(children[1].get_leaf_names())))
                    print("node leaves "+str(sorted(leaf_names5)))
                    print("leafost5 "+str(sorted(leafost5)))
                    dictdfs[tuple(sorted(children[0].get_leaf_names()))] = "black"
                    dictdfs[tuple(sorted(children[1].get_leaf_names()))] = "black"
                    sp1 = random.choice(children[0].get_leaf_names())
                    sp2 = random.choice(children[1].get_leaf_names())
                    print("sp1 "+str(sp1))
                    print("sp2 "+str(sp2))
                    w.write(str(sp1)+"\n")
                    w.write(str(sp2)+"\n")
                sister = node.get_sisters()
                sisterbr = sister[0].get_leaf_names()
                sisterost = list(set(sisterbr) ^ set(leafost5))
                #if tuple(sorted(sisterbr)) not in dictdfs.keys() and tuple(sorted(sisterost)) not in dictdfs.keys():
                print("sisterleaf "+str(sorted(sister[0].get_leaf_names())))
                print("sister "+str(sister[0]))
                #dictdfs[tuple(sorted(sister[0].get_leaf_names()))] = "black" 
                sp3 = random.choice(sisterbr)
                sp4 = random.choice(sisterost)
                #dictdfs[tuple(sisterost)] = "black"
                print("sp3 "+str(sp3))
                print("sp4 "+str(sp4))
                w.write(str(sp3)+"\n")
                w.write(str(sp4)+"\n")
                print("\n\n")
                w.write("\n")

                  

                 
                   
                    
                    
    
print("dictdfs "+str(dictdfs))
'''
for key in dictdfs:
    print(str(key)+" "+str(dictdfs[key]))
'''                      
w.close()
os.chdir("/data/nooroka/ML")
                    
