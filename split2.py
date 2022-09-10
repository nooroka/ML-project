import os
from collections import Counter, defaultdict
#import sklearn.cross_validation.train_test_split
import numpy as np
list1 = ['Actinobacteria-60-genus','Archaea-60-genus','Archaea-60-species','Ascomycota-60-genus', 'Chordata-60-family', 'Eukaryota-60-altbalanced','Eukaryota-60-balanced', 'Eukaryota-60-class','Firmicutes-60-genus','Fungi-60-balanced','Fungi-60-family','Metazoa-60-altbalanced','Metazoa-60-balanced','OtherBacteria-60','Proteobacteria-60-balanced','Proteobacteria-60-family','Protista-60-species','Streptophyta-60-genus']
dictall = defaultdict(list)
listname = []

for i in range(len(list1)):
      list2 = os.listdir("/mnt/work/phylobench/2020-Pfam33/{}/Alignments".format(str(list1[i])))
      for k in range(len(list2)):
          if 'afa' in str(list2[k]):
              name = list2[k].split("_")[1]
              dictall[name].append(list2[k])
              listname.append(name)
#print(len(set(listname)))
setname = set(sorted(listname))
print(len(setname))
train,test,validate, no = np.split(list(setname),[3750,4544,5357])
#print(train)
#print(test)
#print(validate)
#print(no)
train2 = []
test2 = []
validate2 = []
for key in train:
    for i in range(len(dictall[key])):
        train2.append(dictall[key][i])
for key2 in test:
    for k in range(len(dictall[key2])):
        test2.append(dictall[key2][k])
for key3 in validate:
    for j in range(len(dictall[key3])):
        validate2.append(dictall[key3][j])
print(str(validate2) + "validate2")

