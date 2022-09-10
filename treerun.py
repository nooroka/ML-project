import os
import sys
os.chdir("/mnt/work/phylobench/2020-Pfam33")
os.system("ls */ncbi-astral.tre > /data/nooroka/ML/astral.log")
os.chdir("/data/nooroka/ML")
op = open("/data/nooroka/ML/astral.log","r")
for line in op:
    line = line.strip()
    a = "/mnt/work/phylobench/2020-pfam33/{}".format(line)
    line2 = line.split("/")
    b1 ="/mnt/work/phylobench/2020-pfam33/"+str(line2[0])+"/"+"Trees/"+"FastME/"
    b2 ="/mnt/work/phylobench/2020-pfam33/"+str(line2[0])+"/"+"Trees/"+"RAxML/"
    b3 ="/mnt/work/phylobench/2020-pfam33/"+str(line2[0])+"/"+"Trees/"+"TNT/"
    b11 = str(line2[0])+"-"+"FastME"
    b22 = str(line2[0])+"-"+"RAxML"
    b33 = str(line2[0])+"-"+"TNT"
    if os.path.exists("{}".format(b1)) and os.path.exists("{}".format(b2)) and os.path.exists("{}".format(b3)):
        for file1 in os.listdir("{}".format(b1)):
            os.system("python /data/nooroka/ML/treequart3.py {} {}{}".format(a,b1,file1))
        for file2 in os.listdir("{}".format(b2)):
            os.system("python /data/nooroka/ML/treequart3.py {} {}{}".format(a,b2,file2))
        for file3 in os.listdir("{}".format(b3)):
            os.system("python /data/nooroka/ML/treequart3.py {} {}{}".format(a,b3,file3))
op.close()

