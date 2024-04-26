import os, re
from operator import itemgetter
from Bio import SeqIO

os.chdir("/Users/Qiushi/Downloads/test_10_debris_mummer/txt")

path = "/Users/Qiushi/Downloads/test_10_debris_mummer/txt"
files= os.listdir(path)

output = open("/Users/Qiushi/Downloads/test_10_debris_mummer/test_10_homology_coordinate.txt", "w")

items = []

for i in files:

    fastafile = "/Users/Qiushi/Downloads/test_10_debris/"+i[8:-16]+"fa"
    for record in SeqIO.parse(fastafile,"fasta"):
        ilength = len(record.seq)
##        print ilength
    f = open(i,"r")
    context = f.readlines()[5:]
    d = {}
    for line in context:
        pattern = re.compile(r'\D*(\d+)\D+(\d+)\D+\d+\D+\d+\D+\d+\D+(\d+)\D+\d+[.]\d+\D+(HiC\D+\d+)')
        match = pattern.match(line)
        
        if match.group(4) in d:
            d[match.group(4)] = [d[match.group(4)][0],int(match.group(2)),d[match.group(4)][2]+int(match.group(3))]
        else:
            d[match.group(4)] = [int(match.group(1)),int(match.group(2)),int(match.group(3))]

##    print d
    if not d:
        pass
    else:
        candidate = sorted(d.items(),key=lambda x: x[1][2])[-1]
##        print candidate
        if candidate[1][2] >= ilength * 0.75:
            items.append(record.id+"\t"+str(candidate[0])+":"+str(candidate[1][0])+"-"+str(candidate[1][1])+"\t"+str(candidate[1][2]))
        else:
            pass
    f.close()

##print items

output.write("\n".join(items))
output.close()

#########################

os.chdir("/Users/Qiushi/Downloads/test_9_debris_mummer/txt")

path = "/Users/Qiushi/Downloads/test_9_debris_mummer/txt"
files= os.listdir(path)

output = open("/Users/Qiushi/Downloads/test_9_debris_mummer/test_9_homology_coordinate.txt", "w")

items = []

for i in files:

    fastafile = "/Users/Qiushi/Downloads/test_9_debris/"+i[7:-16]+"fa"
    for record in SeqIO.parse(fastafile,"fasta"):
        ilength = len(record.seq)
##        print ilength
    f = open(i,"r")
    context = f.readlines()[5:]
    d = {}
    for line in context:
        pattern = re.compile(r'\D*(\d+)\D+(\d+)\D+\d+\D+\d+\D+\d+\D+(\d+)\D+\d+[.]\d+\D+(HiC\D+\d+)')
        match = pattern.match(line)
        
        if match.group(4) in d:
            d[match.group(4)] = [d[match.group(4)][0],int(match.group(2)),d[match.group(4)][2]+int(match.group(3))]
        else:
            d[match.group(4)] = [int(match.group(1)),int(match.group(2)),int(match.group(3))]

##    print d
    if not d:
        pass
    else:
        candidate = sorted(d.items(),key=lambda x: x[1][2])[-1]
##        print candidate
        if candidate[1][2] >= ilength * 0.75:
            items.append(record.id+"\t"+str(candidate[0])+":"+str(candidate[1][0])+"-"+str(candidate[1][1])+"\t"+str(candidate[1][2]))
        else:
            pass
    f.close()

##print items

output.write("\n".join(items))
output.close()

