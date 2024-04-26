### when headers in R1.fasta and R2.fasta are identical, keep the paired reads, and rename as "species_rna_tissue_readID /1" and "species_rna_tissue_readID /2"
import sys
R1_file = open(sys.argv[1],'r')
R2_file = open(sys.argv[2],'r')

id = []

tissue = sys.argv[1][12:-21]

R1_dict = {}
for line in R1_file:
    if line.startswith(">"):
        name = line[1:-1]
        R1_dict[name] = ""
        id.append(name)
    else:
        R1_dict[name] += line[:-1]

R2_dict = {}
for line in R2_file:
    if line.startswith(">"):
        name = line[1:-1]
        R2_dict[name] = ""
        id.append(name)
    else:
        R2_dict[name] += line[:-1]

from collections import Counter
id_dict = Counter(id)

def get_paired_id(dict,val):
    L = []
    for key,val in dict.items():
        if val == 2:
            L.append(key)
    return L

paired_id = get_paired_id(id_dict,2)

R1_paired_output = open("stickleback_"+tissue+"_R1.fasta","w")
R2_paired_output = open("stickleback_"+tissue+"_R2.fasta","w")

for i in paired_id:
    R1_paired_output.write(">stickleback_"+tissue+"_"+i+" /1"+"\n"+R1_dict[i]+"\n")
    R2_paired_output.write(">stickleback_"+tissue+"_"+i+" /2"+"\n"+R2_dict[i]+"\n")

R1_paired_output.close()
R2_paired_output.close()

        


    

    
