## python 2.7.14
##calculate the overall blast rate for each gene from the blast result.

import os
os.getcwd()
'/Users/Qiushi/Documents'
os.chdir('/Users/Qiushi/Documents/stickleback/lists_20181218')
f = open("Peichel_job_3.hints.pep.reformat.uniq.full_vs_repeatpep.best.blast","r")

###a = {jg10096:{length:341,ranges:[ (1,189),(243,340),(191,238) ] },jg11557:{length:384,ranges:[(28,383)]}}
out_blast_rate = open('/Users/Qiushi/Documents/stickleback/lists_20181218/blast_rate.txt',"w")
a ={}
for line in f.readlines():
     if line.split("\t")[0] not in a.keys():
          a[line.split("\t")[0]] = {}
          a[line.split("\t")[0]]["length"] = line.split("\t")[3][:-1]
          a[line.split("\t")[0]]["ranges"] = []
          a[line.split("\t")[0]]["ranges"].append((line.split("\t")[1],line.split("\t")[2]))
     else:
          a[line.split("\t")[0]]["ranges"].append((line.split("\t")[1],line.split("\t")[2]))

pos = []
for i in a:
     mask = [0]*int(a[i]["length"])
     for j in a[i]["ranges"]:
          for k in range(int(j[0]),int(j[1])):
               mask[k] = 1
##or int/int get 0 or 1
     blast_rate = float(sum(mask))/len(mask)
##display two digits
     pos.append((i,round(blast_rate,2)))

###print pos

for i in pos:
     out_blast_rate.write(i[0] + "\t" + str(i[1]) + "\n")

f.close()
out_blast_rate.close()
     


     
      

