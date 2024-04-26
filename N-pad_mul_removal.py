##/Library/Frameworks/Python.framework/Versions/3.6/bin/python3.6
##remove N-pads longer than given lenghth from the input fasta file, generate collapsed contigs.

import sys,getopt
opts,args = getopt.getopt(sys.argv[1:],"hi:l:o:")

input_file = ""
length_threshold = ""
output_file = ""

for op, value in opts:
    if op == "-i":
        input_file = value
    elif op == "-l":
        length_threshold = value
    elif op == "-o":
        output_file = value

import re
from Bio import SeqIO
import progressbar

output = open(output_file,"w")

contig_number = len(list(SeqIO.parse(input_file,"fasta")))

bar = progressbar.ProgressBar(maxval = contig_number).start()
count = 0
for record in SeqIO.parse(input_file,"fasta"):
     a = re.split("N{1000,100000}|(?<!N)N{100}(?!N)",str(record.seq))
     for i,r in enumerate(a):
          output.write(">%s_%s\n" % (record.id,i))
          output.write(r)
          output.write("\n")
          output.flush()
     count += 1
     bar.update(count)

bar.finish()
output.close()

