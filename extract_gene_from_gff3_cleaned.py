from Bio import SeqIO

record_dict = SeqIO.index("stickleback_broad.fasta", "fasta")
gene_sequence = open('gene.fasta', 'w')

gene = {}

with open('stb_stickleback_broad.gff3', 'r') as f:
    for line in f:
        line1 = line.strip().split()
        chr = line1[0]
        feature = line1[2]
        start = line1[3]
        end = line1[4]
        direction = line1[6]
        name = line1[8].split(";")[1][5:]
        if feature == 'gene':
            gene[name] = (chr, start, end, direction)

# get gene sequence include introns
for key, value in gene.items():
    if value[3] == '+':
        gene_sequence.write('>%s\n%s\n' % (key, record_dict[value[0]][int(value[1])-1:int(value[2])].seq))
    if value[3] == '-':
        gene_sequence.write('>%s\n%s\n' % (key, record_dict[value[0]][int(value[1]) - 1:int(value[2])].seq.reverse_complement()))


gene_sequence.close()
