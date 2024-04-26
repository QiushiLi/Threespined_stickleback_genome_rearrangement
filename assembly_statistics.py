# N50 calculation python3

print('Python and Biopython needed for running this script!')

print("Script for calculating N50 of assembly")

fasta = input('Enter fasta name: ')

# N50 calculation

BaseSum,valueNo,ValueSum,N50 = 0,0,0,0

Length = []


from Bio import SeqIO

for record in SeqIO.parse(open(fasta), "fasta"):
    
    BaseSum += len(record.seq)
        
    Length.append(len(record.seq))
                


N50_pos = BaseSum/2

Length.sort()

Length.reverse()

for value in Length:
    
    ValueSum += value
    valueNo += 1
        
    if N50_pos <= ValueSum:
            
        N50 = value
                
        break

print('Sequences No.:'+str(len(Length)))

print('Sequences Min.:'+str(min(Length)))

print('Sequences Max.:'+str(max(Length)))

print('N50 Size:' + str(N50))

print('N50 No.:' + str(valueNo))

print('Base count:' + str(BaseSum))
