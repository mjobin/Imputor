#Fastablasta

""" Fastablasta
    
    Author - Matthew Jobin
    """


import argparse
from argparse import RawTextHelpFormatter
import random

def mutit(a, tstv):
    ts = True
    if random.random() < (1.0/(tstv+1.0)):
        ts = False
    if a == "A":
        if ts == True:
            return "G"
        else:
            if random.random() < 0.5:
                return "C"
            else:
                return "T"
    if a == "G":
        if ts == True:
            return "A"
        else:
            if random.random() < 0.5:
                return "C"
            else:
                return "T"
    if a == "C":
        if ts == True:
            return "T"
        else:
            if random.random() < 0.5:
                return "A"
            else:
                return "G"
    if a == "T":
        if ts == True:
            return "C"
        else:
            if random.random() < 0.5:
                return "A"
            else:
                return "G"
    return "-"

print "\n\n***FASTABLASTA ***\n\n"

parser = argparse.ArgumentParser(description="This script does: \n\n\t" \
                                             "- .\n\t" \
                                             "- ",formatter_class=RawTextHelpFormatter)

parser.add_argument('-file',metavar='<file>',help='output file: .fasta', required=True)
parser.add_argument('-l', metavar='<length>', help='Sequence length.', default='1000')
parser.add_argument('-c', metavar='<count>', help='Number of sequences.', default='1000')
parser.add_argument('-mu', metavar='<mutrate>', help='Mutation rate.', default='8.71e-10')
parser.add_argument('-miss', metavar='<miss>', help='Missing rate.', default='0.0')
parser.add_argument('-tstv', metavar='<tstv>', help='Transition/travsersion ratio.', default='2.0')
parser.add_argument('-infile', metavar='<infile>', help='FASTA input file.')



args = parser.parse_args()
outfile = args.file
length = int(args.l)
count = int(args.c)
mu = float(args.mu)
miss = float(args.miss)
tstv = float(args.tstv)
infile = args.infile


workseqs = []

if infile is None:
    ref = []
    for x in range(0, length):
        r = random.randint(0,3)
        if r == 0:
            ref.append("A")
        elif r == 1:
            ref.append("C")
        elif r == 2:
            ref.append("G")
        elif r == 3:
            ref.append("T")

    theseq = ("".join(ref))
    for y in range(0, count):
        newseq = []
        newseq.append(str(y))
        newseq.append(theseq)
        workseqs.append(newseq)


else:
    file_data = open(infile, 'r')
    for file_line in file_data:
        if file_line[:1] is '>':
            newseq = []
            newseq.append(file_line[1:].rstrip())
            newseq.append(next(file_data).rstrip())
            workseqs.append(newseq)

# print workseqs
outseqs = []
#Mutation, missing data
for seq in workseqs:
    outseq = []
    seq_line = seq[1]
    out_line = []
    # tot = 0
    # missed = 0
    for base in seq_line:
        if random.random() < miss:
            out_line.append("-")
            # missed = missed + 1
        elif random.random() < mu:
            out_line.append(mutit(base,tstv))
        else:
            out_line.append(base)
        # tot = tot +1
    # print "miss: ", missed
    # print "tot: ", tot
    outseq.append(seq[0])
    outseq.append("".join(out_line))
    outseqs.append(outseq)




#Write to a FASTA file so that it can be read in as a
outfilename = outfile + ".fasta"
outfile = open(outfilename, 'w')
lenchk = -1


for x in outseqs:
    outfile.write(">")
    outfile.write(x[0])
    outfile.write("\n")
    outfile.write(x[1])
    outfile.write("\n")
outfile.close()



			


		

		
	



