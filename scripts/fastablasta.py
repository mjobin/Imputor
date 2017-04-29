#Fastablasta

""" Fastablasta

	Author - Matthew Jobin
	"""


import argparse
from argparse import RawTextHelpFormatter
import random

def mutit(x, tstv):
	ts = True
	if random.random() < (1.0/(tstv+1.0)):
		ts = False
	if x == "A":
		if ts:
			return "G"
		else:
			if random.random() < 0.5:
				return "C"
			else:
				return "T"
	if x == "G":
		if ts:
			return "A"
		else:
			if random.random() < 0.5:
				return "C"
			else:
				return "T"
	if x == "C":
		if ts:
			return "T"
		else:
			if random.random() < 0.5:
				return "A"
			else:
				return "G"
	if x == "T":
		if ts:
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
parser.add_argument('-reffile',metavar='<reffile>',help='reference sequence file')
parser.add_argument('-l', metavar='<length>', help='Sequence length.', default='1000')
parser.add_argument('-c', metavar='<count>', help='Number of sequences.', default='1000')
parser.add_argument('-mu', metavar='<mutrate>', help='Mutation rate.', default='8.71e-10')
parser.add_argument('-miss', metavar='<miss>', help='Missing rate.', default='0.0')
parser.add_argument('-tstv', metavar='<tstv>', help='Transition/travsersion ratio.', default='2.0')
parser.add_argument('-infile', metavar='<infile>', help='FASTA input file.')
parser.add_argument('-indelsize', metavar='<indelsize>', help='Max indel size', default=0)
parser.add_argument('-indelrate', metavar='<indelrate>', help='Indel rate.', default=0.0)
parser.add_argument('-skew', metavar='<skew>', help='Skew towards REF calls.', default=0.0)



args = parser.parse_args()
outfile = args.file
reffile = args.reffile
length = int(args.l)
count = int(args.c)
mu = float(args.mu)
miss = float(args.miss)
tstv = float(args.tstv)
infile = args.infile
idsize = int(args.indelsize)
idrate = float(args.indelrate)
skew = float(args.skew)

acgt = {'A', 'C', 'G', 'T'}
missing = {'.', '-', 'N'}

workseqs = []
ref = []
if infile is None:

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
			seqname = file_line[1:].rstrip()
			seq = next(file_data).rstrip()
			newseq.append(seqname)
			newseq.append(seq)
			workseqs.append(newseq)

# print workseqs
outseqs = []
#Mutation, missing data
for seq in workseqs:
	outseq = []
	seq_line = seq[1]
	out_line = []
	tot = 0
	missed = 0
	for i in xrange(len(seq_line)):
#     for base in seq_line:
		if random.random() < miss:
			out_line.append("-")
			missed = missed + 1
		elif random.random() < mu:
			x = mutit(seq_line[i],tstv)
			print x, " ", str(ref[i])
			if seq_line[i] == ref[i]:
				if random.random() < skew:
					out_line.append(x)
				else:
					out_line.append(seq_line[i])
			else:
				out_line.append(x)
		else:
			out_line.append(seq_line[i])
		tot = tot +1
	print "miss: ", missed
	print "tot: ", tot
	for i in xrange(len(out_line)):
		if random.random() < idrate:
			indel = random.randint(-idsize, idsize)
			if indel < 0:
				for x in range(0, indel):
					del out_line[i]
			elif indel > 0:
				for x in range(0, indel):
					insb = "-"
					r = random.randint(0,3)
					if r == 0:
						insb = "A"
					elif r == 1:
						insb = "C"
					elif r == 2:
						insb = "G"
					elif r == 3:
						insb = "T"
					out_line.insert(i, insb)


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



			


		

		
	



