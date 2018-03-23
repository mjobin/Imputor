#!/usr/bin/python

""" Fastablasta

    Author - Matthew Jobin, Department of Anthropology,
        University of California, Santa Cruz.
    """


import argparse
from argparse import RawTextHelpFormatter
import random
import os

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

parser = argparse.ArgumentParser(description="This script: \n\n\t" \
											 "- takes a fasta infile.\n\t" \
											 "- creates randoml;y-altered copies.\n\t" \
											 "with missing or altered bases .\n\t" \
											 "- ",formatter_class=RawTextHelpFormatter)

# parser.add_argument('-outfile',metavar='<file>',help='output file: .fasta', required=True)
parser.add_argument('-l', metavar='<length>', help='Sequence length.', default='1000')
parser.add_argument('-c', metavar='<count>', help='Number of sequences.', default='1000')
parser.add_argument('-mu', metavar='<mutrate>', help='Mutation rate.', default='0.0')
parser.add_argument('-miss', metavar='<miss>', help='Missing rate.', default='0.0')
parser.add_argument('-tstv', metavar='<tstv>', help='Transition/travsersion ratio.', default='2.0')
parser.add_argument('-infile', metavar='<infile>', help='FASTA input file.', required=True)
parser.add_argument('-indelsize', metavar='<indelsize>', help='Max indel size', default=0)
parser.add_argument('-indelrate', metavar='<indelrate>', help='Indel rate.', default=0.0)
parser.add_argument('-copyrate', metavar='<copyrate>', help='Copy rate.', default=0.01)
parser.add_argument('-reps', metavar='<reps>', help='Replicates.', default=10)


args = parser.parse_args()
# outfile = args.outfile
length = int(args.l)
count = int(args.c)
mu = float(args.mu)
miss = float(args.miss)
tstv = float(args.tstv)
infile = args.infile
idsize = int(args.indelsize)
idrate = float(args.indelrate)
reps = int(args.reps)


acgt = {'A', 'C', 'G', 'T'}
missing = {'.', '-', 'N'}



batchlist = []

if os.path.isdir(infile):
	for bbfile in os.listdir(infile):
		if bbfile.endswith(".fasta"):
			batchlist.append(infile + "/" + bbfile)
elif os.path.isfile(infile):
	if infile.endswith(".fasta"):
		batchlist.append(infile)
	else:
		print "Input file must be .fasta"
		exit()
else:
	print "Input file must be fasta"
	exit()

for infile in batchlist:
# 	inbasename = os.path.basename(afile)
	infilebase, infileext = os.path.splitext(infile)

	bases = []

	workseqs = []
	ref = []

	firstline = True
	firstseq = None
	file_data = open(infile, 'r')
	for file_line in file_data:
		if file_line[:1] is '>':
			newseq = []
			seqname = file_line[1:].rstrip()
			seq = next(file_data).rstrip()
			newseq.append(seqname)
			newseq.append(seq)
			workseqs.append(newseq)
			if firstline:
				firstline = False
				firstseq = seq
				for i in xrange(len(firstseq)):
					bases.append(set())
			else:	
				if len(seq) != len(firstseq):
					print "Error! All lines must be same length!"
					exit()
		

	#work out what bases are at which sites. only these will be errors
	#

	for seq in workseqs:
		seq_line = seq[1]
		for i in xrange(len(seq_line)):
			bases[i].add(seq_line[i])
	
	for rep in range(reps):
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
				elif random.random() < mu and len(bases[i])>1:
					while(True):
						muted = random.choice(tuple(bases[i]))
						if muted is not seq_line[i]:
							out_line.append(muted)
							break
		# 			out_line.append(mutit(seq_line[i],tstv))
				else:
					out_line.append(seq_line[i])
				tot = tot +1
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
		missbase = str(miss)
		misses = missbase.split(".")
		missout = misses[len(misses)-1]
		mutbase = str(mu)
		muts = mutbase.split(".")
		mutout = muts[len(muts)-1]

	
		outfilename = infilebase + "-miss" + missout + "-mu" + mutout + "-" + str(rep) + ".fasta"
		outfile = open(outfilename, 'w')
		lenchk = -1


		for x in outseqs:
			outfile.write(">")
			outfile.write(x[0])
			outfile.write("\n")
			outfile.write(x[1])
			outfile.write("\n")
		outfile.close()
exit()



			


		

		
	



