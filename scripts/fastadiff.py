# Fastadiff

""" Fastadiff

    Author - Matthew Jobin
    """



import argparse
from argparse import RawTextHelpFormatter

print "\n\n***FASTADIFF ***\n\n"

parser = argparse.ArgumentParser(description="This script compares \n\n\t" \
                                             "- two FASTA files.\n\t" \
                                             "- ",formatter_class=RawTextHelpFormatter)


parser.add_argument('-afile',metavar='<afile>',help='input file A: .fasta', required=True)
parser.add_argument('-bfile',metavar='<bfile>',help='input file B: .fasta')

args = parser.parse_args()
afile = args.afile
bfile = args.bfile

twofile = True

if afile[-5:] != 'fasta':
    print "Error: FASTA files with a .fasta extension only."
    quit()

if not bfile:
	twofile = False
	bfile = afile


if bfile[-5:] != 'fasta':
    print "Error: FASTA files with a .fasta extension only."
    quit()

	

aseq = {}
afile_data = open(afile, 'r')
for file_line in afile_data:
    if file_line[:1] is '>':
        id = file_line[1:].rstrip()
        data_line = next(afile_data).rstrip()
        aseq[id] = data_line

print "Number of read sequences in file" , " " ,len(aseq)


bseq = {}
bfile_data = open(bfile, 'r')
for file_line in bfile_data:
    if file_line[:1] is '>':
        id = file_line[1:].rstrip()
        data_line = next(bfile_data).rstrip()
        bseq[id] = data_line

print "Number of read sequences in file" , " " ,len(aseq)


totdiffs = 0
totmissdiffs = 0
totbases = 0
seqdiffs = {}

if twofile:
#Comparison
	for k in aseq:
		if k in bseq:
			print "****************** " , k , " **********************"
			xaseq = aseq[k]
			xbseq = bseq[k]
			if len(xaseq) != len(xbseq):
				print "Error: lines should be same size for sequences: ", k , " Skipping."
				continue
			diffs = [i for i in xrange(len(xaseq)) if xaseq[i] != xbseq[i]]
			totdiffs += len(diffs)
			totbases += len(xaseq)
			seqdiffs[k] = (float(len(diffs)) / float(len(xaseq)))
			for diff in diffs:
				print xaseq[diff] , " : " , xbseq[diff]
else:
	for k in aseq:
		for l in bseq:
			if k == l:
				continue
			xaseq = aseq[k]
			xbseq = bseq[l]
			if len(xaseq) != len(xbseq):
				print "Error: lines should be same size for sequences: ", k , " Skipping."
				continue
			diffs = [i for i in xrange(len(xaseq)) if xaseq[i] != xbseq[i]]
			print k, "\t" , l, "\t", (float(len(diffs)) / float(len(xaseq)))
			
		
			


print "ID\tDiff"
for x in sorted(seqdiffs.keys()):
	print x, "\t", seqdiffs[x]

print "\nMean Difference: " , float(totdiffs) / float(totbases)










