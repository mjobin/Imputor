# Fastadiff

""" Fastadiff

    Author - Matthew Jobin
    """



import argparse
from argparse import RawTextHelpFormatter

print "\n\n***FASTADIFF ***\n\n"

parser = argparse.ArgumentParser(description="This script does: \n\n\t" \
                                             "- .\n\t" \
                                             "- ",formatter_class=RawTextHelpFormatter)


parser.add_argument('-afile',metavar='<afile>',help='input file A: .fasta', required=True)
parser.add_argument('-bfile',metavar='<bfile>',help='input file B: .fasta', required=True)

args = parser.parse_args()
afile = args.afile
bfile = args.bfile

if afile[-5:] != 'fasta' or bfile[-5:] != "fasta":
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
seqdiffs = []

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
        totdiffs = totdiffs + len(diffs)
        totbases = totbases + len(xaseq)
        seqdiffs.append(float(len(diffs)) / float(len(xaseq)))
        for diff in diffs:
            print xaseq[diff] , " : " , xbseq[diff]


print "Differences: " , float(totdiffs) / float(totbases)










