#!/usr/bin/python

""" VCFstrip - script for stripping second half od diploid VCF cells. Use with caution!

    Author - Matthew Jobin, Department of Anthropology, Santa Clara University and Department of Anthropology,
        University of California, Santa Cruz.
    """

import argparse
from argparse import RawTextHelpFormatter
import random
import os

print "\n\n***vcfstrip ***\n\n"

parser = argparse.ArgumentParser(description="This script does: \n\n\t" \
                                             "- .\n\t" \
                                             "- ", formatter_class=RawTextHelpFormatter)

parser.add_argument('-infile', metavar='<infile>', help='VCF input file.', required=True)
args = parser.parse_args()
infile = args.infile

# The eight mandatory columns of a VCF file. Here for clarity in functions below.
vcf_chrom = 0
vcf_pos = 1
vcf_id = 2
vcf_ref = 3
vcf_alt = 4
vcf_qual = 5
vcf_filter = 6
vcf_info = 7
acgt = ['A', 'C', 'G', 'T']
missing = ['-', 'N', '?']

file_data = open(infile, 'r')
raw_data = []
for file_line in file_data:
    if len(file_line.rstrip()) > 0:  # Strip blank lines
        raw_data.append(file_line.rstrip())
file_data.close()
outlines = []

for file_line in raw_data:
    cols = file_line.split('\t')
    if (file_line[:1] == "#") or (cols[vcf_chrom] == '\n' or cols[vcf_info + 1][:2] != "GT"):
        outlines.append(file_line)
    else:
        for j in range(vcf_info + 2, len(cols)):
            ggcol = cols[j].split("|")
            cols[j] = ggcol[0]
        jline = "\t".join(cols)
        outlines.append(jline)

filebase, fileext = os.path.splitext(infile)
outfilename = filebase + "-s.vcf"
outfile = open(outfilename, 'w')
for outline in outlines:
    outfile.write(outline)
    outfile.write("\n")
outfile.close()
exit()
