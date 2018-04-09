#!/usr/bin/python

""" Fastadiff - teensy tiny script for checking NEWICK tree

    Author - Matthew Jobin, Department of Anthropology, Santa Clara University and Department of Anthropology,
        University of California, Santa Cruz.
    """

import argparse
from argparse import RawTextHelpFormatter
from Bio import Phylo

if __name__ == "__main__":


    parser = argparse.ArgumentParser(description="IMPUTsahsahashdn of: \n\n\t"\
                                     "mutations.\n\t"\
                                     "- ",formatter_class=RawTextHelpFormatter)
        
    parser.add_argument('-file',metavar='<file>',help='input file: newick tree', required=True)
    
    args = parser.parse_args()
    inputfile = args.file
    
    
    tree = Phylo.read(inputfile, "newick")
    Phylo.draw_ascii(tree)