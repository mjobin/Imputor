#Imputor

""" Imputor
    
    
    """


import os
import sys
import random
import bz2
import glob
import time
import argparse
from argparse import RawTextHelpFormatter
import pickle
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.Applications import RaxmlCommandline


class InData(object):
    """Input data from Fasta or VCF file converted internally to sequence data.
    """
    
    def __init__(self):
        self.ref_seq = None #Reference sequence used to construct sequence data from VCF
        self.sequence = [] #Sequence data
        self.variantset = set() #Set of locations and states of all variants from reference sequence
        self.variants = {} #Dictionary of each sample and its variation from the reference sequence
        self.tree = None #Phylogenetic tree to be loaded or constructed from data. Newick format.
        self.filebase = None
        self.treeparents = {}
        
        
        #The eight mandatory columns of a VCF file
        self.vcf_chrom = 0
        self.vcf_pos = 1
        self.vcf_id = 2
        self.vcf_ref = 3
        self.vcf_alt = 4
        self.vcf_qual = 5
        self.vcf_filter = 6
        self.vcf_info = 7
        
        return
        
    def load_ref_seq(self, reffile = None):
        """Load a reference sequence for the reconstruction of sequence from
            VCF files.
            
            Keyword arguments:
            reffile -- reference sequence file
            """

        if reffile[-3:]=='obj':
            refobj = open(reffile,'rb')
            self.ref_seq = pickle.load(refobj)
        elif reffile[-3:]=='txt':
            ref_data = open(reffile,'r')
            self.ref_seq = ref_data.read()
        else:
            print "Error. Not a known reference sequence type."
        return

    def load_input_data(self, inputfile = None):
        """ Load input data from a file.
            
            Keyword arguments:
            inputfilile -- Input data. Acceptable formats: FASTA, VCF.
        """
        
        
        filebase, fileext = os.path.splitext(inputfile)
        self.filebase = filebase
        
        if inputfile[-5:]=='fasta':
            self.sequence = AlignIO.read(inputfile, 'fasta')
            self.variants_from_sequence()
        elif inputfile[-3:]=='vcf':

            file_data = open(inputfile, 'r')
            raw_data = []
            for file_line in file_data:
                raw_data.append(file_line.rstrip())
            
            #Split multi-allelic sites
            expanded_file_data = self.vcf_expand_multi_allele(raw_data)

            #Prune non-SNPs
            snps_data = self.vcf_snp_prune(expanded_file_data)
            
            #Generate variants
            for snp_line in snps_data:
                cols = snp_line.split('\t')
                self.variantset.add(cols[self.vcf_pos]+cols[self.vcf_alt])
        
            #Generate sequence from only those areas with any polymorphism
            self.seq_from_variants(raw_data)


        return

    def vcf_expand_multi_allele(self, in_data = None):
        """ Processes multiple ALT alleles in a VCF file. Returns expanded data in list form.
            
            Keyword arguments:
            in_data-- VCF file data.
        """
        expanded_file_data = []
        for file_line in in_data:
            cols = file_line.split('\t')
            if (file_line[:1] == "#") or (cols[self.vcf_chrom]=='\n'): #If the second character is a (meta-info line) or a blank line, ignore
                continue
            if "," in cols[self.vcf_alt]:
                multi_allele  = cols[self.vcf_alt].split(",")
                for allele in multi_allele:
                    allele_cols = cols #A copy of the original line
                    allele_cols[self.vcf_alt] = allele #Replace the ALT column with just this allele
                    expanded_file_data.append("\t".join(allele_cols)) #place line line with just this allele in the expanded data
            else:
                expanded_file_data.append(file_line)
        return expanded_file_data

    def vcf_snp_prune(self, in_data = None):
        """ Returns a VCF file including only lines with SNPs.
            
            Keyword arguments:
            in_data-- VCF file data.
        """
        snps_data = []
        for file_line in in_data:
            cols = file_line.split('\t')
            if (file_line[:1] == "#") or (cols[self.vcf_chrom]=='\n'): #If the second character is a (meta-info line) or a blank line, ignore
                continue
            if len(cols[self.vcf_ref]) > 1 or len(cols[self.vcf_alt]) > 1:  # if not a snp
                continue
            elif cols[self.vcf_ref]=='-' or cols[self.vcf_alt]=='-': # if not a snp
                continue
            elif cols[self.vcf_ref]=='.' or cols[self.vcf_alt]=='.': # if not a snp
                continue
            snps_data.append(file_line)
        return snps_data
    

    def variants_from_sequence(self):
        """ Returns a list of variants from a reference sequence.
            

        """
        for seq_line in self.sequence:
            print seq_line
            if len(seq_line) > len(self.ref_seq):
                print "Error! A sequence line is longer than the reference sequence!"
            diffs = [i for i in xrange(len(self.ref_seq)) if self.ref_seq[i] != seq_line[i]] #All indexes where the sequence does not match
            print diffs
            curdiffs = []
            for diff_pos in diffs:
                self.variantset.add(str(diff_pos)+seq_line[diff_pos]) #Each addition to variantset will thus be unique, since using a set
                curdiffs.append(str(diff_pos) + seq_line[diff_pos])
            self.variants[seq_line.name] = curdiffs
        return


    def seq_from_variants(self, raw_data = None):
        """ Generates sequence data from the variants drived from a VCF file.
            Sequence generated includes only polymorphic sites for the sake of brevity.
            
            Keyword arguments:
            raw_data-- VCF file data.
        """
        for file_line in raw_data:
            cols = file_line.split('\t')
                
            #Locate header line and read genotype names
            if (cols[self.vcf_chrom]=='#CHROM'): #Header line of VCF file
                if(cols[self.vcf_info+1]=='FORMAT'): #On the header line, is there a FORMAT column next to the fixed columns?
                    genotype_names = cols[self.vcf_info+2:] #If so, remaining columns are the genotypes
                    genotype_sequence = {} #Dictionary, thus unique keys. Diploid/triploid sequences will have unique names [name]-#
                    ref_seq_list =[] #The reference sequence is converted to a list so individual sites can be altered
                    for ref_seq_char in self.ref_seq:
                        ref_seq_list.append(ref_seq_char)
                else:
                    print "Error. VCF file with no genotype. Cannot create sequence data."
                    return
        #Step through data lines, constructing list of variants

        #Step through data lines, reconstrucing the sequence of each individual
        snps_data = self.vcf_snp_prune(raw_data) #Ensure only SNPs are being processed
        for file_line in snps_data:
            cols = file_line.split('\t')
            if (file_line[:1] == "#") or (cols[self.vcf_chrom]=='\n' or cols[self.vcf_info+1][:2] != "GT"): #If the second character is a (meta-info line) or a blank line or if GT does not start the format line, ignore
                continue
            indiv_genotypes = cols[self.vcf_info+2:] #Assumes all rows same length, as per VCF standard
            for position, indiv_genotype in enumerate(indiv_genotypes): #Iterates through that row of genotypes for this site
                assigned_alleles = indiv_genotype.split("[/|]+") #Split genotype entry on either character phased or unphased
                changed_genotype_names = []
                for allele_pos, assigned_allele in enumerate(assigned_alleles): #Iterates through the alleles
                    changed_genotype_name = genotype_names[position]
                    if len(assigned_alleles)>1: #Only append to genotype name if not haploid
                        changed_genotype_name = changed_genotype_name + str(allele_pos)
                    changed_genotype_names.append(changed_genotype_name)
                for changed_genotype_name in changed_genotype_names:
                    if changed_genotype_name in genotype_sequence:
                        pass
                    else:
                        genotype_sequence[changed_genotype_name] = [] #Each genotype begins with the reference sequence placed in all valueus
                            
            
                alt_alleles  = cols[self.vcf_alt].split(",") #List of ALT alleles for this row
                for allele_pos, assigned_allele in enumerate(assigned_alleles): #Iterates through the alleles
                    if int(assigned_allele) == 0: #Assigned_allele will be 0 for REF and >0 for any ALT
                        genotype_sequence[changed_genotype_names[allele_pos]].append(cols[self.vcf_ref])
                    else:
                        genotype_sequence[changed_genotype_names[allele_pos]].append(alt_alleles[int(assigned_allele)-1])
                        if changed_genotype_names[allele_pos] in self.variants: #Keys added to self.variants here
                            self.variants[changed_genotype_names[allele_pos]].append( #to avoid empty entries
                                cols[self.vcf_pos]+alt_alleles[int(assigned_allele) - 1])
                        else:
                            self.variants[changed_genotype_names[allele_pos]] = []
                            self.variants[changed_genotype_names[allele_pos]].append(
                                cols[self.vcf_pos]+alt_alleles[int(assigned_allele) - 1])


                        


        for geno in genotype_sequence.keys():
            genotype_sequence[geno] = ''.join(genotype_sequence[geno])
        
        
        #Write to a FASTA file so that it can be read in as a
        outfile = open('vcf_seq_temp.fasta', 'w')
        lenchk = -1
        for geno in genotype_sequence.keys():
            outfile.write(">")
            outfile.write(geno)
            outfile.write("\n")
            outfile.write(genotype_sequence[geno])
            if (lenchk >= 0 and lenchk != len(genotype_sequence[geno])):
                print "mismatch" + lenchk + " " + len(genotype_sequence[geno])
            outfile.write("\n")
        outfile.close()
        self.sequence = AlignIO.read('vcf_seq_temp.fasta', 'fasta')


    def all_parents(self, tree):
        parents = {}
        for clade in tree.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
        return parents

    def input_tree(self, treetype = None, alpha = None, bootstrap = None, rmodel = None):
        """ Generates sequence data from the variants drived from a VCF file.
            Sequence generated includes only polymorphic sites for the sake of brevity.
            
            Keyword arguments:
            treetype -- type of tree to be input or constructed
        """
        
        self.treetype = treetype
        self.tree = None
        
        if self.treetype[-3:]=='xml':
            self.tree = Phylo.read(treetype, "phyloxml")
        elif self.treetype=='parsimony':
            self.parsimony_tree()
        elif self.treetype=='RAxML':
            self.raxml_tree(rmodel)
        else:
            self.phyml_tree(alpha, bootstrap)

        self.treeparents = self.all_parents(self.tree)





    def parsimony_tree(self):
        """ Constructs a tree via maximum parsimony using Biopython's ParsimonyTreeConstructor.
            
        """
        scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher)
        self.tree = constructor.build_tree(self.sequence)



    def raxml_tree(self, rmodel = None):
        """ Constructs a tree via maximum likelihood by invoking external software RAxML.
            See docs for RAxML installation and setup.
        
        """
        
        #Erase RaXML files??? FIXME
        print "Erasing old RAxML output to allow RAxML to run."
        raxml_glob = glob.glob('RAxML_*')
        for delfile in raxml_glob:
            print delfile
            os.remove(delfile)

        
        
        #Output sequence to a temp FASTA file
        tempfastafile = self.filebase + "_fastatmp.fasta"
        AlignIO.write(self.sequence, tempfastafile, "fasta")
        rng = random.SystemRandom() #Uses /dev/urandom
        raxml_cline = RaxmlCommandline(sequences=tempfastafile, model=rmodel, name="imputor",  parsimony_seed = rng.randint(0, sys.maxint))
        out_log, err_log = raxml_cline()
        print err_log
        print out_log
    
        #Import best tree (newick format)
        self.tree  = Phylo.read("RAxML_bestTree.imputor", "newick")
        
        
        #Erase RaXML intermediate files
        for delfile in raxml_glob:
            print delfile
            os.remove(delfile)
        


    def phyml_tree(self, alpha = None, boostrap = None):
        """ Constructs a tree via maximum likelihood by invoking external software PhyML.
            See docs for PhyML installation and setup.
            
        """
        #Output sequence to a temp FASTA file
        tempfastafile = self.filebase + "_fastatmp.fasta"
        AlignIO.write(self.sequence, tempfastafile, "fasta")
        tempphyfile = self.filebase + "_phytmp.phy"
        AlignIO.convert(tempfastafile, "fasta", tempphyfile, "phylip-relaxed")
        cmdline = PhymlCommandline(input=tempphyfile, alpha='e', bootstrap=100)
        print "Commandline for PhyML: " + str(cmdline)
        out_log, err_log = cmdline()
        print err_log
        print out_log
        phytreefile = tempphyfile + "_phyml_tree.txt"
        self.tree  = Phylo.read(phytreefile, "newick")

    def assign_to_haplogroups(self):
        """ Assigns samples to haplogroups.

        """
        print "assign"



    def impute_missing(self):
        """ Imputes missing mutations.

        """
        #get parents
        parents = self.all_parents(self.tree) #These are stored as a dictionary of clades
        print type(parents)
        print parents.keys()
        #Step through all clades

        # for clade1 in self.tree.find_clades(order='level'):
        #     if clade1.name: #If a named clade and thus one that should have a sequence
        #         print "Clade: " + str(clade1) + " " + str(type(clade1))
        #         print "****"
        #         print "Parent: " + str(parents[clade1])
        #         #find that name in sequences ...  SLOW and stupid for now
        #         for chkseq in self.sequence:
        #             if chkseq.name == str(clade1):
        #                 print chkseq.name + " " + chkseq.seq
        #         for clade2 in self.tree.find_clades():
        #             if (str(clade2) == str(clade1)):
        #                 continue
        #             if clade2.name: #If a named clade and thus one that should have a sequence
        #                 cladistance = self.tree.distance(clade1, clade2)
        #                 print str(clade2) + " " + str(cladistance)
        #                 if parents[clade1] == parents[clade2]: #if the two clades share a parent
        #                     print "same parent"
        #         #now sort by distance
        #     print "\n\n"
        #
        # preterms = []
        # #Find all preterminal clades
        # for clade in self.tree.find_clades(order='level'):
        #     if clade.is_preterminal():
        #         preterms.append(clade)
        # print "\n\npretrms"
        # for clade in preterms:
        #     print "Preterm:" + str(clade)
        #     for child in clade:
        #         print child
        #
        # print "\n\nNON"
        # nonterms = self.tree.get_nonterminals() #Get all internal nodes on tree
        # for nonterm in nonterms:
        #     print nonterm
        #     print "***"
        #     for child in nonterm: #
        #         print str(child) + str(child.is_terminal())
        #     print "\n"

        x = 3 #TEST
        for clade1 in self.tree.find_clades(order='level'):
            print "Starting at:" + str(clade1) + " for depth " + str(x)
            empty = []
            print "****"
            kids = self.collect_kids(clade1, empty, 0, x)
            #Build a set of all the variants possessed by any of these samples
            curvariantset = set()
            for child in kids:
                if str(child) in self.variants:
                    for curvariant in self.variants[str(child)]:
                        curvariantset.add(curvariant)
            for child in kids: #Check if a sample does not have a variant others do
                for chkvariant in curvariantset:
                    if str(child) in self.variants:
                        if chkvariant in self.variants[str(child)]:
                            print ""
                        else:
                            total = 0
                            alts = 0
                            for other in kids:
                                if child != other:
                                    if str(other) in self.variants:
                                        if chkvariant in self.variants[str(other)]:
                                            alts = alts +1
                                    total = total + 1
                            print str(alts) + " have the variant "+ str(chkvariant) + " out of a total of " + str(total) + " in this group"


            print "\n"
            #With this list of neighbors, compare sequence



    def collect_kids(self, clade, kids, depth, maxdepth): #Recursive function for collecting all children to specified depth
        if depth < maxdepth:
            for child in clade:
                if child.name:
                    kids.append(child)
                self.collect_kids(child, kids, depth+1, maxdepth)
        return kids



    def all_parents(self, tree):
        parents = {}
        for clade in self.tree.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
        return parents
















if __name__ == "__main__":
    print "\n\n***IMPUTOR ***\n\n"
    
    parser = argparse.ArgumentParser(description="This script does: \n\n\t"\
                                     "- .\n\t"\
                                     "- ",formatter_class=RawTextHelpFormatter)
        
    parser.add_argument('-file',metavar='<file>',help='input file: .fasta, .vcf or .var', required=True)
    parser.add_argument('-ref',metavar='<ref>',help='reference sequence, .txt or .obj', required=True)
    parser.add_argument('-tree',metavar='<tree>',help='tree type; <treefilename.xml>, pars, RAxML, PhyML', default='pars')
    parser.add_argument('-alpha',metavar='<alpha>',help='Value of gamma shape parameter.', default='e')
    parser.add_argument('-boot',metavar='<boot>',help='Number of bootstrap replicates for PhyML.', default='100')
    parser.add_argument('-rmodel',metavar='<rmodel>',help='Model type for RaXML.', default='GTRCAT')
    

    args = parser.parse_args()
    inputfile = args.file
    reffile = args.ref
    treetype = args.tree
    alpha = args.alpha
    bootstrap = args.boot
    rmodel = args.rmodel
    

    print "working in" + os.getcwd() + " on " + inputfile + " and " + reffile + " using " + treetype
    


    indata = InData()
    indata.load_ref_seq(reffile = reffile)
    indata.load_input_data(inputfile = inputfile)
    print "\n****************\nVARIANTS\n****************\n\n"
    #print "all variants"
    #print indata.variantset
    #print "each sample variant"
    #for i in indata.variants.keys():
    #    print i
     #   print indata.variants[i]
    #print "\n****************\nSEQUENCE\n****************\n\n"
    #print indata.sequence
    indata.input_tree(treetype = treetype, alpha = alpha, bootstrap = bootstrap, rmodel = rmodel)
    #print "\n****************\nTREE\n****************\n\n"
    #print indata.tree
    Phylo.draw_ascii(indata.tree)



    #print "\n****************\nHAPLOGROUPS\n****************\n\n"
    #indata.assign_to_haplogroups()

    print "\n****************\nIMPUTATION\n****************\n\n"
    indata.impute_missing()







else:
    print("IMPUTOR is being imported into another module. Not yet implemented.")