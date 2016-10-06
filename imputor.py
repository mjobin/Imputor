#Imputor

""" Imputor - software for imputing missing mutations in sequencing data by the use
    of phylogenetic trees.
    
    Author - Matthew Jobin
    """


import os
import sys
import random
import glob
import argparse
from argparse import RawTextHelpFormatter
import pickle
import threading
import multiprocessing
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment


class InData(object):
    """Input data from FASTA or VCF file converted internally to sequence data.
        Prepares a dictionary of variants and a list of sequences for processing
        and imputation.
    """
    
    def __init__(self):
        self.ref_seq = None #Reference sequence used to construct sequence data from VCF
        self.sequence = [] #Sequence data
        self.variantset = set() #Set of locations and states of all variants from reference sequence
        self.variants = {} #Dictionary of each sample and its variation from the reference sequence
        self.filebase = None

        # The eight mandatory columns of a VCF file. Here for clarity in functions below.
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

        if reffile[-3:] == 'obj':
            refobj = open(reffile, 'rb')
            self.ref_seq = pickle.load(refobj)
        elif reffile[-3:] == 'txt':
            ref_data = open(reffile, 'r')
            self.ref_seq = ref_data.read()
        else:
            print "Error. Not a known reference sequence type."
        return

    def load_input_data(self, inputfile = None):
        """ Load input data from a file.
            
            Keyword arguments:
            inputfile -- Input data. Accepted formats: FASTA, VCF.
        """
        
        filebase, fileext = os.path.splitext(inputfile)
        self.filebase = filebase
        
        if inputfile[-5:] == 'fasta':
            self.sequence = AlignIO.read(inputfile, 'fasta')
            self.variants_from_sequence()
        elif inputfile[-3:] == 'vcf':

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

    def vcf_expand_multi_allele(self, in_data=None):
        """ Processes multiple ALT alleles in a VCF file. Returns expanded data in list form.
            
            Keyword arguments:
            in_data-- VCF file data.
        """
        expanded_file_data = []
        for file_line in in_data:
            cols = file_line.split('\t')
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n'): #If the second character is a (meta-info line) or a blank line, ignore
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

    def vcf_snp_prune(self, in_data=None):
        """ Returns a VCF file including only lines containing SNPs.
            
            Keyword arguments:
            in_data-- VCF file data.
        """
        snps_data = []
        for file_line in in_data:
            cols = file_line.split('\t')
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n'): #If the second character is a (meta-info line) or a blank line, ignore
                continue
            if len(cols[self.vcf_ref]) > 1 or len(cols[self.vcf_alt]) > 1:  # if not a snp
                continue
            elif cols[self.vcf_ref] == '-' or cols[self.vcf_alt] == '-': # if not a snp
                continue
            elif cols[self.vcf_ref] == '.' or cols[self.vcf_alt] == '.': # if not a snp
                continue
            snps_data.append(file_line)
        return snps_data

    def variants_from_sequence(self):
        """ Returns a list of variants from a reference sequence.
            Intended for use with FASTA input, but will work with any AlignIO object or
            list of sequence data.
        """
        for seq_line in self.sequence:
            if len(seq_line) > len(self.ref_seq):
                print "Error! A sequence line is longer than the reference sequence!"
            diffs = [i for i in xrange(len(self.ref_seq)) if self.ref_seq[i] != seq_line[i]] # All indexes where the sequence does not match
            curdiffs = []
            for diff_pos in diffs:
                self.variantset.add(str(diff_pos)+seq_line[diff_pos]) # Each addition to variantset will thus be unique, since using a set
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
                
            # Locate header line and read genotype names
            if (cols[self.vcf_chrom]=='#CHROM'): #H eader line of VCF file
                if(cols[self.vcf_info+1]=='FORMAT'): # On the header line, is there a FORMAT column next to the fixed columns?
                    genotype_names = cols[self.vcf_info+2:] # If so, remaining columns are the genotypes
                    genotype_sequence = {} # Dictionary, thus unique keys. Diploid/triploid sequences will have unique names [name]-#
                    ref_seq_list =[] # The reference sequence is converted to a list so individual sites can be altered
                    for ref_seq_char in self.ref_seq:
                        ref_seq_list.append(ref_seq_char)
                else:
                    print "Error. VCF file with no genotype. Cannot create sequence data."
                    return

        # Step through data lines, constructing list of variants. Add an empty list for every genotype.
        for genotype_name in genotype_names:
            self.variants[genotype_name] = []

        # Step through data lines, reconstructing the sequence of each individual
        snps_data = self.vcf_snp_prune(raw_data) #Ensure only SNPs are being processed
        for file_line in snps_data:
            cols = file_line.split('\t')
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n' or cols[self.vcf_info+1][:2] != "GT"): #If the second character is a (meta-info line) or a blank line or if GT does not start the format line, ignore
                continue
            indiv_genotypes = cols[self.vcf_info+2:] # Assumes all rows same length, as per VCF standard
            for position, indiv_genotype in enumerate(indiv_genotypes): # Iterates through that row of genotypes for this site
                assigned_alleles = indiv_genotype.split("[/|]+") # Split genotype entry on either character phased or unphased
                changed_genotype_names = []
                for allele_pos, assigned_allele in enumerate(assigned_alleles): # Iterates through the alleles
                    changed_genotype_name = genotype_names[position]
                    if len(assigned_alleles)>1: # Only append to genotype name if not haploid
                        changed_genotype_name = changed_genotype_name + str(allele_pos)
                    changed_genotype_names.append(changed_genotype_name)
                for changed_genotype_name in changed_genotype_names:
                    if changed_genotype_name in genotype_sequence:
                        pass
                    else:
                        genotype_sequence[changed_genotype_name] = [] # Each genotype begins with the reference sequence placed in all valueus
                            
            
                alt_alleles  = cols[self.vcf_alt].split(",") #List of ALT alleles for this row
                for allele_pos, assigned_allele in enumerate(assigned_alleles): #Iterates through the alleles
                    if assigned_allele == "0": #Assigned_allele will be 0 for REF and >0 for any ALT
                        genotype_sequence[changed_genotype_names[allele_pos]].append(cols[self.vcf_ref])
                    elif assigned_allele == ".": #VCF format code for missing allele
                        genotype_sequence[changed_genotype_names[allele_pos]].append("N")
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


class PhyloTree(object):
    """A phylogenetic tree either input from phyloxml format or constructed from sequence
    """

    def __init__(self):
        self.tree = None  # Phylogenetic tree to be loaded or constructed from data. Newick format.
        self.treeparents = {}

    def input_tree(self, treetype=None, alpha=None, bootstrap=None, rmodel=None, indata=None):
        """ Takes input tree file or sequence data.

            Keyword arguments:
            treetype -- type of tree to be input or constructed
        """

        self.treetype = treetype
        self.tree = None

        if self.treetype[-3:] == 'xml':
            self.tree = Phylo.read(treetype, "phyloxml")
        elif self.treetype == 'parsimony':
            self.parsimony_tree()
        elif self.treetype == 'RAxML':
            self.raxml_tree(rmodel)
        else:
            self.phyml_tree(alpha, bootstrap)

        self.treeparents = self.all_parents(self.tree)

    def all_parents(self, tree):
        parents = {}
        for clade in tree.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
        return parents

    def collect_kids(self, clade, kids, depth, maxdepth):
        """ Recursive function for collecting all children to specified depth

            Keyword arguments:
            clade -- current node of tree
            kids --lost of all collected children
            depth -- current depth of search
            maxdepth -- maximum depth of search
        """
        #
        if depth < maxdepth:
            for child in clade:
                if child.name:
                    kids.append(child)
                self.collect_kids(child, kids, depth+1, maxdepth)
        return kids

    def collect_all_kids(self, clade, kids):
        """ Recursive function for collecting all children of clade

            Keyword arguments:
            clade -- current node of tree
            kids -- lost of all collected children
        """
        for child in clade:
            if child.name:
                kids.append(child)
            self.collect_all_kids(child, kids)
        return kids

    def parsimony_tree(self):
        """ Constructs a tree via maximum parsimony using Biopython's ParsimonyTreeConstructor.

        """
        scorer = ParsimonyScorer()
        searcher = NNITreeSearcher(scorer)
        constructor = ParsimonyTreeConstructor(searcher)
        self.tree = constructor.build_tree(indata.sequence)

    def raxml_tree(self, rmodel=None):
        """ Constructs a tree via maximum likelihood by invoking external software RAxML.
            See docs for RAxML installation and setup.

            Keyword arguments:
            rmodel -- model type for input into RAxML.

        """

        # Erase RaXML intermediate files from previous runs
        raxml_glob = glob.glob('RAxML_*')
        for delfile in raxml_glob:
            os.remove(delfile)

        # Output sequence to a temp FASTA file
        tempfastafile = indata.filebase + "_fastatmp.fasta"
        AlignIO.write(indata.sequence, tempfastafile, "fasta")
        rng = random.SystemRandom()  # Uses /dev/urandom
        raxml_cline = RaxmlCommandline(sequences=tempfastafile, model=rmodel, name="imputor",
                                       parsimony_seed=rng.randint(0, sys.maxint))
        out_log, err_log = raxml_cline()
        #print err_log
        #print out_log

        # Import best tree (newick format)
        self.tree = Phylo.read("RAxML_bestTree.imputor", "newick")

        # Erase RaXML intermediate files
        for delfile in raxml_glob:
            os.remove(delfile)

    def phyml_tree(self, alpha=None, boostrap=None):
        """ Constructs a tree via maximum likelihood by invoking external software PhyML.
            See docs for PhyML installation and setup.

        """
        # Output sequence to a temp FASTA file
        tempfastafile = indata.filebase + "_fastatmp.fasta"
        AlignIO.write(indata.sequence, tempfastafile, "fasta")
        tempphyfile = indata.filebase + "_phytmp.phy"
        AlignIO.convert(tempfastafile, "fasta", tempphyfile, "phylip-relaxed")
        cmdline = PhymlCommandline(input=tempphyfile, alpha='e', bootstrap=100)
        print "Commandline for PhyML: " + str(cmdline)
        out_log, err_log = cmdline()
        #print err_log
        #print out_log
        phytreefile = tempphyfile + "_phyml_tree.txt"
        self.tree = Phylo.read(phytreefile, "newick")


class Imputation(object):
    """Imputation of missing mutations given input data and a phylogenetic tree.
    """

    def __init__(self, indata, tree):
        self.cpucount = multiprocessing.cpu_count()
        self.workseq = {}
        self.imputedseq = MultipleSeqAlignment([])

    def impute(self, depth):
        """ Sets up multiprocessing of imputation function for all terminal nodes.

            Keyword arguments:
            depth -- Depth of search up and down tree to find neighbours.
        """
        for seq in indata.sequence:
            self.workseq[seq.name] = list(str(seq.seq))  # Begin with output sequence matching input
        impute_threads = []
        terms = phytree.tree.get_terminals()  # Get all internal nodes on tree. These are the ones with samples.
        random.shuffle(terms)  # Randomize list so no ordering effects
        for term in terms:
            t = threading.Thread(target=self.impute_missing, args=(term, depth,))
            t.start()
            impute_threads.append(t)
        for thread in impute_threads:  # Block until all complete
            thread.join()
        self.process_imputed()

    def impute_missing(self, term, depth):
        """Imputes missing mutations.

            Keyword arguments:
            term -- Terminal node to be compared to neighbours.
            depth -- Depth of search up and down tree to find neighbours.
        """

        termvars = indata.variants[str(term)]  # extract all the variants
        theparent = None
        curnode = term
        for x in range(0, depth):  #Descend tree to user-specified depth
            if curnode in phytree.treeparents:  #Will not go past the root clade
                partemp = theparent
                theparent = phytree.treeparents[curnode]  #Should be an internal clade with no associated sample
                curnode = partemp
        if theparent:  # If we have found a valid parent clade at specififed depth
            empty = []
            allkids = phytree.collect_kids(theparent, empty, 0, depth)
            neighbour_variantset = set()  # Set of all variants for this particular group of neighbours
            for kid in allkids:
                thesevars = indata.variants[str(kid)]  # Extract all the variants
                for thisvar in thesevars:
                    neighbour_variantset.add(thisvar)
            for curvar in neighbour_variantset:
                present = 0
                curpresent = False
                for kid in allkids:
                    thesevars = indata.variants[str(kid)]  # extract all the variants
                    if curvar in thesevars:
                        present = present + 1
                        if term == kid:
                            curpresent = True
                if curpresent:
                    pass
                else:
                   if len(allkids) > 2 and float(present) / float(len(allkids)) > 0.5:
                        self.workseq[str(term)][int(curvar[:-1])] = curvar[-1:]

    def process_imputed(self):
        for key, value in self.workseq.iteritems():
            seqrec = SeqRecord(Seq("".join(value)), id=key, name=key, description="Imputed Sequence")
            self.imputedseq.append(seqrec)


class ChromStats(object):
    """Data from object or text files for assessment of false negative rates.
    """

    def __init__(self, indata, haplogroupfile, weightsfile, polysfile):
        self.haplogroups = {}  # Trusted haplogroups
        self.weights = {}
        self.polys = {}

        if haplogroupfile[-3:] == 'obj': #object file
            haplogroupobj = open(haplogroupfile, 'rb')
            self.haplogroups = pickle.load(haplogroupobj)
        elif haplogroupfile[-3:] == 'txt': #text file
            file_data = open(haplogroupfile, 'r')
            raw_data = []
            for file_line in file_data:
                raw_data.append(file_line.rstrip())

        if weightsfile[-3:] == 'obj': #object file
            weightsobj = open(weightsfile, 'rb')
            self.weights = pickle.load(weightsobj)
        elif weightsfile[-3:] == 'txt': #text file
            file_data = open(weightstfilefile, 'r')
            raw_data = []
            for file_line in file_data:
                raw_data.append(file_line.rstrip())

        if polysfile[-3:] == 'obj': #object file
            polysobj = open(polysfile, 'rb')
            self.polys = pickle.load(polysobj)
        elif polysfile[-3:] == 'txt': #text file
            file_data = open(polysfile, 'r')
            raw_data = []
            for file_line in file_data:
                raw_data.append(file_line.rstrip())


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
    parser.add_argument('-depth', metavar='<depth>', help='Depth of search toward root for collecting neeighbors.', default=2)
    parser.add_argument('-hapobj', metavar='<hapobj>', help='Haplogroups file.')
    parser.add_argument('-wtobj', metavar='<wtobj>', help='Weights file.')
    parser.add_argument('-polyobj', metavar='<polyobj>', help='Polymorphisms file.')


    args = parser.parse_args()
    inputfile = args.file
    reffile = args.ref
    treetype = args.tree
    alpha = args.alpha
    bootstrap = args.boot
    rmodel = args.rmodel
    hapobj = args.hapobj
    wtobj = args.wtobj
    polyobj = args.polyobj
    depth = int(args.depth)
    

    print "Working in" + os.getcwd() + " on " + inputfile + " and " + reffile + " using " + treetype

    indata = InData()
    indata.load_ref_seq(reffile = reffile)
    indata.load_input_data(inputfile = inputfile)
    print "\n****************\nVARIANTS\n****************\n\n"
    for i in indata.variants.keys():
        print i + ": " + str(indata.variants[i])

    print "\n****************\nSEQUENCE\n****************\n\n"
    print indata.sequence

    print "\n****************\nTREE\n****************\n\n"
    phytree = PhyloTree()
    phytree.input_tree(treetype = treetype, alpha = alpha, bootstrap = bootstrap, rmodel = rmodel)
    print phytree.tree
    # Phylo.draw_ascii(phytree.tree)

    print "\n****************\nIMPUTATION\n****************\n\n"
    impute = Imputation(indata, phytree)
    impute.impute(depth)
    impute.imputedseq.sort()
    print impute.imputedseq

    # print "\n****************\nFALSE NEGATIVES\n****************\n\n"
    # stats = ChromStats(indata, hapobj, wtobj, polyobj)
    # print stats.haplogroups
    # print stats.weights
    # print stats.polys


else:
    print("IMPUTOR is being imported into another module. Not yet implemented.")