#!/usr/bin/python

""" Imputor - software for imputing missing and non-missing mutations in sequence data by the use
    of phylogenetic trees.
    
    Author - Matthew Jobin, Department of Anthropology, Santa Clara University and Department of Anthropology,
        University of California, Santa Cruz.
    """

import argparse
import glob
import multiprocessing
import os
import gzip
import random
import sys
import re
from argparse import RawTextHelpFormatter
import progressbar
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Phylo.TreeConstruction import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import operator


class InData(object):
    """Input data from FASTA or VCF file converted internally to sequence data.
        Prepares a dictionary of variants and a list of sequences for processing
        and imputation.
    """

    def __init__(self, minfile, mreffile):
        self.inputfile = minfile
        self.reffile = mreffile
        self.fullsequence = MultipleSeqAlignment([])  # Raw sequence input including non-segregating sites
        self.fullvariantset = set()
        self.fullvariants = {}
        self.sequence = MultipleSeqAlignment([])
        self.variantset = set()
        self.reflist = []
        self.variants = {}  # Dictionary of each sample and its variation from the reference sequence
        self.idxvariants = []
        self.filebase = None
        self.orig_vcf_pos = []  # Original listed positions of variants
        self.geninfo = {}  # Dict of list of dicts
        self.gqdict = {}
        self.addict = {}
        self.maxseqlength = 0
        self.pruned_to_full = []  # list where position in list is pruned location
        self.chroms = []

        # The eight mandatory columns of a VCF file. Here for clarity in functions below.
        self.vcf_chrom = 0
        self.vcf_pos = 1
        self.vcf_id = 2
        self.vcf_ref = 3
        self.vcf_alt = 4
        self.vcf_qual = 5
        self.vcf_filter = 6
        self.vcf_info = 7
        self.acgt = ['A', 'C', 'G', 'T']
        self.missing = ['-', 'N', '?']

        self.load_input_data()

        return

    def load_input_data(self):
        """ Load input data from a file.

        :return:
        """

        filecols = self.inputfile.split(".")
        self.filebase = filecols[0]

        if self.inputfile[-5:] == 'fasta':
            self.sequence = AlignIO.read(self.inputfile, 'fasta')
            if self.reffile:
                self.seqref = AlignIO.read(self.reffile, 'fasta')
                self.reflist = self.seqref[0].seq
            self.variants_from_sequence()
            self.prune_non_seg()
        elif self.inputfile[-3:] == 'vcf':
            file_data = open(self.inputfile, 'r')
            raw_data = []
            for file_line in file_data:
                if len(file_line.rstrip()) > 0:  # Strip blank lines
                    raw_data.append(file_line.rstrip())
            self.seq_from_variants(raw_data)
        elif self.inputfile[-6:] == 'vcf.gz':
            file_data = gzip.open(self.inputfile, 'r')
            raw_data = []
            for file_line in file_data:
                if len(file_line.rstrip()) > 0:  # Strip blank lines
                    raw_data.append(file_line.rstrip())

            self.seq_from_variants(raw_data)
        else:
            print "Input file must be either .fasta, .vcf, or .vcf.gz"
            exit()

        self.sequence.sort()

        if len(self.reflist) > 0:  # if there is a reference sequence, find variants
            for ref in self.reflist:
                self.idxvariants.append(list(ref))

            for seq in self.sequence:
                for i in xrange(len(seq.seq)):
                    if seq.seq[i] not in self.idxvariants[i]:
                        self.idxvariants[i].append(seq.seq[i])

        if verbose:
            if outtype == "vcf" and len(self.reflist) > 0:
                self.output_as_vcf()
            else:
                outseqfile = self.filebase
                if not seqonly:
                    outseqfile = outseqfile + "-indata"
                outseqfile = outseqfile + ".fasta"
                outfile = open(outseqfile, 'w')
                outseq = {}
                for seq in self.fullsequence:
                    outseq[seq.id] = str(seq.seq)
                for x in sorted(outseq.keys()):
                    outfile.write(">")
                    outfile.write(str(x))
                    outfile.write("\n")
                    outfile.write(outseq[x])
                    outfile.write("\n")
                outfile.close()
                if self.inputfile[-3:] == 'vcf':
                    outreffilename = self.filebase + "-indata-ref.fasta"
                    outreffile = open(outreffilename, 'w')
                    outreffile.write(">REF")
                    outreffile.write("\n")
                    outreffile.write("".join(self.reflist))
                    outreffile.write("\n")
                    outreffile.close()
        if rej:
            self.rej_infile(self.inputfile)

        print "Finished input."
        return

    def output_as_vcf(self):
        print "what"
        outseqfile = self.filebase
        if not seqonly:
            outseqfile = outseqfile + "-indata"
        outseqfile = outseqfile + ".vcf"
        outfile = open(outseqfile, 'w')
        outfile.write("##fileformat=VCFv4.1\n")
        outfile.write("##source=IMPUTORv1.0\n")
        outfile.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	")
        for seq in self.sequence:
            outfile.write(str(seq.name))
            outfile.write("\t")
        outfile.write("\n")
        for i in xrange(0, len(self.idxvariants)):
            if len(self.idxvariants[i]) > 1:
                outfile.write("0")
                outfile.write("\t")
                outfile.write(str(i))
                outfile.write("\t.\t")
                outfile.write(self.idxvariants[i][0])
                outfile.write("\t")
                for j in xrange(1, len(self.idxvariants[i])):
                    if j > 1:
                        outfile.write(",")
                    outfile.write(self.idxvariants[i][j])
                outfile.write("\t.\t.\t.\tGT\t")
                for seq in self.sequence:
                    outfile.write(str(self.idxvariants[i].index(seq.seq[i])))
                    outfile.write("\t")
                outfile.write("\n")


        outfile.close()


    def vcf_snp_prune(self, in_data=None):
        """ Returns data including only lines containing SNPs.

        :param in_data:
        :return:
        """
        snps_data = []
        print "\nPruning non-SNP entries..."
        bar = progressbar.ProgressBar(redirect_stdout=True)
        for i in bar(range(len(in_data))):
            file_line = in_data[i]
            cols = file_line.split('\t')

            # If the second character is a (meta-info line) or a blank line, ignore
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n') or len(file_line) < 1:
                continue
            cols[self.vcf_ref] = cols[self.vcf_ref].upper()
            cols[self.vcf_alt] = cols[self.vcf_alt].upper()
            if len(cols[self.vcf_ref]) > 1:  # if not a snp
                continue
            elif(cols[self.vcf_ref] not in self.acgt) and (cols[self.vcf_ref] not in self.missing): # if not a snp
                continue
            else:
                alt_alleles = cols[self.vcf_alt].split(",")  # List of ALT alleles for this row
                goodalt = True
                for allele_pos, chk_allele in enumerate(alt_alleles):  # Iterates through the alleles
                    if len(chk_allele) > 1:
                        goodalt = False
                    if chk_allele in self.missing:
                        alt_alleles[allele_pos] = "."
                if goodalt:
                    cols[self.vcf_alt] = ",".join(alt_alleles)
                    clean_file_line = "\t".join(cols)
                    snps_data.append(clean_file_line)
        return snps_data

    def variants_from_sequence(self):
        """ Returns a list of variants from a reference sequence.
            Intended for use with FASTA input, but will work with any AlignIO object or
            list of sequence data.

        :return:
        """
        print "\nSetting up sequence difference sets..."
        firstseq = self.sequence[0]
        diffsets = []
        for i in range(len(firstseq)):
            diffsets.append(set())
        bar = progressbar.ProgressBar(redirect_stdout=True)
        for i in bar(range(len(self.sequence))):
            seq_line = self.sequence[i]
            if len(seq_line) != len(firstseq):
                print "Error! A sequence line is a different size", len(seq_line), "than the first sequence!", len(
                    firstseq)
                print seq_line
                exit()
            for j in range(len(seq_line)):
                diffsets[j].add(seq_line[j])

        print "\nGenerating variants from sequence..."
        bar = progressbar.ProgressBar(redirect_stdout=True)
        for i in bar(range(len(self.sequence))):
            seq_line = self.sequence[i]
            if len(seq_line) != len(firstseq):
                print "Error! A sequence line is a different size", len(seq_line), "than the first sequence!", len(
                    firstseq)
                print seq_line
                return
            diffs = [i for i in xrange(len(firstseq)) if firstseq[i] != seq_line[i]]
            curdiffs = []
            for diff_pos in diffs:
                self.variantset.add(diff_pos)
                curdiffs.append(diff_pos)
            self.variants[seq_line.name] = curdiffs
        return

    def prune_non_seg(self):
        """ Strips out non-segregating sites from a sequence alignment.
            Uses self.variantset, which must be filled first.

        :return:
        """
        self.fullsequence = self.sequence  # First back up the original sequence
        self.fullvariantset = self.variantset
        self.fullvariants = self.variants
        self.sequence = MultipleSeqAlignment([])  # Blank the sequence to be worked on

        print "\nPruning non-segregating sites..."
        locs = []
        for curvar in self.variantset:
            locs.append(curvar)
        locs.sort()

        stripped = {}
        seqnames = []
        for seq in self.fullsequence:
            stripped[seq.name] = []
            seqnames.append(seq.name)

        for i in xrange(len(locs)):
            loc = locs[i]
            self.pruned_to_full.append(loc)
            seqbits = self.fullsequence[:, loc]
            name = 0
            for seqbit in seqbits:
                stripped[seqnames[name]].append(seqbit)
                name += 1

        for strip in stripped.keys():
            self.sequence.append(SeqRecord(Seq(''.join(stripped[strip])), name=strip, id=strip))

        self.variantset = set()
        self.variants = {}
        self.variants_from_sequence()  # Re-run on stripped sequence

    def seq_from_variants(self, raw_data=None):
        """ Create sequence using a list of variants.

        :param raw_data:
        :return:
        """
        print "\nCollecting genotype names..."
        genotype_names = []
        genotype_sequence = {}
        bar = progressbar.ProgressBar()
        for i in bar(range(len(raw_data))):
            file_line = raw_data[i]
            cols = file_line.split('\t')

            # Locate header line and read genotype names
            if cols[self.vcf_chrom] == '#CHROM':  # Header line of VCF file
                if cols[self.vcf_info + 1] == 'FORMAT':  # On header line, a FORMAT column next to the fixed columns?
                    genotype_names = cols[self.vcf_info + 2:]  # If so, remaining columns are the genotypes
                else:
                    print "Error. VCF file with no genotype. Cannot create sequence data."
                    return
        for genotype_name in genotype_names:  # Step through data lines, constructing list of variants
            self.variants[genotype_name] = []

        snps_data = self.vcf_snp_prune(raw_data)  # Ensure only SNPs are being processed
        var_count = 0
        print "\nGenerating sequence..."
        bar = progressbar.ProgressBar(redirect_stdout=True)
        for i in bar(range(len(snps_data))):
            file_line = snps_data[i]

            cols = file_line.split('\t')

            if int(cols[self.vcf_pos]) > self.maxseqlength:
                self.maxseqlength = int(cols[self.vcf_pos])
            self.orig_vcf_pos.append(cols[self.vcf_pos])
            self.reflist.append(cols[self.vcf_ref])
            self.chroms.append(cols[self.vcf_chrom])
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n' or cols[self.vcf_info + 1][:2] != "GT"):
                continue
            formatcols = cols[self.vcf_info + 1].split(":")
            indiv_genotypes = cols[self.vcf_info + 2:]  # Assumes all rows same length, as per VCF standard
            for position, indiv_genotype in enumerate(
                    indiv_genotypes):  # Iterates through that row of genotypes for this site
                genotypecols = indiv_genotype.split(":")
                assigned_alleles = re.split(  # VCF standard GT always first column
                    "[/|]+", genotypecols[0])  # Split genotype entry on either character phased or unphased

                changed_genotype_names = []
                for allele_pos, assigned_allele in enumerate(assigned_alleles):  # Iterates through the alleles
                    changed_genotype_name = genotype_names[position]
                    if len(assigned_alleles) > 1:  # Only append to genotype name if not haploid
                        changed_genotype_name += "-"
                        changed_genotype_name += str(allele_pos)
                    changed_genotype_names.append(changed_genotype_name)
                for changed_genotype_name in changed_genotype_names:
                    if changed_genotype_name not in genotype_sequence:
                        genotype_sequence[changed_genotype_name] = []

                alt_alleles = cols[self.vcf_alt].split(",")  # List of ALT alleles for this row
                for aa in range(len(alt_alleles)): #Convert other missing symbols to "N"
                    if alt_alleles[aa] == "." or alt_alleles[aa] == "-":
                        alt_alleles[aa] = "N"

                for allele_pos, assigned_allele in enumerate(assigned_alleles):  # Iterates through the alleles
                    if assigned_allele == "0":  # Assigned_allele will be 0 for REF and >0 for any ALT
                        genotype_sequence[changed_genotype_names[allele_pos]].append(cols[self.vcf_ref])
                    elif assigned_allele == ".":  # VCF format code for missing allele
                        genotype_sequence[changed_genotype_names[allele_pos]].append("N")
                    else:
                        genotype_sequence[changed_genotype_names[allele_pos]].append(alt_alleles[int(assigned_allele) - 1])
                        self.variantset.add(var_count)
                        if changed_genotype_names[allele_pos] in self.variants:  # Keys added to self.variants here
                            self.variants[changed_genotype_names[allele_pos]].append(var_count)
                        else:
                            self.variants[changed_genotype_names[allele_pos]] = []
                            self.variants[changed_genotype_names[allele_pos]].append(var_count)
                # Now dictionary of all genotype info
                for fi in range(len(formatcols)):
                    if fi < len(genotypecols):
                        for changed_genotype_name in changed_genotype_names:
                            finame = changed_genotype_name + "-" + str(var_count)
                            if formatcols[fi] == "GQ":
                                self.gqdict[finame] = genotypecols[fi]
                            if formatcols[fi] == "AD":
                                adcols = genotypecols[fi].split(",")
                                if len(adcols) == (len(alt_alleles) + 1):
                                    colc = 0
                                    ad = ""
                                    for adcol in adcols:
                                        if colc == 0:
                                            ad += cols[self.vcf_ref]
                                        else:
                                            ad += alt_alleles[colc - 1]
                                        ad += "-"
                                        ad += adcol
                                        colc += 1
                                        if colc < len(adcols):
                                            ad += ","
                                    self.addict[finame] = ad
            var_count += 1

        for geno in genotype_sequence.keys():
            genotype_sequence[geno] = ''.join(genotype_sequence[geno])

        self.sequence = MultipleSeqAlignment([])  # Blank the sequence to be worked on
        for geno in genotype_sequence.keys():
            self.sequence.append(SeqRecord(Seq(''.join(genotype_sequence[geno])), name=geno, id=geno))

        self.fullsequence = self.sequence
        self.fullvariantset = self.variantset
        self.fullvariants = self.variants

    def rej_infile(self, linfile=None):
        """ Output a REJECTOR2 input file

        :param linfile:
        :return:
        """

        rejfilename = self.filebase + "-rej.txt"

        rejfile = open(rejfilename, 'w')
        rejfile.write("/--Data\n")
        rejfile.write("Vnaught 0\n\n")
        rejfile.write("Loci\tSNP\n")
        rejfile.write("Ancestral\t-1\n")
        rejfile.write("RecombRt\t0\n")
        rejfile.write("NumLoci\t")
        rejfile.write(str(len(self.sequence[0].seq)))
        rejfile.write("\n")
        rejfile.write("Length\t1\n")
        rejfile.write("\n")
        rejfile.write("\n")

        rejfile.write("Tag\t")
        rejfile.write("Population\n")

        outseq = {}
        for seq in self.sequence:
            outseq[seq.id] = str(seq.seq)
        for x in sorted(outseq.keys()):
            rejfile.write(str(x))
            rejfile.write("\t")
            rejfile.write("X")
            rejfile.write("\t")
            for y in list(outseq[x]):
                rejfile.write(y)
                rejfile.write("\t")
            rejfile.write("\n")

        rejfile.close()


class PhyloTree(object):
    """A phylogenetic tree either input from phyloxml format or constructed from sequence.


    """

    def __init__(self, mindata, mtreetype, mtimestamp, malpha, mrmodel, mstarttreename, mmaxthreads=4, mmaxdepth=3,
                 mmaxheight=2, mmaxn=5):

        self.indata = mindata
        self.treetype = mtreetype
        self.impname = "imp" + str(mtimestamp)
        self.starttreename = None
        self.rmodel = mrmodel
        self.starttreename = mstarttreename
        self.maxthreads = mmaxthreads
        self.maxdepth = mmaxdepth
        self.maxheight = mmaxheight
        self.maxn = mmaxn
        self.tree = None  # Phylogenetic tree to be loaded or constructed from data. Newick format.
        self.treeparents = {}
        self.starttree = None  # Phylogenetic tree used as starting tree in RAxML
        self.btrees = []  # Bootstrap replicates of trees. Newick format.
        self.btreeparents = []
        self.raxmlalgs = {'a', 'd', 'o'}
        self.alpha = malpha

        self.input_tree()

    def input_tree(self):
        """
        Takes input tree file or sequence data.

        :return:
        """

        if self.starttreename:
            if self.starttreename[-3:] == 'xml':
                self.starttree = Phylo.read(self.starttreename, "phyloxml")
            elif self.starttreename[-6:] == 'newick':
                self.starttree = Phylo.read(self.starttreename, "newick")

        print "Generating phylogenetic tree..."

        if self.treetype[-3:] == 'xml':
            self.tree = Phylo.read(self.treetype, "phyloxml")
        elif self.treetype[-3:] == 'nwk':
            self.tree = Phylo.read(self.treetype, "newick")
        elif self.treetype == 'pars':
            self.parsimony_tree()
        elif self.treetype == 'PhyML':
            self.phyml_tree()
        else:
            self.raxml_tree()

        self.tree.collapse_all(lambda c: c.branch_length <= 0.0)
        self.treeparents = self.all_parents(self.tree)
        for btree in self.btrees:
            btree.collapse_all(lambda c: c.branch_length <= 0.0)
            self.btreeparents.append(self.all_parents(btree))

    def output_tree(self, outputtreetype):
        """ Outputs tree to file.

        :param outputtreetype: type of tree to output
        :return:
        """
        if outputtreetype == 'phyloxml':
            outfile = self.indata.filebase + "-outtree.xml"
            Phylo.write(self.tree, outfile, "phyloxml")
        elif outputtreetype == 'nexus':
            outfile = self.indata.filebase + "-outtree.nexus"
            Phylo.write(self.tree, outfile, "nexus")
        else:  # Default newick
            outfile = self.indata.filebase + "-outtree.nwk"
            Phylo.write(self.tree, outfile, "newick")

    @staticmethod
    def all_parents(tree):
        """ Find all parents of all children in the tree.

        :param tree: Input tree.
        :return: List of parents for each clade.
        """
        parents = {}
        for clade in tree.find_clades(order='level'):
            for child in clade:
                parents[child] = clade
        return parents

    def collect_kids(self, clade, kids, depth, locmaxdepth):
        """ Recursive function for collecting all children to specified depth

        :param locmaxdepth: local instance of maximum depth
        :param clade: current node of tree
        :param kids: list of all collected children
        :param depth: current depth of search
        :return:
        """
        if depth < locmaxdepth:
            for child in clade:
                if child.name:
                    kids.append(child)
            for child in clade:
                self.collect_kids(child, kids, depth + 1, locmaxdepth)
        return kids

    def collect_all_kids(self, clade, kids):
        """ Recursive function for collecting all children of clade

        :param clade: current node of tree
        :param kids: lost of all collected children
        :return:
        """
        for child in clade:
            if child.name:
                kids.append(child)
            self.collect_all_kids(child, kids)
        return kids

    def neighbors_by_distance(self, term, terms, tsize):
        """ Function for collecting all children of clade

        :param term: current node of tree
        :param terms: lost of all collected children
        :param tsize: threshold number of neighbors
        :return: set of closest neighbors
        """
        neighbors = set()
        ndist = {}
        for neb in terms:
            if neb is not term:
                ndist[neb] = self.tree.distance(term, neb)
        sorted_neb = sorted(ndist.items(), key=operator.itemgetter(1))
        for i in xrange(len(sorted_neb)):
            if i >= tsize:
                break
            if i > 0:
                if sorted_neb[i][1] > 0.0 and sorted_neb[i - 1][1] > 0.0:
                    if sorted_neb[i][1] / sorted_neb[i - 1][1] > maxjump:
                        break
            neighbors.add(sorted_neb[i][0])
        return neighbors

    def neighbors_by_mono(self, term, ctree, parents, tsize):
        """ Collect neighbors of minimum sized monophyletic clade including term.

        :param term: current node of tree
        :param ctree: current tree
        :param parents: list of parents of nodes
        :param tsize: threshold number of neighbors
        :return: set of closest neighbors
        """
        neighbors = set()
        monn = set()
        monn.add(term)
        curnode = term
        while len(neighbors) < tsize:
            if curnode not in parents:  # will not go past the root
                break
            curparent = parents[curnode]
            allkids = self.collect_all_kids(curparent, [])
            for kid in allkids:
                if kid is not term:
                    neighbors.add(kid)
                    monn.add(kid)
                    if len(neighbors) >= tsize and ctree.is_monophyletic(monn):
                        return neighbors
                    if len(neighbors) > self.maxn:
                        return set()
            curnode = curparent
        return set()

    def neighbors_by_hops(self, term, ctree, parents, tsize):
        """ Collect neighbors up to a maximum number of hops along tree.

        :param term: current node of tree
        :param ctree: current tree
        :param parents: list of parents of nodes
        :param tsize: threshold number of neighbors
        :return: set of closest neighbors
        """
        workneighbors = set()
        neighbors = set()
        monn = set()
        monn.add(term)

        height = 0
        while len(workneighbors) <= tsize:
            curnode = term
            for i in xrange(height):
                if curnode not in parents:  # will not go past the root
                    break
                curnode = parents[curnode]
                allkids = self.collect_kids(curnode, [], 0, height + 1)
                for kid in allkids:
                    if kid is not term:
                        workneighbors.add(kid)
            height += 1
        ndist = {}
        for neb in workneighbors:
            if len(ctree.trace(term, neb)) <= maxhops:
                ndist[neb] = ctree.distance(term, neb)
        sorted_neb = sorted(ndist.items(), key=operator.itemgetter(1))
        for i in xrange(len(sorted_neb)):
            if i >= tsize:
                break
            monn.add(sorted_neb[i][0])
            neighbors.add(sorted_neb[i][0])
        return neighbors

    def neighbors_by_rootward(self, term, parents, height, tsize, ctree):
        """

        :param ctree:
        :param term: current node of tree
        :param parents: list of parents of nodes
        :param height:
        :param tsize: threshold number of neighbors
        :return: set of closest neighbors
        """
        neighbors = set()
        curnode = term
        # Collect closest neighbors on tree

        # print "TARGET: ", str(term)
        while len(neighbors) < tsize:
            if curnode not in parents:  # will not go past the root
                break
            if height > self.maxheight:
                break
            curparent = parents[curnode]

            # print "\tPARENT: ", str(curparent)
            allkids = self.collect_kids(curparent, [], 0, self.maxdepth)

            # ndist = {}
            # for neb in allkids:
            #     if len(ctree.trace(term, neb)) <= maxhops:
            #         ndist[neb] = ctree.distance(term, neb)
            # sorted_neb = sorted(ndist.items(), key=operator.itemgetter(1))

            addedkids = 0
            for kid in allkids:
                if kid is not term:
                    # print "\t\tADD KID: ", str(kid[0])
                    addedkids += 1
                    neighbors.add(kid)
                    if len(neighbors) >= tsize:
                        break
            if orphanchk and addedkids == 0:
                break
            curnode = curparent
            height += 1
        return neighbors

    def parsimony_tree(self):
        """ Constructs a tree via maximum parsimony using Biopython's ParsimonyTreeConstructor.

        """
        print "Generating maximum parsimony tree.."
        if runs > 0 or boot > 0:
            print "ERROR: Bootstrap and multiple runs not compatible with -tree pars option."
            exit()
        cpus = multiprocessing.cpu_count()
        if cpus > maxthreads:
            cpus = maxthreads
        # Erase RaXML intermediate files from previous runs
        raxml_glob = glob.glob('RAxML_*')
        for delfile in raxml_glob:
            os.remove(delfile)

        # Output sequence to a temp FASTA file
        tempfastafile = self.indata.filebase + self.impname + "_fastatmp.fasta"
        reducedtempfastafile = self.indata.filebase + self.impname + "_fastatmp.fasta.reduced"
        AlignIO.write(self.indata.sequence, tempfastafile, "fasta")

        raxml_args = {"sequences": tempfastafile, "model": self.rmodel, "name": self.impname,
                      "parsimony_seed": rng.randint(0, sys.maxint), "threads": cpus, "parsimony": True,
                      "algorithm": "d"}

        raxmlstarttreename = "RAxML_" + self.impname + "_starttree.newick"
        if self.starttree:
            Phylo.write(self.starttree, raxmlstarttreename, "newick")
            raxml_args["starting_tree"] = raxmlstarttreename

        if exlocal:
            raxml_cline = RaxmlCommandline(cmd='./raxmlHPC', **raxml_args)
        else:
            raxml_cline = RaxmlCommandline(**raxml_args)

        print "Invoking RAxML with ", raxml_cline

        out_log, err_log = raxml_cline()
        if verbose:
            print err_log
            print out_log
        raxmlparstreename = "RAxML_parsimonyTree." + self.impname
        self.tree = Phylo.read(raxmlparstreename, "newick")

        # Erase RaXML intermediate files
        if not verbose:
            raxml_glob = glob.glob('RAxML_*')
            for delfile in raxml_glob:
                os.remove(delfile)

        try:
            os.remove(tempfastafile)
        except OSError:
            pass

        try:
            os.remove(reducedtempfastafile)
        except OSError:
            pass

    def raxml_tree(self):
        """ Constructs a tree via maximum likelihood by invoking external software RAxML.
            See docs for RAxML installation and setup.

        """
        cpus = multiprocessing.cpu_count()
        if cpus > maxthreads:
            cpus = maxthreads

        # Erase RaXML intermediate files from previous runs
        raxml_glob = glob.glob('RAxML_*')
        for delfile in raxml_glob:
            os.remove(delfile)

        # Output sequence to a temp FASTA file
        tempfastafile = self.indata.filebase + self.impname + "_fastatmp.fasta"
        reducedtempfastafile = self.indata.filebase + self.impname + "_fastatmp.fasta.reduced"
        AlignIO.write(self.indata.sequence, tempfastafile, "fasta")

        raxml_args = {"sequences": tempfastafile, "model": self.rmodel, "name": self.impname,
                      "parsimony_seed": rng.randint(0, sys.maxint), "threads": cpus}

        if boot > 0:
            raxml_args["num_replicates"] = boot
            raxml_args['rapid_bootstrap_seed'] = rng.randint(0, sys.maxint)
            raxml_args['algorithm'] = "a"
        if runs > 0:
            print "Multi run"
            raxml_args["num_replicates"] = runs
            raxml_args['algorithm'] = "d"
        else:
            raxml_args['algorithm'] = "d"

        raxmlstarttreename = "RAxML_" + self.impname + "_starttree.newick"
        if self.starttree:
            Phylo.write(self.starttree, raxmlstarttreename, "newick")
            raxml_args["starting_tree"] = raxmlstarttreename

        if exlocal:
            raxml_cline = RaxmlCommandline(cmd='./raxmlHPC', **raxml_args)
        else:
            raxml_cline = RaxmlCommandline(**raxml_args)

        print "Invoking RAxML with ", raxml_cline

        out_log, err_log = raxml_cline()
        if verbose:
            print err_log
            print out_log

        if runs > 0:
            for i in range(runs):
                raxmlrunname = "RAxML_result." + self.impname + ".RUN." + str(i)
                self.btrees.append(Phylo.read(raxmlrunname, "newick"))
        if boot > 0:
            raxmlbsname = "RAxML_bootstrap." + self.impname
            self.btrees = list(Phylo.parse(raxmlbsname, "newick"))
        else:
            raxmlbesttreename = "RAxML_bestTree." + self.impname
            self.tree = Phylo.read(raxmlbesttreename, "newick")

        # Erase RaXML intermediate files
        raxml_glob = glob.glob('RAxML_*')
        for delfile in raxml_glob:
            os.remove(delfile)

        try:
            os.remove(tempfastafile)
        except OSError:
            pass

        try:
            os.remove(reducedtempfastafile)
        except OSError:
            pass

    def phyml_tree(self):
        """ Constructs a tree via maximum likelihood by invoking external software PhyML.
            See docs for PhyML installation and setup.

        """
        print "Invoking PhyML..."
        if runs > 0 or boot > 0:
            print "ERROR: Bootstrap and multiple runs not yet implemented for PhyML."
            print "Try using RAxML."
            exit()
        # Output sequence to a temp FASTA file
        tempfastafile = self.indata.filebase + "_" + self.impname + "_fastatmp.fasta"
        AlignIO.write(self.indata.sequence, tempfastafile, "fasta")
        tempphyfile = self.indata.filebase + "_" + self.impname + "_phytmp.phy"
        AlignIO.convert(tempfastafile, "fasta", tempphyfile, "phylip-relaxed")

        phyml_args = {"input": tempphyfile, "alpha": "e"}
        phystarttreename = "PhyML_imp", self.impname, "starttree.newick"
        if self.starttree:
            Phylo.write(self.starttree, phystarttreename, "newick")
            phyml_args["input_tree"] = phystarttreename

        if exlocal:
            cmdline = PhymlCommandline(cmd='./PhyML', **phyml_args)
        else:
            cmdline = PhymlCommandline(**phyml_args)

        print "Commandline for PhyML: " + str(cmdline)
        out_log, err_log = cmdline()
        if verbose:
            print err_log
            print out_log
        phytreefile = tempphyfile + "_phyml_tree.txt"
        self.tree = Phylo.read(phytreefile, "newick")
        if not verbose:
            phyml_globname = self.indata.filebase + "_" + self.impname + "*"
            phyml_glob = glob.glob(phyml_globname)
            for delfile in phyml_glob:
                os.remove(delfile)


class Imputation(object):
    """Imputation of mutations given input data and a phylogenetic tree.
    
            Keyword arguments:
            indata -- Input data.
            phytree -- Phylogenetic tree
            mutrate -- Mutation rate.
    """

    def __init__(self, mindata, mphytree, mmutrate=8.71e-10, mtstv=2.0):
        self.workseq = {}
        self.imputedseq = MultipleSeqAlignment([])
        self.indata = mindata
        self.phytree = mphytree
        self.mu = mmutrate
        self.tstv = mtstv
        self.imputelist = []
        self.multi = multi
        self.acgt = {'A', 'C', 'G', 'T'}
        self.missing = {'.', '-', 'N'}
        self.indivimputes = {}
        self.neighbors = {}
        self.missinglist = {}
        self.newvariants = []
        # self.backmutchks = []
        self.bootreps = {}

        for seq in indata.sequence:
            self.workseq[seq.name] = list(str(seq.seq))  # Begin with output sequence matching input
            self.indivimputes[seq.name] = []

    def impute(self):
        """ Sets up imputation function for all terminal nodes.

        """

        if ncollect == 'hops':
            print "HOPS"
        elif ncollect == 'distance':
            print "DISTANCE"
        else:
            print "ROOTWARD"

        terms = self.phytree.tree.get_terminals()  # Get all internal nodes on tree. These are the ones with samples.
        if boot > 0 or runs > 0:  # Bootstrap replicates or multiple runs
            for i in range(passes):
                print "\nPass", i
                bpar = iter(phytree.btreeparents)
                bar = progressbar.ProgressBar(redirect_stdout=True)
                for j in bar(range(len(self.phytree.btrees))):
                    # for btree in self.phytree.btrees:
                    btree = self.phytree.btrees[j]
                    terms = btree.get_terminals()
                    random.shuffle(terms)
                    bparents = next(bpar)
                    for term in terms:
                        if ncollect == 'hops':
                            nneighbors = phytree.neighbors_by_hops(term, btree, bparents,
                                                                   nsize)
                            mneighbors = phytree.neighbors_by_hops(term, btree, bparents,
                                                                   msize)
                        elif ncollect == 'distance':
                            nneighbors = phytree.neighbors_by_distance(term, terms, nsize)
                            mneighbors = phytree.neighbors_by_distance(term, terms, msize)
                        elif ncollect == 'mono':
                            nneighbors = phytree.neighbors_by_mono(term, btree, bparents,
                                                                   nsize)
                            mneighbors = phytree.neighbors_by_mono(term, btree, bparents,
                                                                   msize)
                        else:
                            nneighbors = phytree.neighbors_by_rootward(term, bparents, 0, nsize, btree)
                            mneighbors = phytree.neighbors_by_rootward(term, bparents, 0, msize, btree)
                        self.impute_bootstrap(term, bparents, str(i + 1), mneighbors, nneighbors)
                for bootrep in self.bootreps:
                    newimpute = bootrep.split(".")
                    impratio = float(self.bootreps[bootrep][0]) / float(self.bootreps[bootrep][1])
                    newimpute.append(str(impratio))
                    if verbose:
                        if impratio > threshold:
                            newimpute.append("T")
                            self.workseq[newimpute[0]][int(newimpute[1])] = newimpute[3]
                            self.indivimputes[newimpute[0]].append(newimpute[1])
                        else:
                            newimpute.append("F")
                        self.imputelist.append(newimpute)
                    else:
                        if impratio > threshold:
                            newimpute.append("T")
                            self.imputelist.append(newimpute)
                            self.workseq[newimpute[0]][int(newimpute[1])] = newimpute[3]
                            self.indivimputes[newimpute[0]].append(newimpute[1])

        else:
            for i in range(passes):
                print "\nPass", i
                random.shuffle(terms)
                bar = progressbar.ProgressBar(redirect_stdout=True)
                for p in bar(range(len(terms))):
                    term = terms[p]
                    # for term in terms:
                    if ncollect == 'hops':
                        nneighbors = phytree.neighbors_by_hops(term, self.phytree.tree, self.phytree.treeparents, nsize)
                        mneighbors = phytree.neighbors_by_hops(term, self.phytree.tree, self.phytree.treeparents, msize)
                    elif ncollect == 'distance':
                        nneighbors = phytree.neighbors_by_distance(term, terms, nsize)
                        mneighbors = phytree.neighbors_by_distance(term, terms, msize)
                    elif ncollect == 'mono':
                        nneighbors = phytree.neighbors_by_mono(term, self.phytree.tree, self.phytree.treeparents,
                                                               nsize)
                        mneighbors = phytree.neighbors_by_mono(term, self.phytree.tree, self.phytree.treeparents,
                                                               msize)
                    else:
                        nneighbors = phytree.neighbors_by_rootward(term, self.phytree.treeparents, 0, nsize,
                                                                   self.phytree.tree)
                        mneighbors = phytree.neighbors_by_rootward(term, self.phytree.treeparents, 0, msize,
                                                                   self.phytree.tree)

                    self.impute_threshold(term, self.phytree.treeparents, str(i + 1),
                                          mneighbors, nneighbors)
                for newimpute in self.imputelist:
                    if newimpute[6] == "T":
                        self.indivimputes[newimpute[0]].append(newimpute[1])
                        self.workseq[newimpute[0]][newimpute[1]] = newimpute[3]
        self.process_imputed()

    def impute_threshold(self, term, parents, thispass, mneighbors, nneighbors):
        """ Imputation on best tree.

        :param parents: list of parents of all nodes in current tree
        :param term: terminal node
        :param thispass: number of current imputation pass
        :param mneighbors: neighbors collected with the msize threshold
        :param nneighbors: neighbors collected with the nsize threshold
        :return:
        """
        termname = str(term)
        for curvar in self.indata.variantset:
            if self.workseq[termname][curvar] in self.missing:
                neighbors = mneighbors
            else:
                neighbors = nneighbors
            nbs = []
            for nb in neighbors:
                nbs.append(str(nb))
            nbname = str(term) + "-" + str(curvar)
            self.neighbors[nbname] = nbs
            newimpute = self.detect_by_parsimony(term, curvar, parents, neighbors, thispass)
            if verbose:
                self.imputelist.append(newimpute)
            else:
                if newimpute[6] == "T":
                    self.imputelist.append(newimpute)

    def impute_bootstrap(self, term, bparents, thispass, mneighbors, nneighbors):
        """ Imputation with bootstrap replicates or multiple runs on a single file.

        :param term: terminal node
        :param bparents: list of parents of all nodes in current tree
        :param thispass: number of current imputation pass
        :param mneighbors: neighbors collected with the msize threshold
        :param nneighbors: neighbors collected with the nsize threshold
        :return:
        """
        termname = str(term)
        for curvar in self.indata.variantset:
            if self.workseq[termname][curvar] in self.missing:
                neighbors = mneighbors
            else:
                neighbors = nneighbors
            newimpute = self.detect_by_parsimony(term, curvar, bparents, neighbors, thispass)
            bnewimpute = [newimpute[0], str(newimpute[1]), newimpute[2], newimpute[3], newimpute[5]]
            bootfront = ".".join(bnewimpute)
            newimpute[1] = str(indata.pruned_to_full[newimpute[1]])
            if bootfront not in self.bootreps:
                self.bootreps[bootfront] = [0, 0]
            if newimpute[6] == "T":
                self.bootreps[bootfront][0] = self.bootreps[bootfront][0] + 1
            self.bootreps[bootfront][1] = self.bootreps[bootfront][1] + 1

    def process_imputed(self):
        """ Restore original locations of sequence when non-segregating sites removed.
            Make imputations to final sequence.

        """
        print "\nProcessing imputed sequences..."
        locs = []
        for curvar in indata.fullvariantset:
            locs.append(curvar)
        locs.sort()

        bar = progressbar.ProgressBar(redirect_stdout=True)
        for p in bar(range(len(indata.fullsequence))):
            fullseq = indata.fullsequence[p]
            tmpseq = list(fullseq)
            segseq = self.workseq[fullseq.id]

            if len(segseq) == len(locs):
                for site, loc in itertools.izip(segseq, locs):  # Relies on original sequence of non-seg sites
                    tmpseq[loc] = site

            seqrec = SeqRecord(Seq("".join(tmpseq)), id=fullseq.id, name=fullseq.id)
            self.imputedseq.append(seqrec)

            if len(indata.reflist) > 0:  # if there is a reference sequence, find variants
                for ref in indata.reflist:
                    self.newvariants.append(list(ref))

                for seq in self.imputedseq:
                    for i in xrange(len(seq.seq)):
                        if seq.seq[i] not in self.newvariants[i]:
                            self.newvariants[i].append(seq.seq[i])



        self.imputedseq.sort()

    def transversionchk(self, a, b):
        """ Check chance whether mutation is a translation or a transversion

        :param a: One site.
        :param b: another site.
        :return: Translation or Transversion rate.
        """
        if a == "A":
            if b == "G":
                return self.tstv
        if a == "G":
            if b == "A":
                return self.tstv
        if a == "C":
            if b == "T":
                return self.tstv
        if a == "T":
            if b == "C":
                return self.tstv
        if a in ("A", "C", "G", "T") and b in ("A", "C", "G", "T"):
            return 1 / self.tstv
        return -1

    def backmutchk(self, term, parents, allneighbours, curvar, origseq, neighborseq):
        """ Check whether there is a reversion somewhere in the tree outside chosen neighbors of target
        :param term: terminal node
        :param parents: list of parents of all nodes in tree
        :param allneighbours: set of closest neighbors chosen by user-selected method
        :param curvar: current variant
        :param origseq: Original state of site.
        :param neighborseq: State of neighbors' sites.
        :return:
        """
        backmut = False
        origneighbors = []
        for nb in allneighbours:
            origneighbors.append(str(nb))

        # print "\n", term
        curnode = term

        while curnode in parents:
            if backmut:
                return True
            curparent = parents[curnode]
            # nextparent = parents[curparent]
            empty = []
            allkids = phytree.collect_all_kids(curparent, empty)
            # print "ak ", len(allkids)
            for kid in allkids:
                kidseq = self.workseq[str(kid)][curvar]
                if kidseq in self.missing:
                    continue
                if kid not in allneighbours and kidseq == origseq:
                    # self.backmutchks.append(
                    #     [str(term), str(indata.pruned_to_full[curvar]), origseq, self.workseq[str(term)][curvar],
                    #      str(origneighbors), neighborseq,
                    #      str(kid), kidseq, "T"])
                    # print "Found BM neighbors: ", len(allneighbours)
                    return True
                allneighbours.add(kid)
            curnode = curparent
        # self.backmutchks.append(
        #     [str(term), str(str(indata.pruned_to_full[curvar])), origseq, self.workseq[str(term)][curvar],
        #      str(origneighbors), neighborseq, "N/A", "N/A", "F"])
        # print term, " NO BM neighbors: ", len(allneighbours)
        return False

    def detect_by_parsimony(self, term, curvar, parents, neighbors, thispass):
        """ Determine whether imputation should occur. Goes randomly through nods and site by site.

        :param term: terminal node
        :param curvar: current variant
        :param parents: list of parents of all nodes in tree
        :param neighbors: set of closest neighbors chosen by user-selected method
        :param thispass: which pass of the software over the data
        :return: list defining whether to impute and information about the site and term
        """

        termname = str(term)
        finame = termname + "-" + str(curvar)

        bmneighbors = set()
        bmneighbors.add(term)
        nearest = set()
        orig = self.workseq[termname][curvar]

        for neighbor in neighbors:
            bmneighbors.add(neighbor)
            nearest.add(self.workseq[str(neighbor)][curvar])

        if len(nearest) > 1:  # Cannot allow non-matching ever to impute sequence, thus this is first check
            if mnm and len(nearest) == msize and "N" in nearest:
                nearest.remove("N")  # Strip the one N, allowing imputation
            else:
                return [termname, curvar, orig, ",".join(nearest), "Neighbors Non-matching", thispass, "F"]
        if len(nearest) < 1:
            return [termname, curvar, orig, ".", "No Neighbors", thispass, "F"]
        only = nearest.pop()
        if orig in self.missing:
            if len(neighbors) < msize:
                return [termname, curvar, orig, only, "Not Enough Neighbors", thispass, "F"]
            elif only in self.missing:
                return [termname, curvar, orig, only, "Neighbors All Missing", thispass, "F"]
            else:
                return [termname, curvar, orig, only, "Imputed Missing", thispass, "T"]
        else:
            if len(neighbors) < nsize:
                return [termname, curvar, orig, only, "Not enough neighbors", thispass, "F"]
            if only in self.missing:
                return [termname, curvar, orig, only, "Neighbors All Missing", thispass, "F"]
            if orig != only:  # If target sequence does not match matching non-missing neighbors
                adtot = 0
                thisad = 0.0
                otherad = 0.0
                if finame in indata.addict:
                    ads = indata.addict[finame].split(",")
                    for ad in ads:
                        a = ad.split("-")
                        adtot += int(a[1])
                        if (a[0]) == orig:
                            thisad += float(a[1])
                        else:
                            otherad += float(a[1])
                    if adtot < mincoverage:
                        return [termname, curvar, orig, only, "Imputed By Min. Cov.", thispass, "T"]
                    if otherad > 0.0:
                        if thisad / otherad <= adthresh:
                            return [termname, curvar, orig, only, "Imputed By AD", thispass, "T"]
                if finame in indata.gqdict:
                    if int(indata.gqdict.get(finame)) < genoqual:
                        return [termname, curvar, orig, only, "Imputed By GQ", thispass, "T"]
                if nobackmutchk:
                    return [termname, curvar, orig, only, "Imputed Non-missing", thispass, "T"]
                else:
                    if self.backmutchk(term, parents, bmneighbors, curvar, orig, only):
                        return [termname, curvar, orig, only, "Imputed Non-missing", thispass, "T"]
                    else:
                        return [termname, curvar, orig, only, "No Reversion", thispass, "F"]
            else:
                return [termname, curvar, orig, only, "Matches Neighbors", thispass, "F"]

    def output_imputed(self, limpout):
        """ Output imputed sequence and auxilliary files.

        :param limpout: switch to determine whether to print information file about imputations
        :return:
        """
        for imputed in self.imputelist:
            if indata.orig_vcf_pos:
                imputed[1] = str(indata.orig_vcf_pos[int(imputed[1])])
            else:
                imputed[1] = str(imputed[1])

        if verbose:
            if len(self.imputelist) > 0:
                print "Imputed Mutations"
                print "SUBJECTID | VAR | FROM | TO | TYPE | IMPUTED | PASS"
                for imputed in sorted(self.imputelist):
                    print " | ".join(imputed)
                print "\n"
            print impute.imputedseq

        if limpout:
            impoutfilename = indata.filebase + "-impout.txt"
            impoutfile = open(impoutfilename, 'w')
            if boot > 0 or runs > 0:
                impoutfile.write("SUBJECTID\t VAR\t FROM\t TO\t PASS\tRATIO\tIMPUTED\n")
            else:
                impoutfile.write("SUBJECTID\t VAR\t FROM\t TO\t TYPE\tPASS\tIMPUTED\n")
            for imputed in self.imputelist:
                impoutfile.write("\t".join(imputed))
                impoutfile.write("\n")
            impoutfile.close()

        indivoutfilename = indata.filebase + "-indivout.txt"
        indivoutfile = open(indivoutfilename, 'w')
        indivoutfile.write("SUBJECTID\tNUM\tVARS\n")
        for indiv in sorted(self.indivimputes.keys()):
            indivoutfile.write(indiv)
            indivoutfile.write("\t")
            indivoutfile.write(str(len(self.indivimputes[indiv])))
            indivoutfile.write("\t")
            for indivar in self.indivimputes[indiv]:
                indivoutfile.write(str(indivar))
                indivoutfile.write(",")
            indivoutfile.write("\n")
        indivoutfile.close()

        if outtype == "vcf" and len(indata.reflist) > 0:
            outseqfile = indata.filebase + "-out.vcf"
            outfile = open(outseqfile, 'w')
            outfile.write("##fileformat=VCFv4.1\n")
            outfile.write("##source=IMPUTORv1.0\n")
            outfile.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	")
            for seq in self.imputedseq:
                outfile.write(str(seq.name))
                outfile.write("\t")
            outfile.write("\n")
            for i in xrange(0, len(self.newvariants)):
                if len(self.newvariants[i]) > 1:
                    outfile.write(indata.chroms[i])
                    outfile.write("\t")
                    outfile.write(indata.orig_vcf_pos[i])
                    outfile.write("\t.\t")
                    outfile.write(self.newvariants[i][0])
                    outfile.write("\t")
                    for j in xrange(1, len(self.newvariants[i])):
                        if j > 1:
                            outfile.write(",")
                        outfile.write(self.newvariants[i][j])
                    outfile.write("\t.\t.\t.\tGT\t")
                    for seq in self.imputedseq:
                        outfile.write(str(self.newvariants[i].index(seq.seq[i])))
                        outfile.write("\t")
                    outfile.write("\n")

        else:  # default to fasta
            outseqfile = indata.filebase + "-seqout.fasta"
            outfile = open(outseqfile, 'w')
            outseq = {}
            for seq in self.imputedseq:
                outseq[seq.id] = str(seq.seq)
            for x in sorted(outseq.keys()):
                outfile.write(">")
                outfile.write(str(x))
                outfile.write("\n")
                outfile.write(outseq[x])
                outfile.write("\n")
            outfile.close()

            # bmfile = open("backmut.txt", 'w')
            # bmfile.write("term\tvar\torigseq\torgseqchk\torigneighbors\tneighborseq\tbmkid\tkidseq\t\T/F\n")
            # for bmchk in self.backmutchks:
            #     bmfile.write("\t".join(bmchk))
            #     bmfile.write("\n")
            #
            # nbfile = open("neighbors.txt", 'w')
            # for nb in self.neighbors.keys():
            #     nbfile.write(str(nb))
            #     nbfile.write("\t:\t")
            #     for nbb in self.neighbors[nb]:
            #         nbfile.write(str(nbb))
            #         nbfile.write("\t")
            #     nbfile.write("\n")


if __name__ == "__main__":

    print "\n****************\nIMPUTOR\n****************\n"

    parser = argparse.ArgumentParser(description="IMPUTOR is a program for the phylogeny-aware imputation of: \n\n\t"
                                                 "mutations.\n\t"
                                                 "- ", formatter_class=RawTextHelpFormatter)

    parser.add_argument('-file', metavar='<file>', help='input file: FASTA or VCF', required=True)
    parser.add_argument('-reffile', metavar='<reffile>', help='reference sequence for FASTA  input file: FASTA')
    parser.add_argument('-tree', metavar='<tree>', help='tree type; <treefilename.xml>, <treefilename.nwk>, pars, RAxML, PhyML',
                        default='RAxML')
    parser.add_argument('-outtree', metavar='<outtree>', help='Output format for tree', default='newick')
    parser.add_argument('-alpha', metavar='<alpha>', help='Value of gamma shape parameter.', default='e')
    parser.add_argument('-rmodel', metavar='<rmodel>', help='Model type for RaXML.', default='GTRCAT')
    parser.add_argument('-maxheight', metavar='<maxheight>',
                        help='Height of search toward root for collecting neighbors.',
                        default=3)
    parser.add_argument('-maxdepth', metavar='<maxdepth>',
                        help='Depth of search down descendent nodes for collecting neighbors.', default=3)
    parser.add_argument('-nsize', metavar='<nsize>',
                        help='Number of neighbors that must be identical in order to impute.',
                        default=3)
    parser.add_argument('-msize', metavar='<msize>',
                        help='Number of neighbors that must be identical in order to impute for missing data.',
                        default=2)
    parser.add_argument('-mutrate', metavar='<mutrate>', help='Mutation rate.', default=8.71e-10)
    parser.add_argument('-tstv', metavar='<tstv>', help='Transition/travsersion ratio.', default='2.0')
    parser.add_argument('-out', metavar='<out>', help='Output file type: fasta or vcf', default='fasta')
    parser.add_argument('-impout', metavar='<impout>', help='Output list of imputed mutations', default=True)
    parser.add_argument('-multi', metavar='<multi>', help='Multiprocessing.', default=True)
    parser.add_argument('-nobackmutchk', dest='nobackmutchk', help='Skip back-mutation check.', action='store_true')
    parser.set_defaults(nobackmutchk=False)
    parser.add_argument('-starttree', metavar='<starttree>', help='Newick or phyloxml starting tree for RAxML')
    parser.add_argument('-verbose', dest='verbose', help='Verbose output.', action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-genoqual', metavar='<genoqual>', help='Genotype Quality threshold for VCF input.', default=30)
    parser.add_argument('-maxthreads', metavar='<maxthreads>', help='Maximum RAxML threads.', default=4)
    parser.add_argument('-adthresh', metavar='<adthresh>', help='Threshold for Allelic Depth.', default=0.66)
    parser.add_argument('-mincoverage', metavar='<mincoverage>', help='Minimum coverage to impute.',
                        default=0)
    parser.add_argument('-passes', metavar='<passes>', help='Number of imputation passes.', default=1)
    parser.add_argument('-mnm', dest='mnm', help='Impute missing when neighbors a missing/non-missing pair.',
                        action='store_true')
    parser.set_defaults(mnm=False)
    parser.add_argument('-rej', dest='rej', help='Create data section of REJECTOR2 infile.', action='store_true')
    parser.set_defaults(rej=False)
    parser.add_argument('-exlocal', dest='exlocal', help='Invoke external programs in local directory.',
                        action='store_true')
    parser.set_defaults(local=False)
    parser.add_argument('-seqonly', dest='seqonly', help='Exit after outputting input seq.', action='store_true')
    parser.set_defaults(seqonly=False)
    parser.add_argument('-boot', metavar='<boot>', help='Number of bootstrap replicates.', default=0)
    parser.add_argument('-threshold', metavar='<threshold>', help='Bootstrap or multiple run acceptance threshold.',
                        default=0.95)
    parser.add_argument('-runs', metavar='<runs>', help='Number of alternate runs.', default=0)
    parser.add_argument('-maxhops', metavar='<maxhops>', help='Number of alternate runs.', default=5)
    parser.add_argument('-ncollect', metavar='<ncollect>', help='rootward, hops, distance, mono',
                        default='rootward')
    parser.add_argument('-maxn', metavar='<maxn>',
                        help='Maximum number of neighbors that will be searched. (For monophyly-based neighbor search.',
                        default=5)
    parser.add_argument('-maxjump', metavar='<maxjump>',
                        help='Maximum increase in branch distance (over previous branch traversed) allowed before terminating neighbor search if -ncollect set to distance or hops.',
                        default=5)
    parser.add_argument('-orphanchk', dest='orphanchk', help='Stop neighbor search if current has no neighbors.',
                        action='store_true')
    parser.set_defaults(orphanchk=False)

    args = parser.parse_args()
    inputfile = args.file
    reffile = args.reffile
    treetype = args.tree
    outtreetype = args.outtree
    alpha = args.alpha
    rmodel = args.rmodel
    maxdepth = int(args.maxdepth)
    maxheight = int(args.maxheight)
    nsize = int(args.nsize)
    msize = int(args.msize)
    mutrate = float(args.mutrate)
    tstv = float(args.tstv)
    outtype = args.out
    impout = args.impout
    multi = args.multi
    nobackmutchk = bool(args.nobackmutchk)
    starttree = args.starttree
    verbose = bool(args.verbose)
    genoqual = int(args.genoqual)
    adthresh = float(args.adthresh)
    maxthreads = int(args.maxthreads)
    mincoverage = int(args.mincoverage)
    passes = int(args.passes)
    mnm = bool(args.mnm)
    rej = bool(args.rej)
    seqonly = bool(args.seqonly)
    exlocal = bool(args.exlocal)
    boot = int(args.boot)
    threshold = float(args.threshold)
    runs = int(args.runs)
    maxhops = int(args.maxhops)
    ncollect = args.ncollect
    maxn = int(args.maxn)
    maxjump = int(args.maxjump)
    orphanchk = bool(args.orphanchk)

    if seqonly:
        verbose = True

    rng = random.SystemRandom()  # Uses /dev/urandom

    sys.setrecursionlimit(2 ** 20)

    batchlist = []

    if os.path.isdir(inputfile):
        for bfile in os.listdir(inputfile):
            if bfile.endswith(".fasta"):
                batchlist.append(inputfile + "/" + bfile)
            elif bfile.endswith(".vcf"):
                batchlist.append(inputfile + "/" + bfile)
            elif bfile.endswith(".vcf.gz"):
                batchlist.append(inputfile + "/" + bfile)
    elif os.path.isfile(inputfile):
        if inputfile.endswith(".fasta"):
            batchlist.append(inputfile)
        elif inputfile.endswith(".vcf"):
            batchlist.append(inputfile)
        elif inputfile.endswith(".vcf.gz"):
            batchlist.append(inputfile)
        else:
            print "Input file must be either .fasta or .vcf or .vcf.gz"
            exit()
    else:
        print "Input file(s) not found. Exiting."
        exit()

    for infile in batchlist:
        randstamp = rng.getrandbits(32)
        print "\n\n****************\nInput file:", infile, "\n****************"
        print "Unique stamp for temp files: ", randstamp

        indata = InData(infile, reffile)

        if seqonly:
            continue

        phytree = PhyloTree(indata, treetype, randstamp, alpha, rmodel, starttree, maxthreads, maxdepth, maxheight,
                            maxn)

        if verbose:
            print "\n****************\nVARIANTS\n****************\n"
            for gi in indata.variants.keys():
                print gi + ": " + str(indata.variants[gi])

            print "\n****************\nSEQUENCE\n****************\n"
            print indata.sequence

            print "\n****************\nTREE\n****************\n"
            Phylo.draw_ascii(phytree.tree)
        phytree.output_tree(outtreetype)

        print "\n****************\nIMPUTATION\n****************\n"
        impute = Imputation(indata, phytree, mutrate, tstv)
        impute.impute()
        impute.output_imputed(impout)
    exit()

else:
    print("IMPUTOR is being imported into another module. Not yet implemented.")

