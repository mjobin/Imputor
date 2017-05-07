# Imputor

""" Imputor - software for imputing missing mutations in sequencing data by the use
    of phylogenetic trees.
    
    Author - Matthew Jobin
    """

import argparse
import glob
import multiprocessing
import os
import random
import sys
from argparse import RawTextHelpFormatter

import progressbar
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.Applications import PhymlCommandline
from Bio.Phylo.Applications import RaxmlCommandline
from Bio.Phylo.TreeConstruction import *
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class InData(object):
    """Input data from FASTA or VCF file converted internally to sequence data.
        Prepares a dictionary of variants and a list of sequences for processing
        and imputation.
    """

    def __init__(self):
        self.fullsequence = []  # Raw sequence input including non-segregating sites
        self.fullvariantset = set()
        self.fullvariants = {}
        self.sequence = []  # Sequence data
        self.variantset = set()
        self.reflist = []
        self.refset = set()
        self.variants = {}  # Dictionary of each sample and its variation from the reference sequence
        self.filebase = None
        self.orig_vcf_pos = []  # Original listed positions of variants
        self.geninfo = {}  # Dict of list of dicts
        self.gqdict = {}
        self.addict = {}
        self.maxseqlength = 0
        self.revrate = -1

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

        return

    def load_input_data(self, inputfile=None):
        """ Load input data from a file.
            
            Keyword arguments:
            inputfile -- Input data. Accepted formats: FASTA, VCF.
        """

        filebase, fileext = os.path.splitext(inputfile)
        self.filebase = filebase
        if inputfile[-5:] == 'fasta':
            self.sequence = AlignIO.read(inputfile, 'fasta')
            self.variants_from_sequence()
            self.prune_non_seg()
        elif inputfile[-3:] == 'vcf':
            file_data = open(inputfile, 'r')
            raw_data = []
            for file_line in file_data:
                raw_data.append(file_line.rstrip())

            # Split multi-allelic sites
            expanded_file_data = self.vcf_expand_multi_allele(raw_data)

            # Generate sequence from only those areas with any polymorphism
            self.seq_from_variants_excl(raw_data)
        else:
            print "Input file must be either .fasta or .vcf"
            exit()

        if verbose:
            outseqfile = filebase + "-indata.fasta"
            outfile = open(outseqfile, 'w')
            outseq = {}
            for seq in self.sequence:
                outseq[seq.id] = str(seq.seq)
            for x in sorted(outseq.keys()):
                outfile.write(">")
                outfile.write(str(x))
                outfile.write("\n")
                outfile.write(outseq[x])
                outfile.write("\n")
            outfile.close()
            outreffilename = filebase + "-indata-ref.fasta"
            outreffile = open(outreffilename, 'w')
            outreffile.write(">REF")
            outreffile.write("\n")
            outreffile.write("".join(self.reflist))
            outreffile.write("\n")
            outreffile.close()
        if rej:
            self.rej_infile(inputfile)

        print "Finished input."
        return

    def vcf_expand_multi_allele(self, in_data=None):
        """ Processes multiple ALT alleles in a VCF file. Returns expanded data in list form.
            
            Keyword arguments:
            in_data-- VCF file data.
        """
        expanded_file_data = []
        print "\nExpanding multi-allele entries..."
        bar = progressbar.ProgressBar(redirect_stdout=True)
        for i in bar(range(len(in_data))):
            file_line = in_data[i]
            cols = file_line.split('\t')

            # If the second character is a (meta-info line) or a blank line, ignore
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n'):
                continue
            if "," in cols[self.vcf_alt]:
                multi_allele = cols[self.vcf_alt].split(",")
                for allele in multi_allele:
                    allele_cols = cols
                    # Replace the ALT column with just this allele
                    allele_cols[self.vcf_alt] = allele
                    # Place line with just this allele in the expanded data
                    expanded_file_data.append("\t".join(allele_cols))
            else:
                expanded_file_data.append(file_line)
        return expanded_file_data

    def vcf_snp_prune(self, in_data=None):
        """ Returns a VCF file including only lines containing SNPs.
            
            Keyword arguments:
            in_data-- VCF file data.
        """
        snps_data = []
        print "\nPruning non-SNP entries..."
        bar = progressbar.ProgressBar()
        for i in bar(range(len(in_data))):
            file_line = in_data[i]

            cols = file_line.split('\t')
            # If the second character is a (meta-info line) or a blank line, ignore
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n'):
                continue
            # if len(cols[self.vcf_ref]) > 1 or len(cols[self.vcf_alt]) > 1:  # if not a snp
            #     continue
            # elif cols[self.vcf_ref] == '-' or cols[self.vcf_alt] == '-': # if not a snp
            #     continue
            # elif cols[self.vcf_ref] == '.' or cols[self.vcf_alt] == '.': # if not a snp
            #     continue
            snps_data.append(file_line)
        return snps_data

    def variants_from_sequence(self):
        """ Returns a list of variants from a reference sequence.
            Intended for use with FASTA input, but will work with any AlignIO object or
            list of sequence data.
        """
        print "\nGenerating variants from sequence..."
        firstseq = self.sequence[0]
        bar = progressbar.ProgressBar()
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
        """
        self.fullsequence = self.sequence  # First back up the original sequence
        self.fullvariantset = self.variantset
        self.fullvariants = self.variants
        self.sequence = MultipleSeqAlignment([])  # Blank the sequence to be worked on

        print "Pruning non-segregating sites..."
        locs = []
        for curvar in self.variantset:
            locs.append(curvar)
        locs.sort()

        stripped = {}
        seqnames = []
        for seq in self.fullsequence:
            stripped[seq.name] = []
            seqnames.append(seq.name)

        for loc in locs:
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

    def seq_from_variants_excl(self, raw_data=None):
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
        bar = progressbar.ProgressBar()
        for i in bar(range(len(snps_data))):
            file_line = snps_data[i]
            # For file_line in snps_data:
            cols = file_line.split('\t')
            if int(cols[self.vcf_pos]) > self.maxseqlength:
                self.maxseqlength = int(cols[self.vcf_pos])
            self.orig_vcf_pos.append(cols[self.vcf_pos])
            self.refset.add(str(var_count) + cols[self.vcf_ref])
            if (file_line[:1] == "#") or (cols[self.vcf_chrom] == '\n' or cols[self.vcf_info + 1][:2] != "GT"):
                continue
            formatcols = cols[self.vcf_info + 1].split(":")

            indiv_genotypes = cols[self.vcf_info + 2:]  # Assumes all rows same length, as per VCF standard
            for position, indiv_genotype in enumerate(
                    indiv_genotypes):  # Iterates through that row of genotypes for this site
                genotypecols = indiv_genotype.split(":")
                assigned_alleles = genotypecols[0].split(  # VCF standard GT always first column
                    "[/|]+")  # Split genotype entry on either character phased or unphased
                changed_genotype_names = []
                for allele_pos, assigned_allele in enumerate(assigned_alleles):  # Iterates through the alleles
                    changed_genotype_name = genotype_names[position]
                    if len(assigned_alleles) > 1:  # Only append to genotype name if not haploid
                        changed_genotype_name += str(allele_pos)
                    changed_genotype_names.append(changed_genotype_name)
                for changed_genotype_name in changed_genotype_names:
                    if changed_genotype_name not in genotype_sequence:
                        genotype_sequence[changed_genotype_name] = []

                alt_alleles = cols[self.vcf_alt].split(",")  # List of ALT alleles for this row
                for allele_pos, assigned_allele in enumerate(assigned_alleles):  # Iterates through the alleles
                    if assigned_allele == "0":  # Assigned_allele will be 0 for REF and >0 for any ALT
                        genotype_sequence[changed_genotype_names[allele_pos]].append(cols[self.vcf_ref])
                    elif assigned_allele == ".":  # VCF format code for missing allele
                        genotype_sequence[changed_genotype_names[allele_pos]].append("N")
                    else:
                        genotype_sequence[changed_genotype_names[allele_pos]].append(
                            alt_alleles[int(assigned_allele) - 1])
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

    def rej_infile(self, inputfile=None):
        """ Outout a REJECTOR2 input file

            Keyword arguments:
            inputfile -- Input data. Accepted formats: FASTA, VCF.
        """

        filebase, fileext = os.path.splitext(inputfile)
        rejfilename = filebase + "-rej.txt"

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
    """A phylogenetic tree either input from phyloxml format or constructed from sequence
    """

    def __init__(self):
        self.impname = None
        self.starttreename = None
        self.treetype = treetype
        self.tree = None  # Phylogenetic tree to be loaded or constructed from data. Newick format.
        self.treeparents = {}
        self.starttree = None  # Phylogenetic tree used as starting tree in RAxML

    def input_tree(self, treetype=None, alpha=None, rmodel=None, starttreename=None, timestamp=None, maxthreads=None):
        """ Takes input tree file or sequence data.

            Keyword arguments:
            treetype -- type of tree to be input or constructed
        """

        self.impname = "imp" + str(timestamp)
        self.starttreename = starttreename

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
        elif self.treetype == 'RAxML':
            self.raxml_tree(rmodel)
        else:
            self.phyml_tree(alpha)

        self.treeparents = self.all_parents(self.tree)

    def output_tree(self, inputfile, outtreetype):
        """Outputs tree to file

            Keyword arguments:
            inputfile -- name if input file
            outtreetype -- format of tree file for output
        """
        filebase, fileext = os.path.splitext(inputfile)
        if outtreetype == 'phyloxml':
            outfile = filebase + "-outtree.xml"
            Phylo.write(self.tree, outfile, "phyloxml")
        elif outtreetype == 'nexus':
            outfile = filebase + "-outtree.nexus"
            Phylo.write(self.tree, outfile, "nexus")
        else:  # Default newick
            outfile = filebase + "-outtree.newick"
            Phylo.write(self.tree, outfile, "newick")

    @staticmethod
    def all_parents(tree):
        parents = {}
        for clade in tree.find_clades():
            for child in clade:
                parents[child] = clade
        return parents

    def collect_kids(self, clade, kids, depth, maxdepth):
        """ Recursive function for collecting all children to specified depth

            Keyword arguments:
            clade -- current node of tree
            kids --list of all collected children
            depth -- current depth of search
            maxdepth -- maximum depth of search
        """
        #
        if depth < maxdepth:
            for child in clade:
                if child.name:
                    kids.append(child)
                self.collect_kids(child, kids, depth + 1, maxdepth)
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

    @staticmethod
    def get_siblings(clade):
        siblings = []
        if clade in phytree.treeparents:  # Will not go past the root clade
            for kid in phytree.treeparents[clade]:
                if kid == clade:
                    pass
                else:
                    siblings.append(kid)
        return siblings

    def collect_kids_rootward(self, term, parents, height, maxheight, maxdepth, maxneighbors):
        neighbors = set()
        curnode = term
        # Collect closest neighbors on tree
        while len(neighbors) < maxneighbors:
            if curnode not in parents:  # will not go past the root
                break
            if height >= maxheight:
                break
            curparent = parents[curnode]

            empty = []
            allkids = self.collect_kids(curparent, empty, 0, maxdepth)
            for kid in allkids:
                if kid is not term:
                    neighbors.add(kid)
                    if len(neighbors) >= maxneighbors:
                        break
            curnode = curparent
            height += 1
        return neighbors

    def parsimony_tree(self):
        """ Constructs a tree via maximum parsimony using Biopython's ParsimonyTreeConstructor.

        """
        print "Generating maximum parsimony tree.."
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
        cpus = multiprocessing.cpu_count()
        if cpus > maxthreads:
            cpus = maxthreads

        # Erase RaXML intermediate files from previous runs
        raxml_glob = glob.glob('RAxML_*')
        for delfile in raxml_glob:
            os.remove(delfile)

        # Output sequence to a temp FASTA file
        tempfastafile = indata.filebase + self.impname + "_fastatmp.fasta"
        reducedtempfastafile = indata.filebase + self.impname + "_fastatmp.fasta.reduced"
        AlignIO.write(indata.sequence, tempfastafile, "fasta")

        raxml_args = {"sequences": tempfastafile, "model": rmodel, "name": self.impname,
                      "parsimony_seed": rng.randint(0, sys.maxint), "threads": cpus}

        raxmlstarttreename = "RAxML_" + self.impname + "_starttree.newick"
        if self.starttree:
            Phylo.write(self.starttree, raxmlstarttreename, "newick")
            raxml_args["starting_tree"] = raxmlstarttreename

        raxml_cline = ""
        if exlocal:
            raxml_cline = RaxmlCommandline(cmd='./raxmlHPC', **raxml_args)
        else:
            raxml_cline = RaxmlCommandline(**raxml_args)

        print "Invoking RAxML with ", raxml_cline

        out_log, err_log = raxml_cline()
        if verbose:
            print err_log
            print out_log
        raxmlbesttreename = "RAxML_bestTree." + self.impname
        self.tree = Phylo.read(raxmlbesttreename, "newick")

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

    def phyml_tree(self, alpha=None):
        """ Constructs a tree via maximum likelihood by invoking external software PhyML.
            See docs for PhyML installation and setup.

        """
        print "Invoking PhyML..."
        # Output sequence to a temp FASTA file
        tempfastafile = indata.filebase + "_" + self.impname + "_fastatmp.fasta"
        AlignIO.write(indata.sequence, tempfastafile, "fasta")
        tempphyfile = indata.filebase + "_" + self.impname + "_phytmp.phy"
        AlignIO.convert(tempfastafile, "fasta", tempphyfile, "phylip-relaxed")

        phyml_args = {"input": tempphyfile, "alpha": "e"}
        phystarttreename = "PhyML_imp", self.impname, "starttree.newick"
        if self.starttree:
            Phylo.write(self.starttree, phystarttreename, "newick")
            phyml_args["input_tree"] = phystarttreename

        cmdline = ""
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
            phyml_globname = indata.filebase + "_" + self.impname + "*"
            phyml_glob = glob.glob(phyml_globname)
            for delfile in phyml_glob:
                os.remove(delfile)


class Imputation(object):
    """Imputation of mutations given input data and a phylogenetic tree.
    """

    def __init__(self, indata, phytree, mutrate, multi):
        self.cpucount = multiprocessing.cpu_count()
        self.workseq = {}
        self.imputedseq = MultipleSeqAlignment([])
        self.indata = indata
        self.phytree = phytree
        self.mu = mutrate
        self.imputelist = []
        self.reversionlist = []
        self.multi = multi
        self.acgt = {'A', 'C', 'G', 'T'}
        self.missing = {'.', '-', 'N'}
        self.indivimputes = {}
        self.neighbors = {}

        for seq in indata.sequence:
            self.workseq[seq.name] = list(str(seq.seq))  # Begin with output sequence matching input
            self.indivimputes[seq.name] = []

    def impute(self):
        """ Sets up imputation function for all terminal nodes.

            Keyword arguments:
            imputetype -- Type of search used to locate possible sites to impute.
            depth -- Depth of search up and down tree to find neighbours.
            neighbors -- Minimum number of identical non-missing neighbors in order to impute.

        """
        terms = self.phytree.tree.get_terminals()  # Get all internal nodes on tree. These are the ones with samples.
        for i in range(passes):
            print "\nPass", i
            random.shuffle(terms)
            bar = progressbar.ProgressBar()
            for p in bar(range(len(terms))):
                # for term in terms:
                term = terms[p]
                neighbors = phytree.collect_kids_rootward(term, self.phytree.treeparents, 0, maxheight, maxdepth,
                                                          maxneighbors)
                self.impute_threshold(term, self.phytree.treeparents, neighbors, str(i + 1))
            for newimpute in self.imputelist:
                if newimpute[5] == "T":
                    self.indivimputes[newimpute[0]].append(newimpute[1])
                    self.workseq[newimpute[0]][newimpute[1]] = newimpute[3]
        self.process_imputed()

    def impute_threshold(self, term, parents, neighbors, thispass):
        """ Imputation on best tree.

            Keyword arguments:
            term -- Terminal node on tree to be examined.
            parents - List of parent nodes of each node in tree.
            neighbors -- Minimum number of identical neighbors to impute.

        """
        for curvar in self.indata.variantset:
            newimpute = self.detect_by_parsimony(term, curvar, parents, neighbors)
            newimpute.append(thispass)
            if verbose:
                self.imputelist.append(newimpute)
            else:
                if newimpute[5] == "T":
                    self.imputelist.append(newimpute)

    def backmutscan(self):
        terms = self.phytree.tree.get_terminals()
        random.shuffle(terms)
        backmuts = 0
        allmuts = 0
        for term in terms:
            termname = str(term)
            neighbors = phytree.collect_kids_rootward(term, self.phytree.treeparents, 0, maxheight, maxdepth,
                                                      maxneighbors)

            self.neighbors[term] = neighbors
            if len(neighbors) < maxneighbors:
                continue
            for curvar in self.indata.variantset:
                orig = self.workseq[termname][curvar]
                if orig in self.missing:
                    continue
                origvar = str(indata.orig_vcf_pos[curvar])

                nearest = set()
                for neighbor in neighbors:
                    nearest.add(self.workseq[str(neighbor)][curvar])
                if len(nearest) == 1:
                    only = nearest.pop()
                    if only in self.missing or orig == only:
                        continue
                    if self.backmutchk(term, self.phytree.treeparents, neighbors, curvar, orig):
                        self.reversionlist.append([termname, origvar, orig, "T"])
                        backmuts += 1
                    else:
                        self.reversionlist.append([termname, origvar, orig, "F"])
                    allmuts += 1
        print "Reversions: ", backmuts, " in ", allmuts

    def process_imputed(self):
        print "\nProcessing imputed sequences..."
        locs = []
        for curvar in indata.fullvariantset:
            locs.append(curvar)
        locs.sort()

        bar = progressbar.ProgressBar()
        for p in bar(range(len(indata.fullsequence))):
            fullseq = indata.fullsequence[p]
            tmpseq = list(fullseq)
            segseq = self.workseq[fullseq.id]

            if len(segseq) == len(locs):
                for site, loc in itertools.izip(segseq, locs):
                    tmpseq[loc] = site

            seqrec = SeqRecord(Seq("".join(tmpseq)), id=fullseq.id, name=fullseq.id)
            self.imputedseq.append(seqrec)

        self.imputedseq.sort()

    @staticmethod
    def transversionchk(a, b, tstv):
        # print "from " + str(a) + " to " + str(b)
        if a == "A":
            if b == "G":
                return tstv
        if a == "G":
            if b == "A":
                return tstv
        if a == "C":
            if b == "T":
                return tstv
        if a == "T":
            if b == "C":
                return tstv
        if a in ("A", "C", "G", "T") and b in ("A", "C", "G", "T"):
            return 1 / tstv
        return -1

    def backmutchk(self, term, parents, allneighbours, curvar, origseq):
        backmut = False
        allneighbours.add(term)  # don't check self
        curparent = parents[term]
        while curparent in phytree.treeparents:  # Search for allele farther away
            if backmut:
                return True
            nextparent = parents[curparent]
            empty = []
            allkids = phytree.collect_all_kids(curparent, empty)
            for kid in allkids:
                kidseq = self.workseq[str(kid)][curvar]
                if kidseq in self.missing:
                    continue
                if kid not in allneighbours and kidseq == origseq:
                    backmut = True
                allneighbours.add(kid)
            curparent = nextparent
        return False

    def detect_by_parsimony(self, term, curvar, parents, neighbors):

        termname = str(term)
        finame = termname + "-" + str(curvar)

        nearest = set()
        orig = self.workseq[termname][curvar]

        for neighbor in neighbors:
            nearest.add(self.workseq[str(neighbor)][curvar])

        if len(nearest) > 1:  # Cannot allow non-matching ever to impute sequence, thus this is first check
            if mnm and orig in self.missing and len(nearest) == 2 and "N" in nearest:
                nearest.remove("N")  # Strip the one N, allowing imputation
            else:
                return [termname, curvar, orig, ",".join(nearest), "Neighbors Non-matching", "F"]
        if len(nearest) < 1:
            return [termname, curvar, orig, ".", "No neighbors", "F"]
        only = nearest.pop()
        if orig in self.missing:
            if len(neighbors) < 2:
                return [termname, curvar, orig, only, "Not enough neighbors", "F"]
            elif only in self.missing:
                return [termname, curvar, orig, only, "Neighbors All Missing", "F"]
            else:
                return [termname, curvar, orig, only, "Imputed Missing", "T"]
        else:
            if len(neighbors) < maxneighbors:
                return [termname, curvar, orig, only, "Not enough neighbors", "F"]
            if only in self.missing:
                return [termname, curvar, orig, only, "Neighbors All Missing", "F"]
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
                        return [termname, curvar, orig, only, "Imputed by Min. Cov.", "T"]
                    if otherad > 0.0:
                        if thisad / otherad < adthresh:
                            return [termname, curvar, orig, only, "Imputed by AD", "T"]
                if finame in indata.gqdict:
                    if int(indata.gqdict.get(finame)) < genoqual:
                        return [termname, curvar, orig, only, "Imputed by GQ", "T"]
                if nobackmutchk:
                    return [termname, curvar, orig, only, "Imputed Non-missing", "T"]
                else:
                    if self.backmutchk(term, parents, neighbors, curvar, orig):
                        return [termname, curvar, orig, only, "Imputed Non-missing", "T"]
                    else:
                        return [termname, curvar, orig, only, "No Back Mutation", "F"]
            else:
                return [termname, curvar, orig, only, "Matches Neighbors", "F"]

    def output_imputed(self, inputfile, out, impout):

        for imputed in self.imputelist:
            if indata.orig_vcf_pos:
                imputed[1] = str(indata.orig_vcf_pos[imputed[1]])
            else:
                imputed[1] = str(imputed[1])

        if verbose:
            if len(self.imputelist) > 0:
                print "Imputed Mutations"
                print "SUBJECTID | VAR | FROM | TO | TYPE | IMPUTED | PASS"
                for imputed in self.imputelist:
                    print " | ".join(imputed)
                print "\n"
            print impute.imputedseq

        filebase, fileext = os.path.splitext(inputfile)
        if impout:
            impoutfilename = filebase + "-impout.txt"
            impoutfile = open(impoutfilename, 'w')
            impoutfile.write("SUBJECTID\t VAR\t FROM\t TO\t TYPE\t IMPUTED\tPASS\n")
            for imputed in self.imputelist:
                impoutfile.write("\t".join(imputed))
                impoutfile.write("\n")
            impoutfile.close()
        indivoutfilename = filebase + "-indivout.txt"
        indivoutfile = open(indivoutfilename, 'w')
        indivoutfile.write("SUBJECTID\t NUM\n")

        revfilename = filebase + "-rev.txt"
        revoutfile = open(revfilename, 'w')
        revoutfile.write("SUBJECTID\t VAR\t SITE\t REV?\n")
        for rev in self.reversionlist:
            revoutfile.write("\t".join(rev))
            revoutfile.write("\n")
        revoutfile.close()

        for indiv in self.indivimputes.keys():
            indivoutfile.write(indiv)
            indivoutfile.write("\t")
            indivoutfile.write(str(len(self.indivimputes[indiv])))
            indivoutfile.write("\t")
            for indivar in self.indivimputes[indiv]:
                indivoutfile.write(str(indivar))
                indivoutfile.write("\t")
            indivoutfile.write("\n")
        indivoutfile.close()
        if out == "vcf":
            outseqfile = filebase + "-seqout.vcf"
            outfile = open(outseqfile, 'w')
            outfile.write("##fileformat=VCFv4.1")
            outfile.write("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	")
        else:  # default to fasta
            outseqfile = filebase + "-seqout.fasta"
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


if __name__ == "__main__":

    print "\n****************\nIMPUTOR\n****************\n"

    parser = argparse.ArgumentParser(description="IMPUTOR is a program for the phylogeny-aware imputation of: \n\n\t" \
                                                 "mutations.\n\t" \
                                                 "- ", formatter_class=RawTextHelpFormatter)

    parser.add_argument('-file', metavar='<file>', help='input file: FASTA or VCF', required=True)
    parser.add_argument('-tree', metavar='<tree>', help='tree type; <treefilename.xml>, pars, RAxML, PhyML',
                        default='RAxML')
    parser.add_argument('-outtree', metavar='<outtree>', help='Output format for tree', default='newick')
    parser.add_argument('-alpha', metavar='<alpha>', help='Value of gamma shape parameter.', default='e')
    parser.add_argument('-rmodel', metavar='<rmodel>', help='Model type for RaXML.', default='GTRCAT')
    parser.add_argument('-maxheight', metavar='<maxheight>',
                        help='Height of search toward root for collecting neighbors.',
                        default=2)
    parser.add_argument('-maxdepth', metavar='<maxdepth>',
                        help='Depth of search down descendent nodes for collecting neighbors.', default=2)
    parser.add_argument('-maxneighbors', metavar='<maxneighbors>',
                        help='Number of neighbors that must be identical in order to impute.',
                        default=3)
    parser.add_argument('-hapobj', metavar='<hapobj>', help='Haplogroups file.')
    parser.add_argument('-wtobj', metavar='<wtobj>', help='Weights file.')
    parser.add_argument('-polyobj', metavar='<polyobj>', help='Polymorphisms file.')
    parser.add_argument('-mutrate', metavar='<mutrate>', help='Mutation rate.', default='8.71e-10')
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
    parser.add_argument('-maxthreads', metavar='<maxthreads>', help='Maximum RAxML Pthreads.', default=4)
    parser.add_argument('-adthresh', metavar='<adthresh>', help='Threshold for Allelic Depth.', default=0.66)
    parser.add_argument('-mincoverage', metavar='<mincoverage>', help='Minimum coverage to ignore back mutation check.',
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

    # BATCH and RANDOM switches
    parser.add_argument('-seqlength', metavar='<seqlength>',
                        help='Length of randomly generated sequence.', default=1000)
    parser.add_argument('-seqcount', metavar='<seqcount>',
                        help='Number of randomly generated sequences to create.', default=1000)
    parser.add_argument('-seqonly', dest='seqonly', help='Exit after outputting input seq.', action='store_true')
    parser.set_defaults(seqonly=False)
    parser.add_argument('-indelsize', metavar='<indelsize>', help='Max indel size', default=0)
    parser.add_argument('-indelrate', metavar='<indelrate>', help='Indel rate.', default=0.0)
    parser.add_argument('-randbatch', metavar='<randbatch>',
                        help='Number of times to run with randomly generated input.', default=100)

    args = parser.parse_args()
    inputfile = args.file
    treetype = args.tree
    outtreetype = args.outtree
    alpha = args.alpha
    rmodel = args.rmodel
    hapobj = args.hapobj
    wtobj = args.wtobj
    polyobj = args.polyobj
    maxdepth = int(args.maxdepth)
    maxheight = int(args.maxheight)
    maxneighbors = int(args.maxneighbors)
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
    idsize = int(args.indelsize)
    exlocal = bool(args.exlocal)

    rng = random.SystemRandom()  # Uses /dev/urandom

    batchlist = []

    if os.path.isdir(inputfile):
        for bfile in os.listdir(inputfile):
            if bfile.endswith(".fasta"):
                batchlist.append(inputfile + "/" + bfile)
            elif bfile.endswith(".vcf"):
                batchlist.append(inputfile + "/" + bfile)
    elif os.path.isfile(inputfile):
        if inputfile.endswith(".fasta"):
            batchlist.append(inputfile)
        elif inputfile.endswith(".vcf"):
            batchlist.append(inputfile)
        else:
            print "Input file must be either .fasta or .vcf"
            exit()
    else:
        print "Input file(s) not found. Exiting."
        exit()

    for infile in batchlist:
        randstamp = rng.getrandbits(32)
        print "\n\n****************\nInput file:", infile, "\n****************"
        print "Unique stamp for temp files: ", randstamp

        indata = InData()
        indata.load_input_data(inputfile=infile)

        if seqonly:
            continue

        phytree = PhyloTree()
        phytree.input_tree(treetype=treetype, alpha=alpha, rmodel=rmodel, starttreename=starttree,
                           timestamp=randstamp,
                           maxthreads=maxthreads)

        if verbose:
            print "\n****************\nVARIANTS\n****************\n"
            for i in indata.variants.keys():
                print i + ": " + str(indata.variants[i])

            print "\n****************\nSEQUENCE\n****************\n"
            print indata.sequence

            print "\n****************\nTREE\n****************\n"
            Phylo.draw_ascii(phytree.tree)
        phytree.output_tree(inputfile, outtreetype)

        print "\n****************\nIMPUTATION\n****************\n"
        impute = Imputation(indata, phytree, mutrate, multi)
        impute.impute()
        impute.output_imputed(inputfile, outtype, impout)

else:
    print("IMPUTOR is being imported into another module. Not yet implemented.")
