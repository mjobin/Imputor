#!/usr/bin/python

""" imprun - script for running Imputor many times and comparing output.

    Author - Matthew Jobin, Department of Anthropology, Santa Clara University and Department of Anthropology,
        University of California, Santa Cruz.
    """

import argparse
import os
import random
import shutil
import glob
import subprocess


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def bash_command(cmd):
    outfile.write(cmd)
    outfile.write("\n\n")
    p = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # p.wait()
    stdout = []
    while True:
        line = p.stdout.readline()
        if verbose:
            print line,
        if line.startswith('Ratio Pairwise Distance for all infiles:'):
            cols = line.split(":")
            ratio = float(cols[1])
            return ratio
        stdout.append(line)
        # print line,
        if line == '' and p.poll() != None:
            break
            print "ERROR Ratio not found!"
            exit()
    return ''.join(stdout)



if __name__ == "__main__":

    print "\n****************\nIMPRUN\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                    "n"

                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-reps', metavar='<reps>',
                        help='reps.',
                        default=10)
    parser.add_argument('-nseqs', metavar='<nseqs>',
                        help='nseqs.',
                        default=50)
    parser.add_argument('-miss', metavar='<miss>',
                        help='miss',
                        default="0.01")
    parser.add_argument('-maxheight', metavar='<maxheight>',
                        help='Height of search toward root for collecting neighbors.',
                        default="3")
    parser.add_argument('-maxdepth', metavar='<maxdepth>',
                        help='Depth of search down descendent nodes for collecting neighbors.', default=3)
    parser.add_argument('-nsize', metavar='<nsize>',
                        help='Number of neighbors that must be identical in order to impute.',
                        default="3")
    parser.add_argument('-msize', metavar='<msize>',
                        help='Number of neighbors that must be identical in order to impute for missing data.',
                        default="2")
    parser.add_argument('-ncollect', metavar='<ncollect>', help='rootward, hops, distance, mono',
                        default='rootward')
    parser.add_argument('-maxhops', metavar='<maxhops>', help='Number of alternate runs.', default="5")


    args = parser.parse_args()
    wd = args.wd
    verbose = bool(args.verbose)
    reps = int(args.reps)
    nseqs = int(args.nseqs)
    miss = args.miss
    maxdepth = args.maxdepth
    maxheight = args.maxheight
    nsize = args.nsize
    msize = args.msize
    ncollect = args.ncollect
    maxhops = args.maxhops

    rng = random.SystemRandom()  # Uses /dev/urandom

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd
    print "nseqs: " + str(nseqs)

    outfile = open("i_cmds", 'w')

    missfrags = miss.split(".")
    missdec = missfrags[1]

    raws = []
    imps = []
    shps = []
    for i in range(reps):
        print str(i+1) + " of " + str(reps)

        # randstamp = rng.getrandbits(32)
        randstamp = "test"
        randbase = "irf" + str(randstamp)
        shrandbase = "shrf" +  str(randstamp)

        for fl in glob.glob(randbase + "*"):
            os.remove(fl)
        for fl in glob.glob(shrandbase + "*"):
            os.remove(fl)



        invcf = randbase + ".vcf"
        infasta = randbase + ".fasta"
        inref = randbase + "-indata-ref.fasta"



        bash_command("sfs_code 1 1 -t 0.01 -L 1 10000 -x 1 -n " + str(nseqs) + " --VCF -P 1 --propFemale 0 -o " + invcf)

        bash_command("python ./imputor.py -verbose -seqonly -file " + invcf)



        curfile = randbase + "-miss" + missdec + "-mu0-0.fasta"

        bash_command("fastablasta.py -infile " + infasta + " -miss " + miss + " -reps 1")
        raw = bash_command("fastadiff.py -silent -afile " + infasta + " -bfile " + curfile)
        bash_command("imputor.py -exlocal -maxthreads 12 -ncollect " + ncollect + " -maxheight " + maxheight + " -maxdepth " + maxdepth + " -passes 1 -msize " + msize + " -nsize " + nsize + " -maxhops " + maxhops + " -file " + curfile)

        raws.append(raw)

        seqoutfile = randbase + "-miss" + missdec + "-mu0-0-seqout.fasta"
        imp = bash_command("fastadiff.py -silent -seqout -afile " + infasta + " -bfile " + seqoutfile)
        imps.append(imp)


        #shapeit

        shutil.copy(randbase + "-miss" + missdec + "-mu0-0.fasta", shrandbase + "-miss" + missdec + "-mu0-0.fasta")
        shcurbase = shrandbase + "-miss" + missdec + "-mu0-0"
        shcurfile = shcurbase + ".fasta"
        bash_command("imputor.py -seqonly -verbose -out vcf -file " + shcurfile + " -reffile " + inref)
        bash_command("plink --vcf " + shcurbase + ".vcf --recode --out " + shcurbase)
        shcurmale = shcurbase + "-male"
        bash_command("awk '{OFS=\"\t\";if($5==\"0\") {$5=\"1\";} else if($5==\"F\") {$5=\"1\";} else {$5=\"1\";}print;}' " + shcurbase + ".ped > " + shcurmale + ".ped")
        shutil.copy(shcurbase + ".map", shcurmale + ".map")
        bash_command("shapeit --chrX --input-ped " + shcurmale + " -O " + shcurmale)
        bash_command("shapeit -convert --input-haps " + shcurmale + " --output-vcf " + shcurmale + ".phased.vcf")

        for fl in glob.glob("shapeit*"):
            os.remove(fl)

        bash_command("python ./vcfstrip.py -infile " + shcurmale + ".phased.vcf")
        bash_command("python ./imputor.py -seqonly -verbose -out fasta -file " + shcurmale + ".phased-s.vcf")
        rfile = shcurmale + ".fasta"

        shp = bash_command("python ./fastadiff.py -silent -seqout -afile " + infasta + " -bfile " + rfile)

        if not is_number(shp):
            print "shp not number " + shrandbase
            print shp
            exit()
        shps.append(shp)

        for fl in glob.glob(randbase + "*"):
            os.remove(fl)
        for fl in glob.glob(shrandbase + "*"):
            os.remove(fl)

    mraw = sum(raws) / float(len(raws))
    mimp = sum(imps) / float(len(imps))
    mshp = sum(shps) / float(len(shps))


    print "RAW Mean: " + str(mraw)
    print "IMPUTED Mean: " + str(mimp)
    print "SHAPEIT Mean: " + str(mshp)


    outfile.close()
