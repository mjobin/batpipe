#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated April 2018    #####
#####       MJJ                 #####
#####################################

# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# Wrapper for pmdtools

import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import progressbar
import datetime
import gzip
import shutil
import fnmatch
import subprocess
from subprocess import Popen, PIPE

def bash_command(cmd):
    cmdfile.write(cmd)
    cmdfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    stdout, stderr = subp.communicate()
    if verbose:
        print stdout
    logfile.write(stdout)
    if verbose:
        print stderr
    logfile.write(stderr)
    return stdout


if __name__ == "__main__":

    print "\n****************\nPMDTools\n****************\n"

    parser = argparse.ArgumentParser(description="# This script is a wrapper for pmdtools:\n"
                                                "MAKE SURE you know the names of your files!\n"
                                                "including the q number, or CRASH AND DIEEEE!\n"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.', required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-pmdref', metavar='<pmdref>', help='',
                        default="hg19")
    parser.add_argument('-q', metavar='<q>', help='BWA min quality. 20 provides a fairly low cutoff',
                        default="20")
    parser.add_argument('-pmd_threshold', metavar='<pmd_threshold>', help='PMDtools threshold',
                        default="3")
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)


    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    pmdref = args.pmdref
    q = args.q
    pmd_threshold = args.pmd_threshold
    overwrite = bool(args.overwrite)
    verbose = bool(args.verbose)

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", wd

    today = datetime.date.today()
    logfilename = wd + "/out.pmd." + str(today) + ".log"
    print "Logging to: ", logfilename






    logfile = open(logfilename, 'w')


    cmdfile = open("4_cmds", 'w')

    pmd_out = wd + "/PMDtools"
    if overwrite:
        if os.path.exists(pmd_out):
            shutil.rmtree(pmd_out)
        os.mkdir(pmd_out)
    else:
        if os.path.exists(pmd_out):
            print "ERROR: Directory " + pmd_out + " exists and overwrite not set to true. Exiting."
            exit(1)
        else:
            os.mkdir(pmd_out)

    bcin = open(bcfile, 'r')
    bc = []
    for bcinline in bcin:
        bc.append(bcinline)

    bclength = len(bc)
    print "Number of entries: ", bclength

    if os.path.isfile("./PMD_temp.pdf"):
        os.remove("./PMD_temp.pdf")

    print "\nRunning PMDtools..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        sample = bccols[1]
        output =wd + "/" + sample
        bwa_output1 = output + "/BWA_" + pmdref
        bo_s = bwa_output1 + "/" + sample
        po_s = pmd_out + "/" + sample

        logfile.write("PMDtools for " + sample)
        ### assign PMD score to each read and outputs all above X pmd_threshold to a new bam file ###
        bash_command("samtools view -h " + bo_s + "_allreads.cf." + pmdref + ".q" + q + ".s.bam | pmdtools.py --threshold " + pmd_threshold + " --header --stats | samtools view -Sb - > " + bo_s + "_allreads.cf." + pmdref + ".q" + q + ".pmds" + pmd_threshold + "filter.bam")

        ### plot damage for unfiltered bam for comparison with pmd filtered bam ##
        bash_command("samtools view " + bo_s + "_allreads.cf." + pmdref + ".q" + q + ".s.bam | pmdtools.py --deamination --range 30 --CpG > PMD_temp.txt")
        bash_command("R CMD BATCH /data/scripts/plotPMD.R")
        shutil.copy("./PMD_plot.pdf", po_s + "_allreads.cf." + pmdref + ".q" + q + ".NOTpmdfiltered.plot.pdf")

        ### plot damage for pmd filtered bams - accounting for UDG treatment (and rename plot) ##
        bash_command("samtools view " + bo_s + "_allreads.cf." + pmdref + ".q" + q + ".pmds" + pmd_threshold + "filter.bam | pmdtools.py --deamination --range 30 --CpG > PMD_temp.txt")
        bash_command("R CMD BATCH /data/scripts/plotPMD.R")
        shutil.copy("./PMD_plot.pdf", po_s + "_allreads.cf." + pmdref + ".q" + q + ".pmds" + pmd_threshold + "filter.CpG.plot.pdf")

        ## plot damage for pmd filtered bams - NOT accounting for UDG treatment (and rename plot) ##
        bash_command("samtools view " + bo_s + "_allreads.cf." + pmdref + ".q" + q + ".pmds" + pmd_threshold + "filter.bam | pmdtools.py --deamination --range 30 > PMD_temp.txt")
        bash_command("R CMD BATCH /data/scripts/plotPMD.R")
        shutil.copy("./PMD_plot.pdf", po_s + "_allreads.cf." + pmdref + ".q" + q + ".pmds" + pmd_threshold + "filter.plot.pdf")

    if os.path.isfile("./PMD_plot.pdf"):
        os.remove("./PMD_plot.pdf")

    logfile.close()
    cmdfile.close()
    print "4_pmdtools.py complete."
    exit(0)

