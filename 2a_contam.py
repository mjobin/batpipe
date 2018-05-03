#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated April 2018     #####
#####       MJJ                 #####
#####################################

# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# from Kircher pipeline used at University of Tuebingen 2013
# M Kircher. Analysis of high-throughput ancient DNA sequencing data. Methods Mol Biol 840:197-228 (2012).
# Original bash wrapper written by Alissa Mitnik modified by Kelly Harkins, UCSC, 2016

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
    outfile.write(cmd)
    outfile.write("\n\n")
    subp = subprocess.Popen(['/bin/bash', '-c', cmd], stdout=PIPE, stderr=PIPE)
    subp.wait()
    theout = subp.stdout.read()
    if verbose:
        print theout
    logfile.write(theout)
    theerr = subp.stderr.read()
    if verbose:
        print theerr
    logfile.write(theerr)
    return theout




if __name__ == "__main__":

    print "\n****************\nCONTAM\n****************\n"

    parser = argparse.ArgumentParser(description="# :\n"
                                                    "from Kircher pipeline used at University of Tuebingen 2013s\n"
                                                    "M Kircher. Analysis of high-throughput ancient DNA sequencing data. Methods Mol Biol 840:197-228 (2012).\n"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.', required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-mt311', metavar='<mt311>', help='mt311',
                        default='/data/genomes/mt311.fna')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)


    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    mt311 = args.mt311
    overwrite = bool(args.overwrite)
    verbose = bool(args.verbose)


    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    today = datetime.date.today()
    logfilename = wd + "/out.contam." + str(today) + ".log"
    print "Logging to: ", logfilename


    logfile = open(logfilename, 'w')

    outfile = open("2a_cmds", 'w')


    newdir = wd + "/ccheck_results"
    if overwrite:
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
    else:
        if os.path.exists(newdir):
            print "ERROR: Directory " + newdir + " exists and overwrite not set to true. Exiting."
            exit(1)
        else:
            os.mkdir(newdir)


    rawfilename = "./ccheck_results/results.raw.contamination.txt"
    forcedfilename = "./ccheck_results/results.forced.raw.contamination.txt"
    rawfile = open(rawfilename, 'w')
    forcedfile = open(forcedfilename, 'w')
    topline = "Sample" + "\t" +  "Strongly diagnostic positions" + "\t" +  "Fragments Orig" + "\t" +  "Fragments Contam" + "\t" + "Fragments total (orig+contam)" + "\t" + "Contamination rate\n"
    rawfile.write(topline)
    forcedfile.write(topline)


    bcin = open(bcfile, 'r')
    bc = []
    for bcinline in bcin:
        bc.append(bcinline)

    bclength = len(bc)
    print "Number of entries: ", bclength
    rawdic = {}
    forceddic = {}
    print "\nWorking..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        sample = bccols[1]
        output =wd + "/" + sample
        mia_output = output + "/MIA_output"

        # get the highest numbered maln file in the MIA output folder

        pattern = sample + "*.maln*"
        files = os.listdir(mia_output)
        malnname = None
        oldmalnint = 0
        for name in files:
            if (fnmatch.fnmatch(name, pattern)):
                malncols = name.split(".")
                malnlen = len(malncols)
                malnlast = malncols[malnlen-1]
                if malnlast.isdigit():
                    malnint = int(malnlast)
                    if not malnname:
                        malnname = name
                    elif malnint > oldmalnint:
                        malnname = name
                        oldmalnint = malnint

        rawout = bash_command("ccheck -a -r " + mt311 + " " + mia_output + "/" + malnname)
        rawlines = rawout.splitlines()

        rl = len(rawlines)
        diag = "N/A"
        fragment_orig = "N/A"
        fragment_contam = "N/A"
        cont_rate = "???"
        fragment_total = "N/A"
        if rl >= 11:
            diag = rawlines[7].split(":")[1]
            fragment_orig = rawlines[9].split(":")[1]
            fragment_contam = rawlines[10].split(":")[1]
            if fragment_orig.isdigit() and fragment_contam.isdigit():
                fragment_total_f = int(fragment_orig) + int(fragment_contam)

        rawoutline = sample + "\t" + diag + "\t" +  fragment_orig + "\t" + fragment_contam + "\t" +  fragment_total + "\t" +  cont_rate + "\n"
        rawfile.write(rawoutline)
        rawdic[sample] = rawoutline

        forcedout = bash_command("ccheck -a -F -r " + mt311 + " " + mia_output + "/" + malnname)
        forcedlines = forcedout.splitlines()

        frl = len(forcedlines)

        diag_f = "N/A"
        fragment_orig_f = "N/A"
        fragment_contam_f = "N/A"
        cont_rate_f = "???"
        fragment_total_f = "N/A"
        if frl >= 11:
            diag_f = forcedlines[7].split(":")[1]
            # # diag_f = cols[1]
            fragment_orig_f = forcedlines[9].split(":")[1]
            fragment_contam_f = forcedlines[10].split(":")[1]
            if fragment_orig_f.isdigit() and fragment_contam_f.isdigit():
                fragment_total_f = int(fragment_orig_f) + int(fragment_contam_f)


        forcedoutline = sample + "\t" + diag_f + "\t" + fragment_orig_f + "\t" + fragment_contam_f + "\t" + fragment_total_f + "\t" +  cont_rate_f + "\n"
        forcedfile.write(forcedoutline)
        forceddic[sample] = forcedoutline


    rawsortfilename = "./ccheck_results/results.raw.contamination.sorted.txt"
    forcedsortfilename = "./ccheck_results/results.forced.raw.contamination.sorted.txt"
    rawsortfile = open(rawsortfilename, 'w')
    forcedsortfile = open(forcedsortfilename, 'w')
    rawsortfile.write(topline)
    forcedsortfile.write(topline)

    for x in sorted(rawdic.keys()):
        rawsortfile.write(x)
    for x in sorted(forceddic.keys()):
        forcedsortfile.write(x)


    logfile.close()
    outfile.close()
    rawfile.close()
    forcedfile.close()
    rawsortfile.close()
    forcedsortfile.close()
    print "2a_contam.py complete."
    exit(0)