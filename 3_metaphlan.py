#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated May 2018      #####
#####       MJJ                 #####
#####################################

# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# This script runs 'metaphlan' on fastq files
# Sequence data should already be trimmed by barcode remover and SeqPrep in SeqPrep output, with complexity of reads filtered

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

    print "\n****************\nMETAPHLAN\n****************\n"

    parser = argparse.ArgumentParser(description="# Thie script:\n"
                                                    "This script runs 'metaphlan' on fastq files\n"
                                                    "Sequence data should already be trimmed by barcode remover and SeqPrep in SeqPrep output, with complexity of reads filtered\n"
                                                 "" "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.', required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.add_argument('-mpa_pkl', metavar='<mpa_pkl>', help='mpa_pkl',
                        default='/data/db/db_v20/mpa_v20_m200.pkl')
    parser.add_argument('-bowtie2db', metavar='<bowtie2db>', help='I thoght this was called Bowie2DB and was momentarily happy',
                        default='/data/db/db_v20')
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                        default='/data/adaptertrimmed')
    parser.add_argument('-seqprep_output_in_output', metavar='<seqprep_output_in_output>', help='Prepend output directory to seqprep_output',
                        default=False)


    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    mpa_pkl = args.mpa_pkl
    bowtie2db = args.bowtie2db
    overwrite = bool(args.overwrite)
    verbose = bool(args.verbose)
    seqprep_output_orig = args.seqprep_output
    seqprep_output_in_output = bool(args.seqprep_output_in_output)

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    newdir = wd + "/Metaphlan"
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
    metadir = newdir


    today = datetime.date.today()
    logfilename = wd + "/out.metaphlan." + str(today) + ".log"
    print "Logging to: ", logfilename


    logfile = open(logfilename, 'w')

    logfile.write("Arguments used:\n")
    logfile.write("__________________________________________:\n")
    for arg in vars(args):
        logfile.write(arg)
        logfile.write("\t")
        logfile.write(str(getattr(args, arg)))
        logfile.write("\n")


    cmdfile = open("3_cmds", 'w')

    bcin = open(bcfile, 'r')
    bc = []
    for bcinline in bcin:
        bc.append(bcinline)

    bclength = len(bc)
    print "Number of entries: ", bclength

    print "\nWorking..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        sample = bccols[1]
        output = wd + "/" + sample
        if seqprep_output_in_output:
            seqprep_output = output + "/" + seqprep_output_orig
        else:
            seqprep_output = seqprep_output_orig
        so_s = seqprep_output + "/" + sample
        mo_s = metadir + "/" + sample

        alllist = []

        if not os.path.isfile(so_s + ".all.SP.fq"):

            #Read in F, R, M, read out as all
            partfilename = so_s + ".R.fq"
            if os.path.isfile(partfilename):
                partfile = open(partfilename, 'r')
                for line in partfile:
                    alllist.append(line)
                partfile.close()

            partfilename = so_s + ".F.fq"
            if os.path.isfile(partfilename):
                partfile = open(partfilename, 'r')
                for line in partfile:
                    alllist.append(line)
                partfile.close()

            partfilename = so_s + ".M.fq"
            if os.path.isfile(partfilename):
                partfile = open(partfilename, 'r')
                for line in partfile:
                    alllist.append(line)
                partfile.close()

            alloutname = so_s + ".all.SP.fq"
            allout = open(alloutname, 'w')
            for line in alllist:
                allout.write(line)
            allout.close()

        bash_command("metaphlan2.py " + so_s + ".all.SP.fq --mpa_pkl " + mpa_pkl + " --bowtie2db " + bowtie2db + " --bt2_ps very-sensitive --bowtie2out " + mo_s + ".bz2  --input_type fastq > " + mo_s + ".profiled_metagenome.txt")


    logfile.close()
    cmdfile.close()
    print "3_metaphlan.py complete."
    exit(0)