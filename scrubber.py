#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated April 2018     #####
#####       MJJ                 #####
#####################################

# For cleaning up intermediate files from batpipe.runs

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
    print "\n****************\nSCRUBBER\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                 "1. cleans up intermediate files from batpipe run\n"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)




    parser.add_argument('-rawreads', metavar='<rawreads>', help='Location of raw reads',
                        default='/data/raw')
    parser.add_argument('-bc_trim', metavar='<bc_trim>', help='Location of barcode trimmed files',
                        default='/data/barcodetrimmed')
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                    default='/data/adaptertrimmed')

    args = parser.parse_args()
    seqprep_output = args.seqprep_output
    rawreads = args.rawreads
    bc_trim = args.bc_trim


    print "\n*******************WARNING*********************\n"
    print "Do not run this script while anyone is running batpipe on the same machine."
    response = raw_input("Enter Y to proceed, any other key to cancel: ")
    if response != 'Y':
        print "Exiting."
        exit()

    pattern = "bcpos-*"
    files = os.listdir(seqprep_output)
    for name in files:
        if (fnmatch.fnmatch(name, pattern)):
            os.remove(seqprep_output + "/" + name)

    pattern = "nobcpos-*"
    files = os.listdir(seqprep_output)
    for name in files:
        if (fnmatch.fnmatch(name, pattern)):
            os.remove(seqprep_output + "/" + name)

    pattern = "bcpos-*"
    files = os.listdir(bc_trim)
    for name in files:
        if (fnmatch.fnmatch(name, pattern)):
            os.remove(bc_trim + "/" + name)

    pattern = "nobcpos-*"
    files = os.listdir(bc_trim)
    for name in files:
        if (fnmatch.fnmatch(name, pattern)):
            os.remove(bc_trim + "/" + name)

    pattern = "bcpos-*"
    files = os.listdir(rawreads)
    for name in files:
        if (fnmatch.fnmatch(name, pattern)):
            os.remove(rawreads + "/" + name)

    pattern = "nobcpos-*"
    files = os.listdir(rawreads)
    for name in files:
        if (fnmatch.fnmatch(name, pattern)):
            os.remove(rawreads + "/" + name)


    exit(0)





