#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated April 2018      #####
#####       MJJ                 #####
#####################################

# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
## contammix-1.0-10 from Philip Johnson:
###    Caveats:
###    --------
### The key assumption is that contamination is <50%, and thus that the consensus reflects the authentic sequence.  As always with MCMC, achieving convergence can be fiddly.  I suggest always graphically checking for convergence in addition to looking at the Gelman-Rubin diagnostic.  For some datasets, I have had to tweak the "alpha" hyperparameter (this is an option to estimate.R).


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

    print "\n****************\nCONTAMMIX\n****************\n"

    parser = argparse.ArgumentParser(description="# :\n"
                                                    "The key assumption is that contamination is <50%, and thus that the consensus reflects the authentic sequence.  As always with MCMC, achieving convergence can be fiddly.  I suggest always graphically checking for convergence in addition to looking at the Gelman-Rubin diagnostic.  For some datasets, I have had to tweak the alpha hyperparameter (this is an option to estimate.R).\n"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.', required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-mt311', metavar='<mt311>', help='mt311',
                        default='/data/genomes/mt311.fna')
    parser.add_argument('-mia_ref', metavar='<mia_ref>', help='mia ref',
                        default='rCRS')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                        default='/data/adaptertrimmed')
    parser.add_argument('-seqprep_output_in_output', metavar='<seqprep_output_in_output>', help='Prepend output directory to seqprep_output',
                        default=False)
    parser.add_argument('-threads', metavar='<threads>', help='To speed up analysis',
                        default="23")


    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    mt311 = args.mt311
    mia_ref = args.mia_ref
    overwrite = bool(args.overwrite)
    verbose = bool(args.verbose)
    seqprep_output_orig = args.seqprep_output
    seqprep_output_in_output = bool(args.seqprep_output_in_output)
    threads = args.threads


    cols311 = mt311.split("/")
    last311 = cols311[len(cols311)-1]
    namecols311 = last311.split(".")
    name311 = namecols311[0]

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    today = datetime.date.today()
    logfilename = wd + "/out.contammix." + str(today) + ".log"
    print "Logging to: ", logfilename


    logfile = open(logfilename, 'w')

    logfile.write("Arguments used:\n")
    logfile.write("__________________________________________:\n")
    for arg in vars(args):
        logfile.write(arg)
        logfile.write("\t")
        logfile.write(str(getattr(args, arg)))
        logfile.write("\n")


    cmdfile = open("4b_cmds", 'w')

    newdir = wd + "/contammix"
    if overwrite:
        if os.path.exists(newdir):
            shutil.rmtree(newdir)
        os.mkdir(newdir)
    else:
        if os.path.exists(newdir):
            print "ERROR: Directory " + newdir + " exists and overwrite not set to true. Exiting."
            exit()
        else:
            os.mkdir(newdir)

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
        in_sample = bccols[0]
        out_sample = bccols[1]
        output = wd + "/" + out_sample
        mia_output = output + "/MIA_output"
        mo_s = mia_output + "/" + out_sample
        if seqprep_output_in_output:
            seqprep_output = output + "/" + seqprep_output_orig
        else:
            seqprep_output = seqprep_output_orig
        so_s = seqprep_output + "/" + in_sample
        cmix_output = wd + "/contammix"
        co_s = cmix_output + "/" + out_sample

        if not os.path.isfile(so_s + ".all.SP.cf.fastq"):
            with gzip.open(so_s + ".all.SP.cf.fastq.gz", 'rb') as f_in, open(so_s + ".all.SP.cf.fastq", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)


        bash_command("cat " + mo_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.mia_consensus.fasta " + mt311 + " > " + mo_s + ".allreads." + mia_ref + "." + name311 + ".fasta")
        bash_command("mafft --auto " + mo_s + ".allreads." + mia_ref + "." + name311 + ".fasta > " + co_s + ".allreads." + mia_ref + "." + name311 + ".mafft.fasta")
        bash_command("bwa index " + mo_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.mia_consensus.fasta")
        bash_command("bwa aln -l 1024 -n 0.02 -t " + threads + " " + mo_s + ".all.SP.cf.rd." + mia_ref +".maln.F.mia_consensus.fasta " + so_s + ".all.SP.cf.fastq > " + co_s + ".all.SP.cf.remapped.sai")
        bash_command("bwa samse " + mo_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.mia_consensus.fasta " + co_s + ".all.SP.cf.remapped.sai " + so_s + ".all.SP.cf.fastq  > " + co_s + ".all.SP.cf.remapped.sam")
        bash_command("samtools view -bSh " + co_s + ".all.SP.cf.remapped.sam > " + co_s + ".all.SP.cf.remapped.bam")
        bash_command("/data/install/contamMix/exec/estimate.R --samFn " + co_s + ".all.SP.cf.remapped.bam --malnFn " + co_s + ".allreads." + mia_ref + "." + name311 + ".mafft.fasta  --figure " + co_s + ".contammix_fig")

        if os.path.isfile(so_s + ".all.SP.cf.fastq") and os.path.isfile(so_s + ".all.SP.cf.fastq" + ".gz"):
            os.remove(so_s + ".all.SP.cf.fastq" + ".gz")
        if os.path.isfile(so_s + ".all.SP.cf.fastq"):
            with open(so_s + ".all.SP.cf.fastq", 'rb') as f_in, gzip.open(so_s + ".all.SP.cf.fastq" + ".gz", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(so_s + ".all.SP.cf.fastq")

    logfile.close()
    cmdfile.close()
    print "2b_contammix.py complete."
