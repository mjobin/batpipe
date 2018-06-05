#!/usr/bin/python


#####################################
#####        HPG Lab            #####
#####    updated April 2018     #####
#####       MJJ                 #####
#####################################


# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# Based on bash scripts by Kelly Harkins and Peter Heintzman
# Script uses two programs:
# 1. bc_bin_clip separates then trims 5' of R1/R2 with correct barcodes
# 2. SeqPrep trim barcodes + adapters from 3' end/inside read when frag shorter than read length

import argparse
from argparse import RawTextHelpFormatter
import progressbar
import os
import sys
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

    print "\n****************\nTRIM BARCODES\n****************\n"

    parser = argparse.ArgumentParser(description="Script uses two programs:\n\t"
                                                "1. bc_bin_clip separates then trims 5' of R1/R2 with correct barcodes\n"
                                                "2. SeqPrep trim barcodes + adapters from 3' end/inside read when frag shorter than read length\n"
                                                "", formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.', required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-bc_bin_clip', metavar='<bcbinclip>', help='Location of BC_BIN_CLIP',
                        default='/data/scripts')
    parser.add_argument('-rawreads', metavar='<rawreads>', help='Location of raw reads',
                        default='/data/raw')
    parser.add_argument('-ref_barcodes', metavar='<ref_barcodes>', help='Location of references barcodes',
                        default='/data/projects/internal_barcodes.txt')
    parser.add_argument('-univ_illumina', metavar='<univ_illumina>', help='Universal Illumina adapter',
                        default='AGATCGGAAG')
    parser.add_argument('-seqprep_min_length', metavar='<seqprep_min_length>', help='Seqprep Min lengthr',
                        default=30)
    parser.add_argument('-seqprep_overlap', metavar='<seqprep_overlap>', help='Seqprep overlap',
                        default=10)
    parser.add_argument('-mismatch', metavar='<mismatch>', help='Mismatch tolerance',
                        default=2)
    parser.add_argument('-threads', metavar='<threads>', help='Max threads',
                        default="8")
    parser.add_argument('-bc_trim', metavar='<bc_trim>', help='Location of barcode trimmed files',
                        default='/data/barcodetrimmed')
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                        default='/data/adaptertrimmed')
    parser.add_argument('-bformat', dest='bformat', help='Barcode format type old (barcodes as sequence not numbers)?.',
                        action='store_true')
    parser.set_defaults(bformat=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)
    parser.add_argument('-fqloc', metavar='<fqloc>', help='Location of fq files',
                        default='/data/adaptertrimmed')
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-nopos', dest='nopos', help='Do not add/check digital positive files,',
                        action='store_true')
    parser.set_defaults(nopos=False)


    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    bcbinclip = args.bc_bin_clip
    rawreads = args.rawreads
    univ_illumina = args.univ_illumina
    ref_barcodes = args.ref_barcodes
    seqprep_min_length = int(args.seqprep_min_length)
    seqprep_overlap = int(args.seqprep_overlap)
    mismatch = int(args.mismatch)
    bc_trim = args.bc_trim
    seqprep_output = args.seqprep_output
    bformat = bool(args.bformat)
    overwrite = bool(args.overwrite)
    verbose = bool(args.verbose)
    fqloc = args.fqloc
    threads = args.threads
    nopos = bool(args.nopos)

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    today = datetime.date.today()
    rightnow = str(datetime.datetime.now().time())



    logfilename = "out.trim.m" + str(mismatch) + "." + str(today) + ".log"
    print "Logging to:", logfilename
    logfile = open(logfilename, 'w')

    cmdfile = None
    if nopos:
        cmdfile = open("0_nopos_cmds", 'w')
    else:
        cmdfile = open("0_cmds", 'w')

    logfile.write("Parameters used: \nDebarcoding mismatch tolerance: " + str(mismatch) + "\n")
    logfile.write("Seqprep Minimum Length: " + str(seqprep_min_length) + "\n")
    logfile.write("Seqprep overlap: " + str(seqprep_overlap) + "\n")
    logfile.write("-------------------------\n")

    refbcs = {}
    refbcin = open(ref_barcodes, 'r')
    for refinline in refbcin:
        refcols = refinline.split()
        if len(refcols) == 3: #only use lines of correct size
            refbcs[refcols[0]] = refcols

    bcin = open(bcfile, 'r')
    bc = []
    barcodespresent = True
    for bcinline in bcin:
        bccols = bcinline.split()
        if len(bccols) > 3:
            barcodespresent = True
        else:
            barcodespresent = False
        bc.append(bcinline)

    if not nopos:
        shutil.copy("/data/raw/bcpos_S00_L00_R1_001.fastq.gz", "/data/raw/bcpos" + rightnow + "_S00_L00_R1_001.fastq.gz")
        shutil.copy("/data/raw/bcpos_S00_L00_R2_001.fastq.gz", "/data/raw/bcpos" + rightnow + "_S00_L00_R2_001.fastq.gz")
        shutil.copy("/data/raw/nobcpos_S00_L00_R1_001.fastq.gz", "/data/raw/nobcpos" + rightnow + "_S00_L00_R1_001.fastq.gz")
        shutil.copy("/data/raw/nobcpos_S00_L00_R2_001.fastq.gz", "/data/raw/nobcpos" + rightnow + "_S00_L00_R2_001.fastq.gz")


        if bformat:
            bcposline = "bcpos"+rightnow+"	bcpos"+rightnow+"	TCGAACA	AGCACAT ATGTGCT TGTTCGA"
        else:
            bcposline = "bcpos"+rightnow+"	bcpos"+rightnow+"	198	205"
        nobcposline = "nobcpos"+rightnow+"	nobcpos"+rightnow+"		"
        if barcodespresent:
            bc.append(bcposline)
        else:
            bc.append(nobcposline)

    bclength = len(bc)
    print "Number of entries: ", bclength

    missingfilepatterns = []

    print "Checking for files..."
    for i in range(bclength):
        filestorezip = []
        bcline = bc[i]
        bccols = bcline.split()
        in_sample = bccols[0]
        r1pattern = in_sample + "_*_L00*_R1_001.fastq"
        logfile.write("Matching to " + r1pattern + "\n")
        files = os.listdir(rawreads)
        r1name = None
        for name in files:
            if (fnmatch.fnmatch(name, r1pattern)):
                continue
        if not r1name:
            r1pattern = in_sample + "_*_L00*_R1_001.fastq.gz"
            for name in files:
                if (fnmatch.fnmatch(name, r1pattern)):
                    nogz = name.rsplit(".", 1)[0]
                    r1name = rawreads + "/" + nogz
        if not r1name:
            print "Cannot find file matching " + rawreads + "/" + r1pattern + "."
            missingfilepatterns.append(r1pattern)
        filestorezip.append(r1name)

        logfile.write("Gunzipping R2 " + in_sample + "\n")
        r2pattern = in_sample + "_*_L00*_R2_001.fastq"
        logfile.write("Matching to " + r2pattern + "\n")
        files = os.listdir(rawreads)
        r2name = None
        for name in files:
            if (fnmatch.fnmatch(name, r2pattern)):
                r2name = rawreads + "/" + name
                continue
        if not r2name:
            r2pattern = in_sample + "_*_L00*_R2_001.fastq.gz"
            for name in files:
                if (fnmatch.fnmatch(name, r2pattern)):
                    nogz = name.rsplit(".", 1)[0]
                    r2name = rawreads + "/" + nogz
        if not r2name:
            print "Cannot find file matching " + rawreads + "/" + r2pattern + "."
            missingfilepatterns.append(r2pattern)
        filestorezip.append(r2name)

    if len(missingfilepatterns) > 0:
        exit(1)


    print "\nTrimming barcodes..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        filestorezip = []
        bcline = bc[i]
        bccols = bcline.split()
        in_sample = bccols[0]
        out_sample = bccols[1]
        barcodespresent = True
        newdir = out_sample
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
        logfile.write(in_sample + "\t" + out_sample + "\n")
        if len(bccols)>3:
            if bformat:
                bc1 = bccols[2]
                if bc1.isdigit():
                    print "ERROR. -bformat invoked but " + p5 + " is numeric. Check what barcode format you are using!\n ABORT! \n ABORT!\n ABORRRRRRRRRT!"
                    exit(1)
                bc2 = bccols[3]
                sp1 = bccols[4]
                sp2 = bccols[5]
                sp_r1 = sp1 + univ_illumina
                sp_r2 = sp2 + univ_illumina
            else:
                p5 = bccols[2]
                if not p5.isdigit():
                    print "ERROR. -bformat NOT invoked but " + p5 + " is not numeric. Check what barcode format you are using! \n ABORT! \n ABORT!\n ABORRRRRRRRRT!"
                    exit(1)
                p7 = bccols[3]
                bc1 = refbcs[p5][1]
                bc2 = refbcs[p7][1]
                sp1 = refbcs[p7][2]
                sp2 = refbcs[p5][2]
                sp_r1 = sp1 + univ_illumina
                sp_r2 = sp2 + univ_illumina
        else:
            barcodespresent = False
            sp_r1 = univ_illumina
            sp_r2 = univ_illumina


#Gunzip only if corresponding fastq file not already there
        logfile.write("Gunzipping R1 " + in_sample + "\n")
        r1pattern = in_sample + "_*_L00*_R1_001.fastq"
        logfile.write("Matching to " + r1pattern + "\n")
        files = os.listdir(rawreads)
        r1name = None
        for name in files:
            if (fnmatch.fnmatch(name, r1pattern)):
                r1name = rawreads + "/" + name
                filestorezip.append(r1name)
                continue
        if not r1name:
            r1pattern = in_sample + "_*_L00*_R1_001.fastq.gz"
            for name in files:
                if (fnmatch.fnmatch(name, r1pattern)):
                    nogz = name.rsplit(".", 1 )[ 0 ]
                    r1name = rawreads + "/" + nogz
                    filestorezip.append(r1name)
                    with gzip.open(rawreads + "/" + name, 'rb') as f_in, open(r1name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

        if not r1name:
            print "Cannot find file matching " + rawreads + "/" + r2pattern + ". Exiting."
            exit(1)

        logfile.write("Gunzipping R2 "+ in_sample + "\n")
        r2pattern = in_sample + "_*_L00*_R2_001.fastq"
        logfile.write("Matching to " + r2pattern + "\n")
        files = os.listdir(rawreads)
        r2name = None
        for name in files:
            if (fnmatch.fnmatch(name, r2pattern)):
                r2name = rawreads + "/" + name
                filestorezip.append(r2name)
                continue
        if not r2name:
            r2pattern = in_sample + "_*_L00*_R2_001.fastq.gz"
            for name in files:
                if (fnmatch.fnmatch(name, r2pattern)):
                    nogz = name.rsplit(".", 1)[0]
                    r2name = rawreads + "/" + nogz
                    filestorezip.append(r2name)
                    with gzip.open(rawreads + "/" + name, 'rb') as f_in, open(r2name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
        if not r2name:
            print "Cannot find file matching " + rawreads + "/" + r2pattern + ". Exiting."
            exit(1)


        if barcodespresent: # if sample has barcodes, use this part
            logfile.write(in_sample + " has barcodes " + bc1 + " and " + bc2 + "\n")
            logfile.write("Starting bc_bin_clip for " + in_sample + "\n")
            bcbinclipcmd = "python " + bcbinclip + "/bc_bin_clip -m " + str(mismatch) + " " + bc1 + " " + bc2 + " " + r1name + " " + r2name + " -p " + bc_trim + "/" + in_sample + ".woBC"
            bash_command(bcbinclipcmd)




            logfile.write("Starting Seqprep2 " + in_sample + "\n")
            bc_os = bc_trim + "/" + in_sample

            sp_os = seqprep_output + "/" + in_sample

            sr1pattern = in_sample + "*_R1.fastq"
            sr1upattern = in_sample + "*_R1_unmatched.fastq"
            files = os.listdir(bc_trim)
            sr1name = None
            for name in files:
                if (fnmatch.fnmatch(name, sr1pattern)):
                    sr1name = bc_trim + "/" + name
                    filestorezip.append(sr1name)
                if (fnmatch.fnmatch(name, sr1upattern)):
                    sr1uname = bc_trim + "/" + name
                    filestorezip.append(sr1uname)

            sr2pattern = in_sample + "*_R2.fastq"
            sr2upattern = in_sample + "*_R2_unmatched.fastq"
            files = os.listdir(bc_trim)
            sr2name = None
            for name in files:
                if (fnmatch.fnmatch(name, sr2pattern)):
                    sr2name = bc_trim + "/" + name
                    filestorezip.append(sr2name)
                if (fnmatch.fnmatch(name, sr2upattern)):
                    sr2uname = bc_trim + "/" + name
                    filestorezip.append(sr2uname)

            seqprepcmd = "SeqPrep2 -f " + sr1name + " -r " + sr2name + " -1 " + sp_os + ".F.fq.gz -2 " + sp_os + ".R.fq.gz -s " + sp_os + ".M.fq.gz -A " + sp_r1 + " -B " + sp_r2 + " -o " + str(seqprep_overlap) + " -L " + str(seqprep_min_length) + " -d 1 -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT >& " + seqprep_output + "/SP." + in_sample + ".stderr.txt"
            bash_command(seqprepcmd)

            logfile.write("SeqPrep complete" + "\n")
            logfile.write("-------------------------" + "\n")





        else:  # if sample has NO barcodes, use this part
            logfile.write("Meyer-Kircher library protocol, no internal barcodes for " + in_sample + "\n")
            rr_is = rawreads + "/" + in_sample

            sp_os = seqprep_output + "/" + in_sample
            seqprepcmd = "SeqPrep2 -f " + rr_is + "_*_L00*_R1_001.fastq -r " + rr_is + "_*_L00*_R2_001.fastq -1 " + sp_os + ".F.fq.gz -2 " + sp_os + ".R.fq.gz -s " + sp_os + ".M.fq.gz -A " + sp_r1 + " -B " + sp_r2 + " -o " + str(seqprep_overlap) + " -L " + str(seqprep_min_length) + " -d 1 -C ATCTCGTATGCCGTCTTCTGCTTG -D GATCTCGGTGGTCGCCGTATCATT >& " + seqprep_output + "/SP." + in_sample + ".stderr.txt"
            bash_command(seqprepcmd)
            logfile.write("SeqPrep complete" + "\n")
            logfile.write("-------------------------" + "\n")

            r1pattern = in_sample + "_*_L00*_R1_001.fastq"
            logfile.write ("Matching to " + r1pattern + "\n")
            files = os.listdir(rawreads)
            r1name = None
            for name in files:
                if (fnmatch.fnmatch(name, r1pattern)):
                    r1name = rawreads + "/" + name
                    if not os.path.isfile(r1name):
                        with open(r1name, 'rb') as f_in, gzip.open(r1name + ".gz", 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)

            r2pattern = in_sample + "_*_L00*_R2_001.fastq"
            logfile.write ("Matching to " + r2pattern + "\n")
            files = os.listdir(rawreads)
            r2name = None
            for name in files:
                if (fnmatch.fnmatch(name, r2pattern)):
                    r2name = rawreads + "/" + name
                    if not os.path.isfile(r2name):
                        with open(r2name, 'rb') as f_in, gzip.open(r2name + ".gz", 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)


        #Check digital positive
        sp_file = seqprep_output + "/SP." + in_sample + ".stderr.txt"
        spc = []
        if os.path.isfile(sp_file):
            spin = open(sp_file, 'r')
            for spinline in spin:
                spc.append(spinline)

        else:
            print "SP No file matching " + sp_file + " found. Exiting."
            exit(1)


        reads_debarcoded = int(str(spc[1].split()[2]).strip())


        if in_sample == "bcpos"+rightnow:
            if reads_debarcoded != 208306:
                print "ERROR barcode positive control reads debarcoded should read 208306 but reads " + str(reads_debarcoded)
                exit(1)
        elif in_sample == "nobcpos"+rightnow:
            if reads_debarcoded != 247087:
                print "ERROR no barcode positive control reads debarcoded should read 247087 " + str(reads_debarcoded)
                exit(1)


        #Rezip files
        for ftrz in filestorezip:
            if os.path.isfile(ftrz) and os.path.isfile(ftrz + ".gz"):
                os.remove(ftrz + ".gz")
            with open(ftrz, 'rb') as f_in, gzip.open(ftrz + ".gz", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(ftrz)



        # makes symbolic links of trimmed files to project working directories
        newdir = out_sample + "/SeqPrep2_output"
        sp_output_s = out_sample + "/SeqPrep2_output"
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

        newdir = out_sample + "/BC_trimmed"
        bc_trim_s = out_sample + "/BC_trimmed"
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

        if barcodespresent:
            pattern = in_sample + "*_R1.fastq.gz"
            files = os.listdir(bc_trim)
            for name in files:
                if (fnmatch.fnmatch(name, pattern)):
                    os.symlink(bc_trim + "/" + name, bc_trim_s + "/" + name)

            pattern = in_sample + "*_R2.fastq.gz"
            files = os.listdir(bc_trim)
            for name in files:
                if (fnmatch.fnmatch(name, pattern)):
                    os.symlink(bc_trim + "/" + name, bc_trim_s + "/" + name)

        pattern = in_sample + "*.R.fq.gz"
        files = os.listdir(seqprep_output)
        for name in files:
            if (fnmatch.fnmatch(name, pattern)):
                os.symlink(seqprep_output + "/" + name, sp_output_s + "/" + name)

        pattern = in_sample + "*.F.fq.gz"
        files = os.listdir(seqprep_output)
        for name in files:
            if (fnmatch.fnmatch(name, pattern)):
                os.symlink(seqprep_output + "/" + name, sp_output_s + "/" + name)

        pattern = in_sample + "*.M.fq.gz"
        files = os.listdir(seqprep_output)
        for name in files:
            if (fnmatch.fnmatch(name, pattern)):
                os.symlink(seqprep_output + "/" + name, sp_output_s + "/" + name)

        pattern = in_sample + "*.stderr.txt"
        files = os.listdir(seqprep_output)
        for name in files:
            if (fnmatch.fnmatch(name, pattern)):
                os.symlink(seqprep_output + "/" + name, sp_output_s + "/" + name)

    if not nopos:
        delfilelist = []
        # Delete barcode positive files
        if barcodespresent:
            delfilelist.append(rawreads + "/bcpos" + rightnow + "_S00_L00_R1_001.fastq.gz")
            delfilelist.append(rawreads + "/bcpos" + rightnow + "_S00_L00_R2_001.fastq.gz")

            delfilelist.append(bc_trim + "/bcpos" + rightnow + ".woBC_R1.fastq.gz")
            delfilelist.append(bc_trim + "/bcpos" + rightnow + ".woBC_R2.fastq.gz")
            delfilelist.append(bc_trim + "/bcpos" + rightnow + ".woBC_R1_unmatched.fastq.gz")
            delfilelist.append(bc_trim + "/bcpos" + rightnow + ".woBC_R2_unmatched.fastq.gz")
            delfilelist.append(seqprep_output + "/bcpos" + rightnow + ".F.fq.gz")
            delfilelist.append(seqprep_output + "/bcpos" + rightnow + ".R.fq.gz")
            delfilelist.append(seqprep_output + "/bcpos" + rightnow + ".M.fq.gz")
        else:
            delfilelist.append(rawreads + "/nobcpos" + rightnow + "_S00_L00_R1_001.fastq.gz")
            delfilelist.append(rawreads + "/nobcpos" + rightnow + "_S00_L00_R2_001.fastq.gz")

            delfilelist.append(seqprep_output + "/nobcpos" + rightnow + ".F.fq.gz")
            delfilelist.append(seqprep_output + "/nobcpos" + rightnow + ".R.fq.gz")
            delfilelist.append(seqprep_output + "/nobcpos" + rightnow + ".M.fq.gz")

            delfilelist.append(bc_trim + "/nobcpos" + rightnow + ".woBC_R1.fastq.gz")
            delfilelist.append(bc_trim + "/nobcpos" + rightnow + ".woBC_R2.fastq.gz")
            delfilelist.append(bc_trim + "/nobcpos" + rightnow + ".woBC_R1_unmatched.fastq.gz")
            delfilelist.append(bc_trim + "/nobcpos" + rightnow + ".woBC_R2_unmatched.fastq.gz")
        delfilelist.append("posdummy.txt")

    logfile.close()
    cmdfile.close()
    print "0_trim_barcodes complete."
    exit(0)
