#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated April 2018     #####
#####       MJJ                 #####
#####################################

# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# print a tab delimited file with basic stats.
# Numbers of raw reads, seqprep results, bwa > summary.txt
# Based off of K Harkins's  bash script
# this works when certain naming conventions apply:

## 1. Seqprep output is named SP2.libraryID.stderr.txt
## 2. raw fastq file is named iPCR##-LibraryNumber_S0_L0*_R1_001.fastq
## 3. sorted BAM files are named: L001.M.cf.hg19.q*.s.bam or fileprefix_allreads.cf.hg19.q*.s.bam

import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import datetime
import gzip
import shutil
import fnmatch
import subprocess
import progressbar
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

    print "\n****************\nRESULTS\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                 "Prints a tab delimited file with basic stats.\n\t"
                                                 "Numbers of raw reads, seqprep results, bwa > summary.txt\n\t"
                                                 "This works when certain naming conventions apply:\n\t"
                                                 "1. Seqprep output is named SP2.libraryID.stderr.txt\n\t"
                                                 "2. raw fastq file is named iPCR##-LibraryNumber_S0_L0*_R1_001.fastq\n\t"
                                                 "3. sorted BAM files are named: L001.M.cf.hg19.q*.s.bam or fileprefix_allreads.cf.hg19.q*.s.bam\n\t"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.',
                        required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-raw', metavar='<raw>', help='Location of raw reads',
                        default='/data/raw')
    parser.add_argument('-bwa_ref', metavar='<bwa_ref>', help='bwa_ref',
                        default='hg19')
    parser.add_argument('-mito_ref', metavar='<mito_ref>', help='mito_ref',
                        default='rCRS')
    parser.add_argument('-wobc_output', metavar='<wobc_output>', help='wobc_output',
                        default='/data/barcodetrimmed')
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                        default='/data/adaptertrimmed')
    parser.add_argument('-seqprep_output_in_output', metavar='<seqprep_output_in_output>',
                        help='Prepend output directory to seqprep_output',
                        default=False)
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-leehom', dest='leehom', help='Use leeHom instead of seqprep.',
                        action='store_true')
    parser.set_defaults(leehom=False)
    parser.add_argument('-nomia', dest='nomia', help='Do not run MIA.',
                        action='store_true')
    parser.set_defaults(nomia=False)
    parser.add_argument('-refs', dest='refs', nargs='+', default=[],
                        help='List of reference sequences other than hg19 and rCRS.')
    parser.add_argument('-bformat', dest='bformat', help='Barcode format type old (barcodes as sequence not numbers)?.',
                        action='store_true')
    parser.set_defaults(bformat=False)
    parser.add_argument('-pmd_threshold', metavar='<pmd_threshold>', help='PMDtools threshold',
                        default="3")

    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    raw = args.raw
    bwa_ref = args.bwa_ref
    mito_ref = args.mito_ref
    seqprep_output_orig = args.seqprep_output
    seqprep_output_in_output = bool(args.seqprep_output_in_output)
    wobc_output = args.wobc_output
    verbose = bool(args.verbose)
    leehom = bool(args.leehom)
    nomia = bool(args.nomia)
    refs = args.refs
    bformat = bool(args.bformat)
    pmd_threshold = args.pmd_threshold

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    today = datetime.date.today()
    logfilename = wd + "/results." + str(today) + ".log"
    print "Logging to: ", logfilename
    logfile = open(logfilename, 'w')

    today = datetime.date.today()

    cmdfile = open("results_cmds", 'w')

    arefs = ["/data/genomes/hg19.fa", "/data/genomes/rCRS.fas"]
    for ref in refs:
        arefs.append(ref)

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


    bclength = len(bc)
    print "Number of entries: ", bclength

    refdic = {}
    logfile.write("Reference sequences used: \n")
    for ref in arefs:
        refname = os.path.basename(ref)
        filebase, fileext = os.path.splitext(refname)
        refdic[filebase] = ref
        logfile.write(ref + "\n")
    logfile.write("-------------------\n")

    for key, value in refdic.iteritems():
        sumfilename = wd + "/summary." + str(today) + "." + key + ".txt"
        print "\nResult for " + key + " to " + sumfilename
        sumfile = open(sumfilename, 'w')

        # 1
        sumfile.write("Date   ")
        sumfile.write("	")
        # 2
        sumfile.write(" SeqID   ")
        sumfile.write("	")
        # 3
        sumfile.write(" Library ID   ")
        sumfile.write("	")
        # 4
        sumfile.write(" Sample   ")
        sumfile.write("	")
        # 5
        sumfile.write(" total reads ")
        sumfile.write("	")
        # 6
        sumfile.write(" total after BC removal   ")
        sumfile.write("	")
        # 7
        sumfile.write(" % reads with BC   ")
        sumfile.write("	")
        # 8
        sumfile.write(" pairs processed   ")
        sumfile.write("	")
        # 9
        sumfile.write(" pairs merged   ")
        sumfile.write("	")
        # 10
        sumfile.write(" pairs with adapters   ")
        sumfile.write("	")
        # 11
        sumfile.write(" pairs discarded   ")
        sumfile.write("	")
        # 12
        sumfile.write(" M reads used   ")
        sumfile.write("	")
        # 13
        sumfile.write(" flagstats M reads mapped  ")
        sumfile.write("	")
        # 14
        sumfile.write(" M reads mapped, rmdup   ")
        sumfile.write("	")
        # 15
        sumfile.write(" ALL reads, rmdup, uniqly mapping   ")
        sumfile.write("	")
        # 16
        sumfile.write(" Avg length all mapped reads   ")
        sumfile.write("	")
        # 17
        sumfile.write(" Quality used   ")
        sumfile.write("	")

        # 18
        sumfile.write("% reads with barcodes (no mismatches) ")
        sumfile.write("	")
        # 19
        sumfile.write("% not duplicate ")
        sumfile.write("	")
        # 20
        sumfile.write("q20 %endogenous ")
        sumfile.write("	")
        # 21
        sumfile.write("% merged ")
        sumfile.write("	")

        # 22
        sumfile.write("5' damage ")
        sumfile.write("	")
        # 23
        sumfile.write("3' damage ")
        sumfile.write("	")
        # 24
        sumfile.write("analysis notes ")
        sumfile.write("	")
        # 25
        sumfile.write("location")
        sumfile.write("	")

        sumfile.write("\n")

        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample

            # Get lib name if its there
            in_samplecols = in_sample.split("-")
            if len(in_samplecols) >= 2:
                lib = in_samplecols[1]
            else:
                lib = "???"

            if seqprep_output_in_output:
                seqprep_output = output + "/" + seqprep_output_orig
            else:
                seqprep_output = seqprep_output_orig
            so_s = seqprep_output + "/" + in_sample

            rawgzpattern = in_sample + "*_L00*_R1_001.fastq.gz"
            files = os.listdir(raw)
            r1name = None
            for name in files:
                if (fnmatch.fnmatch(name, rawgzpattern)):
                    r1name = raw + "/" + name
                    r1outname = os.path.splitext(r1name)[0]
                    if not os.path.isfile(r1outname):
                        with gzip.open(r1name, 'rb') as f_in, open(r1outname, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)

            # Will take the FIRST matching file
            rawfile = None
            rawpattern = in_sample + "*_L00*_R1_001.fastq"
            files = os.listdir(raw)
            for name in files:
                if (fnmatch.fnmatch(name, rawpattern)):
                    rawfile = raw + "/" + name
                    break
            if not os.path.isfile(rawfile):
                print "No file matching " + rawpattern + " found. Exiting."
                exit(1)

            rawin = open(rawfile, 'r')
            rawc = []
            for rawcinline in rawin:
                rawc.append(rawcinline)

            rawclength = len(rawc)
            rawreads = str(float(rawclength) / 4.0)

            if leehom:
                lh_file = seqprep_output + "/LH." + in_sample + ".stderr.txt"
                lhc = []
                if os.path.isfile(lh_file):
                    lhin = open(lh_file, 'r')
                    for lhinline in lhin:
                        lhc.append(lhinline)

                else:
                    print "No file matching " + lh_file + " found. Exiting."
                    exit(1)

                rdbname = wobc_output + "/" + in_sample + ".woBC_R1.fastq"
                if not os.path.isfile(rdbname):
                    with gzip.open(rdbname + ".gz", 'rb') as f_in, open(rdbname, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                reads_debarcoded = int(bash_command("wc -l " + rdbname).split()[0]) / 4

                mname = so_s + ".M.fq"
                if not os.path.isfile(mname):
                    with gzip.open(mname + ".gz", 'rb') as f_in, open(mname, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                merged = int(bash_command("wc -l " + mname).split()[0]) / 4
                noperc = lhc[1].split("\t")[0]
                splitspc = noperc.split(" ")
                wadapters = splitspc[2]  #
                discarded = "N/A"

            else:

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
                merged = int(spc[2].split()[2].strip())
                wadapters = int(spc[3].split()[3].strip())
                discarded = int(spc[4].split()[2].strip())

            rdbf = float(reads_debarcoded)
            if rdbf == 0:
                percentmerged = "N/A"
                print "Reads_debarcoded is 0. Are you sure you have the right indexes?"
            else:
                percentmerged = str((100.0 * float(merged) / rdbf))

            rrf = float(rawreads)
            if rrf == 0:
                percentbarcodes = "N/A"
            else:
                percentbarcodes = str((100.0 * rdbf) / rrf)

            bwa_output1 = output + "/BWA_" + key
            bo1_s = bwa_output1 + "/" + out_sample

            mpattern = out_sample + ".M.cf.*.q*.s.bam"
            files = os.listdir(bwa_output1)
            mname = None
            for name in files:
                if (fnmatch.fnmatch(name, mpattern)):
                    mname = name

            if not mname:
                print "ERROR: Cannot find match to" + mpattern + " Exiting. "
                exit(1)

            allpattern = out_sample + "_allreads.cf.*.q*.s.bam"
            files = os.listdir(bwa_output1)
            allname = None
            for name in files:
                if (fnmatch.fnmatch(name, allpattern)):
                    allname = name
            if not allname:
                print "ERROR: Cannot find " + allpattern + " Exiting. "
                exit(1)

            merged_filtered_bam = bash_command("samtools view -c " + bo1_s + ".M.cf.*.q*.s.bam").strip()
            if merged_filtered_bam:
                merged_filtered_bam = int(merged_filtered_bam)
            rd_merged_filtered_bam = bash_command("samtools view -c " + bo1_s + ".M.cf.*.q*.s.bam").strip()
            if rd_merged_filtered_bam:
                rd_merged_filtered_bam = int(rd_merged_filtered_bam)
            all_filtered_bam = bash_command("samtools view -c " + bo1_s + "_allreads.cf.*.q*.s.bam").strip()
            if all_filtered_bam:
                all_filtered_bam = int(all_filtered_bam)
            rd_all_filtered_bam = bash_command("samtools view -c " + bo1_s + "_allreads.cf.*.q*.s.bam").strip()
            if rd_all_filtered_bam:
                rd_all_filtered_bam = int(rd_all_filtered_bam)

            uniq_all_filtered_bam = bash_command(
                "samtools view -c " + bo1_s + "_allreads.cf.*.q*.s.rd.uniq.bam").strip()
            if uniq_all_filtered_bam:
                uniq_all_filtered_bam = int(uniq_all_filtered_bam)

            if merged_filtered_bam == 0:
                percent_m_dup = "NaN"
            else:
                percent_m_dup = str(float(rd_merged_filtered_bam) / float(merged_filtered_bam))
            if merged == 0:
                percent_m_endog = "NaN"
            else:
                percent_m_endog = str((100.0 * float(merged_filtered_bam)) / float(merged))
            if reads_debarcoded == 0:
                percent_all_endog = "NaN"
            else:
                percent_all_endog = str((100 * float(all_filtered_bam)) / float(rawreads))

            chkpipe = subprocess.Popen(
                ['/bin/bash', '-c', "samtools view " + bo1_s + "_allreads.cf.*.q*.s.rd.bam "], stdout=PIPE)
            chkout = chkpipe.communicate()[0]
            if chkout:
                avglen_all_mapped = bash_command(
                    "samtools view " + bo1_s + "_allreads.cf.*.q*.s.rd.bam | awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'").strip()
            else:
                avglen_all_mapped = "NaN"

            bwa_q_pattern = out_sample + ".M.cf." + bwa_ref + ".q*.s.bam"
            files = os.listdir(bwa_output1)
            bwa_q = "N/A"
            for name in files:
                if (fnmatch.fnmatch(name, bwa_q_pattern)):
                    bwaqname = name
                    bwaqcols = bwaqname.split(".")
                    bwa_q = bwaqcols[4].strip()

            # MapDamage summary
            firstpos = "N/A"
            ct_damagefile = wd + "/mapDamage/mapDamage_" + out_sample + "_" + bwa_ref + "/5pCtoT_freq.txt"
            ctc = []
            if os.path.isfile(ct_damagefile):
                ctin = open(ct_damagefile, 'r')
                for ctinline in ctin:
                    ctc.append(ctinline)
                firstpos = ctc[1].split()[1].strip()

            else:
                print "No file matching " + ct_damagefile + " found."

            ga_damagefile = wd + "/mapDamage/mapDamage_" + out_sample + "_" + bwa_ref + "/3pGtoA_freq.txt"
            gac = []
            lastpos = "N/A"
            if os.path.isfile(ga_damagefile):
                gain = open(ga_damagefile, 'r')
                for gainline in gain:
                    gac.append(gainline)
                lastpos = gac[1].split()[1].strip()
            else:
                print "No file matching " + ga_damagefile + " found."

            percent_firstpos = bash_command("echo \"scale=2; 100 * " + str(firstpos) + "/1\" | bc").strip()
            percent_lastpos = bash_command("echo \"scale=2; 100 * " + str(lastpos) + "/1\" | bc").strip()

            # 1
            sumfile.write(str(today))
            sumfile.write(str("	"))
            # 2
            sumfile.write(str(in_sample))
            sumfile.write(str("	"))
            # 3
            sumfile.write(str(lib))
            sumfile.write(str("	"))
            # 4
            sumfile.write(str(out_sample))
            sumfile.write(str("	"))
            # 5
            sumfile.write(str(rawreads))
            sumfile.write(str("	"))
            # 6
            sumfile.write(str(reads_debarcoded))
            sumfile.write(str("	"))
            # 7
            sumfile.write(str(percentbarcodes))
            sumfile.write(str("	"))
            # 8
            sumfile.write(str(reads_debarcoded))
            sumfile.write(str("	"))
            # 9
            sumfile.write(str(merged))
            sumfile.write(str("	"))
            # 10
            sumfile.write(str(wadapters))
            sumfile.write(str("	"))
            # 11
            sumfile.write(str(discarded))
            sumfile.write(str("	"))
            # 12
            sumfile.write(str(merged))
            sumfile.write(str("	"))
            # 13
            sumfile.write(str(merged_filtered_bam))
            sumfile.write(str("	"))
            # 14
            sumfile.write(str(rd_merged_filtered_bam))
            sumfile.write(str("	"))
            # 15
            sumfile.write(str(uniq_all_filtered_bam))
            sumfile.write(str("	"))
            # 16
            sumfile.write(str(avglen_all_mapped))
            sumfile.write(str("	"))
            # 17
            sumfile.write(str(bwa_q))
            sumfile.write(str("	"))

            # 18
            sumfile.write(str(percentbarcodes))
            sumfile.write(str("	"))
            # 19
            sumfile.write(str(percent_m_dup))
            sumfile.write(str("	"))
            # 20
            sumfile.write(str(percent_m_endog))
            sumfile.write(str("	"))
            # 21
            sumfile.write(str(percentmerged))
            sumfile.write(str("	"))

            # 22
            sumfile.write(str(percent_firstpos))
            sumfile.write(str("	"))
            # 23
            sumfile.write(str(percent_lastpos))
            sumfile.write(str("	"))

            # 24
            sumfile.write(str("	"))

            # 25
            sumfile.write(str(cwd))
            sumfile.write(str("	"))

            sumfile.write("\n")
        sumfile.close()

    # Original sumfile
    sumfilename = wd + "/summary." + str(today) + ".txt"
    print "Summary to: ", sumfilename
    sumfile = open(sumfilename, 'w')
    # 1
    sumfile.write("Date   ")
    sumfile.write("	")
    # 2
    sumfile.write(" SeqID   ")
    sumfile.write("	")
    # 3
    sumfile.write(" Library ID   ")
    sumfile.write("	")
    # 4
    sumfile.write(" Sample   ")
    sumfile.write("	")
    # 5
    sumfile.write(" total reads ")
    sumfile.write("	")
    # 6
    sumfile.write(" total after BC removal   ")
    sumfile.write("	")
    # 7
    sumfile.write(" % reads with BC   ")
    sumfile.write("	")
    # 8
    sumfile.write(" pairs processed   ")
    sumfile.write("	")
    # 9
    sumfile.write(" pairs merged   ")
    sumfile.write("	")
    # 10
    sumfile.write(" pairs with adapters   ")
    sumfile.write("	")
    # 11
    sumfile.write(" pairs discarded   ")
    sumfile.write("	")
    # 12
    sumfile.write(" M reads used   ")
    sumfile.write("	")
    # 13
    sumfile.write(" flagstats M reads mapped  ")
    sumfile.write("	")
    # 14
    sumfile.write(" M reads mapped, rmdup   ")
    sumfile.write("	")
    # 15
    sumfile.write(" ALL reads, rmdup, uniqly mapping   ")
    sumfile.write("	")
    # 16
    sumfile.write(" Avg length all mapped reads   ")
    sumfile.write("	")
    # 17
    sumfile.write(" Quality used   ")
    sumfile.write("	")
    # 18
    sumfile.write(" Ref used   ")
    sumfile.write("	")
    # 19
    sumfile.write(" flagstats M reads mapped   ")
    sumfile.write("	")
    # 20
    sumfile.write("M reads mapped, rmdup   ")
    sumfile.write("	")
    # 21
    sumfile.write(" M reads, mapped, rmdup, uniq   ")
    sumfile.write("	")
    # 22
    sumfile.write(" Avg length all mapped mito reads   ")
    sumfile.write("	")
    # 23
    sumfile.write(" Q used   ")
    sumfile.write("	")
    # 24
    sumfile.write(" Ref used   ")
    sumfile.write("	")
    # 25
    sumfile.write("reads mapped   ")
    sumfile.write("	")
    # 26
    sumfile.write("avg coverage  ")
    sumfile.write("	")
    # 27
    sumfile.write("haplotype  ")
    sumfile.write("	")
    # 28
    sumfile.write("% reads with barcodes (no mismatches) ")
    sumfile.write("	")
    # 29
    sumfile.write("% not duplicate ")
    sumfile.write("	")
    # 30
    sumfile.write("q20 %endogenous ")
    sumfile.write("	")
    # 31
    sumfile.write("% merged ")
    sumfile.write("	")
    # 32
    sumfile.write("Average M fragment length ")
    sumfile.write("	")
    # 33
    sumfile.write("Sex estimate ")
    sumfile.write("	")
    # 34
    sumfile.write("5' damage ")
    sumfile.write("	")
    # 35
    sumfile.write("3' damage ")
    sumfile.write("	")
    # 36
    sumfile.write("analysis notes ")
    sumfile.write("	")
    # 37
    sumfile.write("location")
    sumfile.write("	")
    # 38
    sumfile.write("RX Sex assignment ")
    sumfile.write("	")
    # 39
    sumfile.write("RX estimate ")
    sumfile.write("	")
    # 40
    sumfile.write("RX conf interval ")
    sumfile.write("	")
    # 41
    sumfile.write("Strongly diagnostic positions ")
    sumfile.write("	")
    # 42
    sumfile.write("Frags Orig. ")
    sumfile.write("	")
    # 43
    sumfile.write("Frags Contam. ")
    sumfile.write("	")
    # 44
    sumfile.write("Frags Total ")
    sumfile.write("	")
    # 45
    sumfile.write("Frags Contam. Rate ")
    sumfile.write("	")
    # 46
    sumfile.write("All reads PMD filtered ")
    sumfile.write("	")
    # 46
    sumfile.write("PMD threshold ")
    sumfile.write("	")

    sumfile.write("\n")

    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        in_sample = bccols[0]
        out_sample = bccols[1]
        output = wd + "/" + out_sample

        # Get lib name if its there
        in_samplecols = in_sample.split("-")
        if len(in_samplecols) >= 2:
            lib = in_samplecols[1]
        else:
            lib = "???"

        if seqprep_output_in_output:
            seqprep_output = output + "/" + seqprep_output_orig
        else:
            seqprep_output = seqprep_output_orig
        so_s = seqprep_output + "/" + in_sample

        rawgzpattern = in_sample + "*_L00*_R1_001.fastq.gz"
        files = os.listdir(raw)
        r1name = None
        for name in files:
            if (fnmatch.fnmatch(name, rawgzpattern)):
                r1name = raw + "/" + name
                r1outname = os.path.splitext(r1name)[0]
                if not os.path.isfile(r1outname):
                    with gzip.open(r1name, 'rb') as f_in, open(r1outname, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

        # Will take the FIRST matching file
        rawfile = None
        rawpattern = in_sample + "*_L00*_R1_001.fastq"
        files = os.listdir(raw)
        for name in files:
            if (fnmatch.fnmatch(name, rawpattern)):
                rawfile = raw + "/" + name
                break
        if not os.path.isfile(rawfile):
            print "No file matching " + rawpattern + " found. Exiting."
            exit(1)

        rawin = open(rawfile, 'r')
        rawc = []
        for rawcinline in rawin:
            rawc.append(rawcinline)

        rawclength = len(rawc)
        rawreads = str(float(rawclength) / 4.0)

        if leehom:
            lh_file = seqprep_output + "/LH." + in_sample + ".stderr.txt"
            lhc = []
            if os.path.isfile(lh_file):
                lhin = open(lh_file, 'r')
                for lhinline in lhin:
                    lhc.append(lhinline)

            else:
                print "No file matching " + lh_file + " found. Exiting."
                exit(1)

            rdbname = wobc_output + "/" + in_sample + ".woBC_R1.fastq"
            if not os.path.isfile(rdbname):
                with gzip.open(rdbname + ".gz", 'rb') as f_in, open(rdbname, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            reads_debarcoded = int(bash_command("wc -l " + rdbname).split()[0]) / 4

            mname = so_s + ".M.fq"
            if not os.path.isfile(mname):
                with gzip.open(mname + ".gz", 'rb') as f_in, open(mname, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            merged = int(bash_command("wc -l " + mname).split()[0]) / 4
            noperc = lhc[1].split("\t")[0]
            splitspc = noperc.split(" ")
            wadapters = splitspc[2]
            discarded = "N/A"

        else:

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
            merged = int(spc[2].split()[2].strip())
            wadapters = int(spc[3].split()[3].strip())
            discarded = int(spc[4].split()[2].strip())

        rdbf = float(reads_debarcoded)
        if rdbf == 0:
            percentmerged = "N/A"
            print "Reads_debarcoded is 0. Are you sure you have the right indexes?"
        else:
            percentmerged = str((100.0 * float(merged) / rdbf))

        rrf = float(rawreads)
        if rrf == 0:
            percentbarcodes = "N/A"
        else:
            percentbarcodes = str((100.0 * rdbf) / rrf)

        bwa_output1 = output + "/BWA_" + bwa_ref
        bo1_s = bwa_output1 + "/" + out_sample

        mpattern = out_sample + ".M.cf.*.q*.s.bam"
        files = os.listdir(bwa_output1)
        mname = None
        for name in files:
            if (fnmatch.fnmatch(name, mpattern)):
                mname = name

        if not mname:
            print "ERROR: Cannot find match to" + mpattern + " Exiting. "
            exit(1)

        allpattern = out_sample + "_allreads.cf.*.q*.s.bam"
        files = os.listdir(bwa_output1)
        allname = None
        for name in files:
            if (fnmatch.fnmatch(name, allpattern)):
                allname = name
        if not allname:
            print "ERROR: Cannot find " + allpattern + " Exiting. "
            exit(1)

        merged_filtered_bam = bash_command("samtools view -c " + bo1_s + ".M.cf.*.q*.s.bam").strip()
        if merged_filtered_bam:
            merged_filtered_bam = int(merged_filtered_bam)
        rd_merged_filtered_bam = bash_command("samtools view -c " + bo1_s + ".M.cf.*.q*.s.rd.bam").strip()
        if rd_merged_filtered_bam:
            rd_merged_filtered_bam = int(rd_merged_filtered_bam)
        all_filtered_bam = bash_command("samtools view -c " + bo1_s + "_allreads.cf.*.q*.s.bam").strip()
        if all_filtered_bam:
            all_filtered_bam = int(all_filtered_bam)
        rd_all_filtered_bam = bash_command("samtools view -c " + bo1_s + "_allreads.cf.*.q*.s.rd.bam").strip()
        if rd_all_filtered_bam:
            rd_all_filtered_bam = int(rd_all_filtered_bam)

        uniq_all_filtered_bam = bash_command("samtools view -c " + bo1_s + "_allreads.cf.*.q*.s.rd.uniq.bam").strip()
        if uniq_all_filtered_bam:
            uniq_all_filtered_bam = int(uniq_all_filtered_bam)

        if merged_filtered_bam == 0:
            percent_m_dup = "NaN"
        else:
            percent_m_dup = str(float(rd_merged_filtered_bam) / float(merged_filtered_bam))
        if merged == 0:
            percent_m_endog = "NaN"
        else:
            percent_m_endog = str((100.0 * float(merged_filtered_bam)) / float(merged))
        if reads_debarcoded == 0:
            percent_all_endog = "NaN"
        else:
            percent_all_endog = str((100 * float(all_filtered_bam)) / float(rawreads))



        chkpipe = subprocess.Popen(
            ['/bin/bash', '-c', "samtools view " + bo1_s + "_allreads.cf.*.q*.s.rd.bam "], stdout=PIPE)
        chkout = chkpipe.communicate()[0]
        if chkout:
            avglen_all_mapped = bash_command(
                "samtools view " + bo1_s + "_allreads.cf.*.q*.s.rd.bam | awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'").strip()
        else:
            avglen_all_mapped = "NaN"
        # summary BWA (mito) statistics
        bwa_output2 = output + "/BWA_rCRS"
        bo2_s = bwa_output2 + "/" + out_sample

        mpattern = out_sample + ".M.cf.*.q*.s.bam"
        files = os.listdir(bwa_output2)
        mname = None
        for name in files:
            if (fnmatch.fnmatch(name, mpattern)):
                mname = name

        if not mname:
            print "ERROR: Cannot find match to" + mpattern + " Exiting. "
            exit(1)

        allpattern = out_sample + "_allreads.cf.*.q*.s.bam"
        files = os.listdir(bwa_output2)
        allname = None
        for name in files:
            if (fnmatch.fnmatch(name, allpattern)):
                allname = name
        if not allname:
            print "ERROR: Cannot find " + allpattern + " Exiting. "
            exit(1)

        m_filtered_bam_mito = bash_command("samtools view -c " + bo2_s + ".M.cf.*.q*.s.bam").strip()
        rd_m_filtered_bam_mito = bash_command("samtools view -c " + bo2_s + ".M.cf.*.q*.s.rd.bam").strip()
        all_filtered_bam_mito = bash_command("samtools view -c " + bo2_s + "_allreads.cf.*.q*.s.bam").strip()
        rd_filtered_bam_mito = bash_command("samtools view -c " + bo2_s + "_allreads.cf.*.q*.s.rd.bam").strip()
        all_uniq_mito = bash_command("samtools view -c " + bo2_s + "_allreads.cf.*.q*.s.rd.uniq.bam").strip()

        chkpipe = subprocess.Popen(
            ['/bin/bash', '-c', "samtools view " + bo2_s + "_allreads.cf.*.q*.s.rd.uniq.bam"], stdout=PIPE)
        chkout = chkpipe.communicate()[0]
        if chkout:
            avglen_all_mito = bash_command(
                "samtools view " + bo2_s + "_allreads.cf.*.q*.s.rd.uniq.bam | awk '{SUM+=length($10);DIV++}END{print SUM/DIV}'").strip()
        else:
            avglen_all_mito = "NaN"

        bwa_q_pattern = out_sample + ".M.cf." + bwa_ref + ".q*.s.bam"
        files = os.listdir(bwa_output1)
        bwa_q = "N/A"
        for name in files:
            if (fnmatch.fnmatch(name, bwa_q_pattern)):
                bwaqname = name
                bwaqcols = bwaqname.split(".")
                bwa_q = bwaqcols[4].strip()

        mitopattern = out_sample + ".M.cf." + mito_ref + ".q*.s.bam"
        files = os.listdir(bwa_output2)
        mito_q = "N/A"
        for name in files:
            if (fnmatch.fnmatch(name, mitopattern)):
                mitoname = name
                mitocols = mitoname.split(".")
                mito_q = mitocols[4].strip()
        if not nomia:
            # MIA summary
            mia_output = output + "/MIA_output"
            # Get FIRST matching MIA output file

            miapattern = out_sample + "*.mia_stats.txt"
            files = os.listdir(mia_output)
            mianame = None
            for name in files:
                if (fnmatch.fnmatch(name, miapattern)):
                    mianame = name
                    break

            mia_output_file = mia_output + "/" + mianame

            miac = []
            if os.path.isfile(mia_output_file):
                miain = open(mia_output_file, 'r')
                for mialine in miain:
                    miac.append(mialine)

            else:
                print "No file matching " + mia_output_file + " found. Exiting."
                exit(1)
            if len(miac) >= 5:  # crashing Mia might leave an empty file
                mia_reads = miac[2].split()[7]
                mia_coverage = miac[4].split()[3]
            else:
                mia_reads = "N/A"
                mia_coverage = "N/A"
        else:
            mia_reads = "N/A"
            mia_coverage = "N/A"

        # MapDamage summary
        ct_damagefile = wd + "/mapDamage/mapDamage_" + out_sample + "_" + bwa_ref + "/5pCtoT_freq.txt"
        ctc = []
        firstpos = "N/A"
        if os.path.isfile(ct_damagefile):
            ctin = open(ct_damagefile, 'r')
            for ctinline in ctin:
                ctc.append(ctinline)
            firstpos = ctc[1].split()[1].strip()

        else:
            print "No file matching " + ct_damagefile + " found."

        ga_damagefile = wd + "/mapDamage/mapDamage_" + out_sample + "_" + bwa_ref + "/3pGtoA_freq.txt"
        gac = []
        lastpos = "N/A"
        if os.path.isfile(ga_damagefile):
            gain = open(ga_damagefile, 'r')
            for gainline in gain:
                gac.append(gainline)
            lastpos = gac[1].split()[1].strip()
        else:
            print "No file matching " + ga_damagefile + " found."

        percent_firstpos = bash_command("echo \"scale=2; 100 * " + str(firstpos) + "/1\" | bc").strip()
        percent_lastpos = bash_command("echo \"scale=2; 100 * " + str(lastpos) + "/1\" | bc").strip()

        sexpattern = out_sample + "*.sex.txt"
        files = os.listdir(output)
        sexname = None
        for name in files:
            if (fnmatch.fnmatch(name, sexpattern)):
                sexname = name
                break

        sexc = []
        if sexname is not None:
            sex_output_file = output + "/" + sexname
            if os.path.isfile(sex_output_file):
                sexin = open(sex_output_file, 'r')
                for sexline in sexin:
                    sexc.append(sexline)

        if len(sexc) >= 2:
            ry_sex = str(sexc[1].split("\t")[6]).strip()
        else:
            ry_sex = "N/A"

        rxfilename = output + "/" + out_sample + ".Rx.txt"
        rxc = []
        if os.path.isfile(rxfilename):
            rxin = open(rxfilename, 'r')
            for rxline in rxin:
                rxc.append(rxline)

        if len(rxc) >= 20:
            rxsex = str(rxc[19].split(":")[1]).strip()
        else:
            rxsex = "N/A"

        if len(rxc) >= 21:
            rxsex_ci = str(rxc[20].split(":")[1]).strip()
            rxsex_ci.replace(" ", "-")
        else:
            rxsex_ci = "N/A"

        if len(rxc) >= 22:
            rxsex_assignment = str(rxc[21].split(":")[1]).strip()
        else:
            rxsex_assignment = "N/A"

        # 2a_contam
        contamfilename = wd + "/ccheck_results/results.raw.contamination.txt"
        cxc = []

        ct_sdpos = "N/A"
        ct_forig = "N/A"
        ct_fcontam = "N/A"
        ct_ftotal = "N/A"
        ct_ct = "N/A"

        if os.path.isfile(contamfilename):
            ctin = open(contamfilename, 'r')
            for ctline in ctin:
                cxc.append(ctline)

        for cx in cxc:
            cxcols = cx.split("\t")
            if cxcols[0] == out_sample:
                ct_sdpos = cxcols[1].strip()
                ct_forig = cxcols[2].strip()
                ct_fcontam = cxcols[3].strip()
                ct_ftotal = cxcols[4].strip()
                ct_ct = cxcols[5].strip()


        # 4 PMDTools
        pmd_filtered_bam = bash_command(
            "samtools view -c " + bo2_s + "_allreads.cf.*.q*.pmds" + pmd_threshold + "filter.bam").strip()
        if pmd_filtered_bam:
            pmd_filtered_bam = int(rd_all_filtered_bam)

        # 1
        sumfile.write(str(today))
        sumfile.write(str("	"))
        # 2
        sumfile.write(str(in_sample))
        sumfile.write(str("	"))
        # 3
        sumfile.write(str(lib))
        sumfile.write(str("	"))
        # 4
        sumfile.write(str(out_sample))
        sumfile.write(str("	"))
        # 5
        sumfile.write(str(rawreads))
        sumfile.write(str("	"))
        # 6
        sumfile.write(str(reads_debarcoded))
        sumfile.write(str("	"))
        # 7
        sumfile.write(str(percentbarcodes))
        sumfile.write(str("	"))
        # 8
        sumfile.write(str(reads_debarcoded))
        sumfile.write(str("	"))
        # 9
        sumfile.write(str(merged))
        sumfile.write(str("	"))
        # 10
        sumfile.write(str(wadapters))
        sumfile.write(str("	"))
        # 11
        sumfile.write(str(discarded))
        sumfile.write(str("	"))
        # 12
        sumfile.write(str(merged))
        sumfile.write(str("	"))
        # 13
        sumfile.write(str(merged_filtered_bam))
        sumfile.write(str("	"))
        # 14
        sumfile.write(str(rd_merged_filtered_bam))
        sumfile.write(str("	"))
        # 15
        sumfile.write(str(uniq_all_filtered_bam))
        sumfile.write(str("	"))
        # 16
        sumfile.write(str(avglen_all_mapped))
        sumfile.write(str("	"))
        # 17
        sumfile.write(str(bwa_q))
        sumfile.write(str("	"))
        # 18
        sumfile.write(str(bwa_ref))
        sumfile.write(str("	"))
        # 19
        sumfile.write(str(m_filtered_bam_mito))
        sumfile.write(str("	"))
        # 20
        sumfile.write(str(rd_m_filtered_bam_mito))
        sumfile.write(str("	"))
        # 21
        sumfile.write(str(all_uniq_mito))
        sumfile.write(str("	"))
        # 22
        sumfile.write(str(avglen_all_mito))
        sumfile.write(str("	"))
        # 23
        sumfile.write(str(mito_q))
        sumfile.write(str("	"))
        # 24
        sumfile.write(str(mito_ref))
        sumfile.write(str("	"))
        # 25
        sumfile.write(str(mia_reads))
        sumfile.write(str("	"))
        # 26
        sumfile.write(str(mia_coverage))
        sumfile.write(str("	"))
        # 27
        sumfile.write(str("	"))
        # 28
        sumfile.write(str(percentbarcodes))
        sumfile.write(str("	"))
        # 29
        sumfile.write(str(percent_m_dup))
        sumfile.write(str("	"))
        # 30
        sumfile.write(str(percent_m_endog))
        sumfile.write(str("	"))
        # 31
        sumfile.write(str(percentmerged))
        sumfile.write(str("	"))

        # 32
        sumfile.write(str("	"))

        # 33
        sumfile.write(str(ry_sex))
        sumfile.write(str("	"))
        # 34
        sumfile.write(str(percent_firstpos))
        sumfile.write(str("	"))
        # 35
        sumfile.write(str(percent_lastpos))
        sumfile.write(str("	"))

        # 36
        sumfile.write(str("	"))

        # 37
        sumfile.write(str(cwd))
        sumfile.write(str("	"))
        # 38
        sumfile.write(str(rxsex_assignment))
        sumfile.write(str("	"))
        # 39
        sumfile.write(str(rxsex))
        sumfile.write(str("	"))
        # 40
        sumfile.write(str(rxsex_ci))
        sumfile.write(str("	"))
        # 41
        sumfile.write(str(ct_sdpos))
        sumfile.write(str("	"))
        # 42
        sumfile.write(str(ct_forig))
        sumfile.write(str("	"))
        # 43
        sumfile.write(str(ct_fcontam))
        sumfile.write(str("	"))
        # 44
        sumfile.write(str(ct_ftotal))
        sumfile.write(str("	"))
        # 45
        sumfile.write(str(ct_ct))
        sumfile.write(str("	"))
        # 46
        sumfile.write(str(pmd_filtered_bam))
        sumfile.write(str("	"))
        # 47
        sumfile.write(str(pmd_threshold))
        sumfile.write(str("	"))

        sumfile.write("\n")
    sumfile.close()
    cmdfile.close()
    logfile.close()
    print "results.py complete."
    exit(0)

