#!/usr/bin/python

#####################################
#####        HPG Lab            #####
#####    updated April 2018     #####
#####       MJJ                 #####
#####################################

# Converted to Python, extended and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# Script for processing shotgun sequencing data on HPG Lab
# Based off of UCSC's Pete Heintzman processing script
# Sequence data should already be trimmed by barcode remover and SeqPrep2 in SeqPrep2 output
# Requires 'barcode' file with fastq prefix, output name and internal barcodes

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


def bwa_align_se(refname, refloc, bso_s, btail, bbwa_output, bsample, bq):
    infile = bso_s + btail + ".fastq"
    bwaalignline = "bwa aln -l " + seed_disable + " -n " + bwamaxedit + " -t " + threads + " " + refloc + " " + infile
    bwasline = "bwa samse " + refloc + " - " + infile
    samviewline = "samtools view -q " + bq + " -F 4 -buh - "
    samsortline = "samtools sort -@ " + threads + " -m 4G -o " + bbwa_output + "/" + bsample + btail + "." + refname + ".q" + bq + ".s.bam"
    seline = bwaalignline + " | " + bwasline + " | " + samviewline + " | " + samsortline
    bash_command(seline)


def bwa_align_pe(refname, refloc, bso_s, btail, bbwa_output, bsample, btail1, btail2, bq):
    infile1 = bso_s + btail1 + ".fastq"
    infile2 = bso_s + btail2 + ".fastq"
    bwasaline = "bwa sampe " + refloc + " <(bwa aln -l " + seed_disable + " -n " + bwamaxedit + " -t " + threads \
                + " " + refloc + " " + infile1 + ") " + "<(bwa aln -l " + seed_disable + " -n 0.01 -t " \
                + threads + " " + refloc + " " + infile2 + ") " + infile1 + " " + infile2
    samviewline = "samtools view -q " + bq + " -F 4 -buh - "
    samsortline = "samtools sort -@ " + threads + " -m 4G -o " + bbwa_output + "/" + bsample + btail + "." + refname + ".q" + bq + ".s.bam"
    peline = bwasaline + " | " + samviewline + " | " + samsortline
    bash_command(peline)


if __name__ == "__main__":

    print "\n****************\nMAP\n****************\n"

    parser = argparse.ArgumentParser(description="# This script:\n"
                                                 "1. uncompresses fastq.gz files\n"
                                                 "2. bwa maps to reference genome(s)\n"
                                                 "3. samtools processing (eg sam --> bam, rmdup)\n"
                                                 "4. Calculates frag length mean and distribution; outputs figure\n"
                                                 "5. MapDamage2 profiles\n"
                                                 "6. MIA maps to mtDNA (rCRS)\n"
                                                 "7. MALT analysis\n"
                                                 "8. Sex estimation > $SAMPLE.sex.txt.\n\t"
                                                 "- ", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.',
                        required=True)
    parser.add_argument('-wd', metavar='<wd>', help='Working directory. Defaults to current.', default='.')
    parser.add_argument('-mia_ref', metavar='<mia_ref>', help='mia_ref',
                        default='rCRS')
    parser.add_argument('-mtdna', metavar='<mtdna>', help='mtDNA mapping',
                        default='/data/genomes/rCRS.fas')
    parser.add_argument('-index_algorithm', metavar='<index_algorithm>',
                        help='If reference is <2Gb use is, if >2Gb use bwtsw',
                        default='is')
    parser.add_argument('-seed_disable', metavar='<seed_disable>',
                        help='Following ancient DNA data processing protocols',
                        default="1024")
    parser.add_argument('-threads', metavar='<threads>',
                        help='Number of concerrent threads to use when invoking software.',
                        default="23")
    parser.add_argument('-q', metavar='<q>', help='BWA min quality. 20 provides a fairly low cutoff',
                        default="20")
    parser.add_argument('-adna_matrix', metavar='<adna_matrix>', help='aDNA matrix',
                        default='/data/scripts/ancient.submat.txt')
    parser.add_argument('-max_misincorp_frequency', metavar='<max_misincorp_frequency>',
                        help=' Use 0.3 if not too badly damaged, use 0.5 if badly damaged',
                        default="0.3")
    parser.add_argument('-read_plot_length', metavar='<read_plot_length>',
                        help='The number of nucleotides to plot at the 5\' and 3\' ends of the read',
                        default="25")
    parser.add_argument('-max_read_length', metavar='<max_read_length>', help='The maximum read length to consider',
                        default="150")
    parser.add_argument('-prinseq_lite', metavar='<prinseq_lite>', help='prinseq_lite',
                        default='/data/scripts/prinseq-lite.pl')
    parser.add_argument('-prinseq_graphs', metavar='<prinseq_graphs>', help='prinseq_graphs',
                        default='/data/scripts/prinseq-graphs.pl')
    parser.add_argument('-combine_paired_end_reads', metavar='<combine_paired_end_reads>',
                        help='combine_paired_end_reads',
                        default='/data/scripts/combinePairedEndReads.pl')
    parser.add_argument('-split_paired_end_reads', metavar='<split_paired_end_reads>', help='split_paired_end_reads',
                        default='/data/scripts/splitPairedEndReads.pl')
    parser.add_argument('-frag_length_r', metavar='<frag_length_r>', help='frag_length_r',
                        default='/data/scripts/frag_len_hist.R')
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                        default='/data/adaptertrimmed')
    parser.add_argument('-bwaindex', dest='bwaindex', help='Need to index if never used the reference genome before.',
                        action='store_true')
    parser.set_defaults(bwaindex=False)
    parser.add_argument('-verbose', dest='verbose', help='Print stdout and stderr to console.',
                        action='store_true')
    parser.set_defaults(verbose=False)
    parser.add_argument('-nomia', dest='nomia', help='Do not run MIA.',
                        action='store_true')
    parser.set_defaults(nomia=False)
    parser.add_argument('-malt', dest='malt', help='Run MALT.',
                        action='store_true')
    parser.set_defaults(malt=False)
    parser.add_argument('-nosexest', dest='nosexest', help='Do not run sex estimation.',
                        action='store_true')
    parser.set_defaults(nosexest=False)
    parser.add_argument('-kraken', dest='kraken', help='Release the kraken!',
                        action='store_true')
    parser.set_defaults(kraken=False)
    parser.add_argument('-overwrite', dest='overwrite', help='Overwrite existing files and directories.',
                        action='store_true')
    parser.set_defaults(overwrite=False)
    parser.add_argument('-bformat', dest='bformat', help='Barcode format type old (barcodes as sequence not numbers)?.',
                        action='store_true')
    parser.set_defaults(bformat=False)
    parser.add_argument('-refs', dest='refs', nargs='+', default=[],
                        help='List of reference sequences other than hg19 and rCRS.')
    parser.add_argument('-sexref', metavar='<sexref>', help="Reference sequence for Skoglunds sex estimation scripts",
                        default='/data/genomes/hg19.fa')
    parser.add_argument('-bwamaxedit', metavar='<bwamaxedit>',
                        help='Maximum edit distance if the value is INT, or the fraction of missing alignments given 2 percent uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths.',
                        default="0.01")
    parser.add_argument('-clipleft', metavar='<clipleft>',
                        help='Number of bases to clip from left end of BAM file after mapDamage run.',
                        default="0")
    parser.add_argument('-clipright', metavar='<clipright>',
                        help='Number of bases to clip from right end of BAM file after mapDamage run.',
                        default="0")
    parser.add_argument('-maltdb', metavar='<maltdb>', help='MALT database.',
                        default="/data/malt/nr/")
    parser.add_argument('-maltblast', metavar='<maltblast>', help='MALT BLAST type.',
                        default="BlastX")
    parser.add_argument('-skipprinseq', dest='skipprinseq', help='Skip the prinseq part of script.',
                        action='store_true')
    parser.set_defaults(overwrite=False)

    args = parser.parse_args()
    wd = args.wd
    bcfile = args.bc_file
    mia_ref = args.mia_ref
    mtdna = args.mtdna
    index_algorithm = args.index_algorithm
    seed_disable = args.seed_disable
    threads = args.threads
    q = args.q
    adna_matrix = args.adna_matrix
    max_misincorp_frequency = args.max_misincorp_frequency
    read_plot_length = args.read_plot_length
    max_read_length = args.max_read_length
    bformat = bool(args.bformat)

    prinseq_lite = args.prinseq_lite
    prinseq_graphs = args.prinseq_graphs
    skipprinseq = bool(args.skipprinseq)
    combine_paired_end_reads = args.combine_paired_end_reads
    split_paired_end_reads = args.split_paired_end_reads

    frag_length_r = args.frag_length_r
    seqprep_output = args.seqprep_output
    bwaindex = bool(args.bwaindex)
    nomia = bool(args.nomia)
    malt = bool(args.malt)
    nosexest = bool(args.nosexest)
    kraken = bool(args.kraken)
    overwrite = bool(args.overwrite)
    verbose = bool(args.verbose)
    refs = args.refs
    sexref = args.sexref
    bwamaxedit = args.bwamaxedit
    clipleft = args.clipleft
    clipright = args.clipright
    maltdb = args.maltdb
    maltblast = args.maltblast

    arefs = ["/data/genomes/hg19.fa", "/data/genomes/rCRS.fas"]
    for ref in refs:
        arefs.append(ref)

    os.chdir(wd)
    cwd = os.getcwd()
    print "Working in: ", cwd

    today = datetime.date.today()
    logfilename = wd + "/out.map." + str(today) + ".log"
    print "Logging to: ", logfilename
    logfile = open(logfilename, 'w')

    refdic = {}
    logfile.write("Reference sequences used: \n")
    for ref in arefs:
        refname = os.path.basename(ref)
        filebase, fileext = os.path.splitext(refname)
        refdic[filebase] = ref
        logfile.write(ref + "\n")
    logfile.write("-------------------\n")

    outfile = open("1_cmds", 'w')

    newdir = wd + "/Frag_lengths"
    if os.path.exists(newdir):
        pass
    else:
        os.mkdir(newdir)

    newdir = wd + "/mapDamage"
    if os.path.exists(newdir):
        pass
    else:
        os.mkdir(newdir)

    logfile.write("Parameters used: \n")
    logfile.write("Prinseq lite, lc_method = dust, threshold 7\n")
    logfile.write("BWA aln -l " + seed_disable + " -n " + bwamaxedit + " -t " + threads + "\n")
    logfile.write("Map quality cutoff = " + q + "\n")
    logfile.write("\n-------------------------------------------------\n")

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
    if bformat:
        bcposline = "bcpos	bcpos	TCGAACA	AGCACAT ATGTGCT TGTTCGA"
    else:
        bcposline = "bcpos	bcpos	198	205"
    nobcposline = "nobcpos	nobcpos		"
    if barcodespresent:
        bc.append(bcposline)
    else:
        bc.append(nobcposline)

    bclength = len(bc)
    print "Number of entries: ", bclength

    flist = []

    filestorezip = []

    print "\nUncompressing..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        in_sample = bccols[0]
        out_sample = bccols[1]
        output = wd + "/" + out_sample
        so_s = seqprep_output + "/" + in_sample
        flist.append(in_sample)

        logfile.write("Uncompressing " + in_sample + "\n")

        if not os.path.isfile(so_s + ".M.fq"):  # merged reads
            with gzip.open(so_s + ".M.fq.gz", 'rb') as f_in, open(so_s + ".M.fq", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        if not os.path.isfile(so_s + ".F.fq"):  # forward reads
            with gzip.open(so_s + ".F.fq.gz", 'rb') as f_in, open(so_s + ".F.fq", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        if not os.path.isfile(so_s + ".R.fq"):  # reverse reads
            with gzip.open(so_s + ".R.fq.gz", 'rb') as f_in, open(so_s + ".R.fq", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        filestorezip.append(so_s + ".M.fq")
        filestorezip.append(so_s + ".F.fq")
        filestorezip.append(so_s + ".R.fq")

    flength = len(flist)
    print "Number of entries: ", flength

    print "\nComplexity filtering using prinseq..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        in_sample = bccols[0]
        out_sample = bccols[1]
        output = wd + "/" + out_sample
        so_s = seqprep_output + "/" + in_sample

        logfile.write('Prinseq removing low complexity from merged ' + in_sample + "\n")
        bash_command(
            "perl " + prinseq_lite + " -fastq " + so_s + ".M.fq -out_good " + so_s + ".M.cf -out_bad null -lc_method dust -lc_threshold 7 -line_width 0")

        bash_command(
            "perl " + combine_paired_end_reads + " " + so_s + ".F.fq " + so_s + ".R.fq " + so_s + ".uM_combined.fastq")

        bash_command(
            "perl " + prinseq_lite + " -fastq " + so_s + ".uM_combined.fastq -out_good " + so_s + ".uM_combined.cf -out_bad null -lc_method dust -lc_threshold 7 -line_width 0")

        bash_command("perl " + split_paired_end_reads + " " + so_s + ".uM_combined.cf.fastq")

        if os.path.isfile(so_s + ".uM_combined.cf.fastq_1"):
            shutil.move(so_s + ".uM_combined.cf.fastq_1", so_s + ".F.cf.fastq")
            filestorezip.append(so_s + ".F.cf.fastq")
        if os.path.isfile(so_s + ".uM_combined.cf.fastq_2"):
            shutil.move(so_s + ".uM_combined.cf.fastq_2", so_s + ".R.cf.fastq")
            filestorezip.append(so_s + ".R.cf.fastq")

        filestorezip.append(so_s + ".M.cf.fastq")

        logfile.write("Removing " + in_sample + "\n")
        sopattern = in_sample + ".uM_combined*.fastq"
        files = os.listdir(seqprep_output)
        soname = None
        for name in files:
            if (fnmatch.fnmatch(name, sopattern)):
                os.remove(seqprep_output + "/" + name)

    print "\nCreating BWA directories..."
    bar = progressbar.ProgressBar()
    for i in bar(range(bclength)):
        bcline = bc[i]
        bccols = bcline.split()
        in_sample = bccols[0]
        out_sample = bccols[1]
        output = wd + "/" + out_sample
        so_s = seqprep_output + "/" + in_sample
        for key, value in refdic.iteritems():
            newdir = output + "/BWA_" + key
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

    ########### BWA - aligning reads to a reference sequence ############
    if bwaindex:
        for key, value in refdic.iteritems():
            logfile.write("Indexing to " + key + "...\n")
            bash_command("bwa index -p " + value + " -a " + index_algorithm + " " + value)

    for key, value in refdic.iteritems():
        print "\nAligning complexity-filtered merged reads to " + key + "..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bwa_align_se(key, value, so_s, ".M.cf", bwa_output, out_sample, q)

    for key, value in refdic.iteritems():
        print "\nAligning complexity-filtered unmerged reads to " + key + "..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bwa_align_pe(key, value, so_s, ".uM.cf.pe", bwa_output, out_sample, ".F.cf", ".R.cf", q)

    for key, value in refdic.iteritems():
        print "\nCombining filtered merged and unmerged BAMs for " + key + "..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample

            logfile.write("-------------------------------------------------" + "\n")
            logfile.write("Combining filtered merged and unmerged BAMs for " + out_sample + "\n")
            bash_command("samtools merge -f " + bo_s + "_allreads.cf." + key + ".q" + q + ".bam " + bo_s + ".M.cf." \
                         + key + ".q" + q + ".s.bam " + bo_s \
                         + ".uM.cf.pe." + key + ".q" + q + ".s.bam")
            # Sort and index allreads BAM file
            bash_command(
                "samtools sort " + bo_s + "_allreads.cf." + key + ".q" + q + ".bam -o " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam -T " + out_sample + ".all")
            bash_command("samtools index " + bo_s + ".M.cf." + key + ".q" + q + ".s.bam")
            bash_command("samtools index " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam")

            # Remove duplicates and index, for all reads and merged reads using Dedup
            #       # All reads
            logfile.write("-------------------------------------------------" + "\n")
            logfile.write("Removing duplicates and indexing BAM for " + out_sample + " allreads" + "\n")

            bash_command(
                "java -jar $DEDUP -m -i " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam -o " + bwa_output + "/")
            shutil.move(bo_s + "_allreads.cf." + key + ".q" + q + ".s_rmdup.bam",
                        bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.bam")

            bash_command("samtools index " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.bam")

            # Merged Reads
            logfile.write("-------------------------------------------------" + "\n")
            logfile.write("Removing duplicates and indexing BAM for " + out_sample + " merged" + "\n")
            bash_command("java -jar $DEDUP -m -i " + bo_s + ".M.cf." + key + ".q" + q + ".s.bam -o " + bwa_output + "/")
            shutil.move(bo_s + ".M.cf." + key + ".q" + q + ".s_rmdup.bam",
                        bo_s + ".M.cf." + key + ".q" + q + ".s.rd.bam")
            bash_command("samtools index " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.bam")

            logfile.write("-------------------------------------------------" + "\n")
            logfile.write("Writing summary flagstat files for " + out_sample + "\n")

            # Get reads that are uniquely mapping
            samfline = "samtools view -h " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.bam"
            samfgreplines = " |grep -v 'XT:A:R'|grep -v 'XA:Z' |grep -v 'XT:A:M' |awk '{if($0~/X1:i:0/||$0~/^@/)print $0}' |"
            samgline = "samtools view -bS - > " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.uniq.bam"
            bash_command(samfline + samfgreplines + samgline)

            samfallline = "samtools view -h " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.bam"
            samgallline = "samtools view -bS - > " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.uniq.bam"
            bash_command(samfallline + samfgreplines + samgallline)

            # print average lengths of mapped reads using awk
            logfile.write(
                "Average length of uniquely mapping MERGED reads " + out_sample + " >> ./ avgmappedlen." + key + ".M.q" + q + ".txt" + "\n")
            bash_command(
                "samtools view " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.uniq.bam | awk '{SUM+=length($10);DIV++}END{print SUM/DIV}' >> ./avgmappedlen." + key + ".M.q" + q + ".txt" + "\n")

            # sorted
            bash_command(
                "samtools flagstat " + bo_s + ".M.cf." + key + ".q" + q + ".s.bam > " + bo_s + ".M.cf." + key + ".q" + q + ".s.flagstat.txt")
            bash_command(
                "samtools flagstat " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam > " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.flagstat.txt")

            # sorted, duplicates removed
            bash_command(
                "samtools flagstat " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.bam > " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.flagstat.txt")
            bash_command(
                "samtools flagstat " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.bam > " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.flagstat.txt")

            # sorted
            bash_command(
                "samtools idxstats " + bo_s + ".M.cf." + key + ".q" + q + ".s.bam > " + bo_s + ".M.cf." + key + ".q" + q + ".s.idxstats")
            bash_command(
                "samtools idxstats " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam > " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.idxstats")

            # sorted, duplicates removed
            bash_command(
                "samtools idxstats " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.bam > " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.idxstats")
            bash_command(
                "samtools idxstats " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.bam > " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.idxstats")

    # Compute fragment length distributions of all reads in the library (mapped)
    for key, value in refdic.iteritems():
        print "\nFragment length distributions from sorted bam file for " + key + "..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample
            frag_lengths = wd + "/Frag_lengths"

            bash_command(
                "samtools view " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam | awk '{print length($10)}'| sort -n | uniq -c | awk '{print $1\"\t\"$2}' | tail -n +2 > " + frag_lengths + "/" + out_sample + "_allreads.cf." + key + ".q" + q + ".s.fraglen.txt")
            rscrline = "Rscript " + frag_length_r + " " + frag_lengths + "/" + out_sample + "_allreads.cf." + key + ".q" + q + ".s.fraglen.txt " + frag_lengths + "/" + out_sample + "_allreads.cf." + key + ".q" + q + ".s.fraglen.pdf"
            bash_command(rscrline)

    for key, value in refdic.iteritems():
        print "\nMaking mapDamage plots for all quality filtered mapped reads " + key + "..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample

            newdir = wd + "/mapDamage/mapDamage_" + out_sample + "_" + key
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

            logfile.write("----------------------------------------------------" + "\n")
            logfile.write(
                "Making mapDamage plots for all quality filtered mapped reads " + out_sample + " " + key + "\n")

            bash_command("	mapDamage -i " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.bam -r " + value + " -l " \
                         + max_read_length + " --merge-reference-sequences --no-stats -d " + wd + "/mapDamage/mapDamage_" + out_sample \
                         + "_" + key + " -y " + max_misincorp_frequency + " -m " + read_plot_length + " -t " + key + "_" + out_sample)

            logfile.write("mapDamage plots finished " + out_sample + " " + key + "\n")
            logfile.write("----------------------------------------------------" + "\n")

    logfile.write("\n\n Clipping bases from BAM files with bamUtils\n")
    for key, value in refdic.iteritems():
        print "\nClipping ends with bamUtil... " + key + "..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):

            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample

            if os.path.isfile(bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.uniq.bam"):
                bash_command(
                    "bam trimBam " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.uniq.bam " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.uniq.clip.bam " + clipleft + " -R " + clipright)
                bash_command("samtools index " + bo_s + "_allreads.cf." + key + ".q" + q + ".s.rd.uniq.bam")
            if os.path.isfile(bo_s + ".M.cf." + key + ".q" + q + ".s.rd.uniq.bam"):
                bash_command(
                    "bam trimBam " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.uniq.bam " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.uniq.clip.bam " + clipleft + " -R " + clipright)
                bash_command("samtools index " + bo_s + ".M.cf." + key + ".q" + q + ".s.rd.uniq.bam")

    if not nomia:
        print "\nMEGAN and MIA -- Initial data processing..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample

            # Concatenate merged and unmerged reads - MEGAN and MIA
            logfile.write("----------------------------------------------------" + "\n")
            logfile.write("Prepping " + out_sample + " for MIA" + "\n")
            bash_command(
                "cat " + so_s + ".M.cf.fastq " + so_s + ".F.cf.fastq " + so_s + ".R.cf.fastq > " + so_s + ".all.SP.cf.fastq")
            filestorezip.append(so_s + ".all.SP.cf.fastq")

            # Convert FASTQ to FASTA - MEGAN and MIA (capture only)
            bash_command("fastq_to_fasta -n -Q33 -r -i " + so_s + ".all.SP.cf.fastq -o " + so_s + ".all.SP.cf.fasta")
            filestorezip.append(so_s + ".all.SP.cf.fasta")

            # Collapse identical reads - forward, reverse, and 5' duplicates - MEGAN and MIA (capture only)
            logfile.write("Collapsing identical reads with prinseq lite" + "\n")
            bash_command(
                "perl " + prinseq_lite + " -fasta " + so_s + ".all.SP.cf.fasta -out_good " + so_s + ".all.SP.cf.rd -out_bad null -derep 124 -line_width 0")
            filestorezip.append(so_s + ".all.SP.cf.rd.fasta")
            logfile.write("----------------------------------------------------" + "\n")

        print "\nStarting MIA with " + mia_ref + " and using merged/unmerged rmdup complexity filtered reads..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample
            mia_output = wd + "/" + out_sample + "/MIA_output"

            newdir = mia_output
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

            mia_s = mia_output + "/" + out_sample
            logfile.write(
                "Starting MIA with" + mia_ref + "and using merged/unmerged rmdup complexity filtered reads" + "\n")
            bash_command(
                "mia -r " + mtdna + " -f " + so_s + ".all.SP.cf.rd.fasta -c -C -U -s " + adna_matrix + " -i -F -k 14 -m  " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln")
            bash_command(
                "ma -M " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.* -f 3 > " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.mia_stats.txt")
            bash_command(
                "ma -M " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.* -f 2 > " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.mia_coverage_per_site.txt")
            bash_command(
                "ma -M " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.* -f 5 > " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.mia_consensus.fasta")
            bash_command(
                "ma -M " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.* -f 41 > " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.inputfornext.txt")
            bash_command(
                "perl /data/scripts/mia_consensus_coverage_filter.pl -c 3 -I < " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.inputfornext.txt > " + mia_s + ".all.SP.cf.rd." + mia_ref + ".maln.F.3xFiltered.consensus.fas")

            # Align mia_ref to MIA output

            # Call haplogroups using everything in mt and/or Y

    if malt:
        print "\nPerforming MALT analysis..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample
            malt_output = output + "/MALT_output"

            newdir = malt_output
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

            malt_s = malt_output + "/" + out_sample
            logfile.write("Performing MALT analysis " + out_sample + "\n")
            bash_command(
                "/data/malt/malt-run -i " + so_s + ".all.SP.cf.rd.fasta -d " + maltdb + " -m " + maltblast + " -o " + malt_output + " -t " + threads + " -L /data/install/MEGAN5-academic-license.txt")
            logfile.write("Done with MALT analysis" + "\n")
            logfile.write("----------------------------------------------------" + "\n")

    if not nosexest:
        print "\nPontus Skoglund's sex estimation..."

        sexreffullname = os.path.basename(sexref)
        sexrefname, fileext = os.path.splitext(sexreffullname)

        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key

            sexref_output = output + "/BWA_" + sexrefname
            bo_s = sexref_output + "/" + out_sample

            # Lets not have a divide by zero stop the groovy train here...
            chkpipe = subprocess.Popen(
                ['/bin/bash', '-c', "samtools view -q 30 " + bo_s + ".M.cf." + sexrefname + ".q" + q + ".s.rd.bam"],
                stdout=PIPE)
            chkout = chkpipe.communicate()[0]
            if chkout:
                logfile.write("Running sex estimation " + out_sample + "\n")
                logfile.write("Sex summary " + out_sample + " >> " + output + "/" + out_sample + ".sex.txt" + "\n")
                bash_command(
                    "samtools view -q 30 " + bo_s + ".M.cf." + sexrefname + ".q" + q + ".s.rd.bam | python /data/scripts/ry_compute.py >> " + output + "/" + out_sample + ".sex.txt")
            else:
                logfile.write(
                    "samtools view -q 30 " + bo_s + ".M.cf." + sexrefname + ".q" + q + ".s.rd.bam gives no output for ry_compute.py. No sex estimation peformed." + "\n")

            bash_command(
                "samtools idxstats " + bo_s + ".M.cf." + sexrefname + ".q" + q + ".s.rd.bam > " + bo_s + ".idxstats")

            bash_command(
                "Rscript /data/scripts/rx_identifier.r " + bo_s + " > " + output + "/" + out_sample + ".Rx.txt")
            if os.path.isfile("Rplots.pdf"):
                os.remove("Rplots.pdf")

    if kraken:
        print "\nDIAMOND and Kraken metagenome analysis, mapping fasta reads to NR NCBI database..."
        bar = progressbar.ProgressBar()
        for i in bar(range(bclength)):
            bcline = bc[i]
            bccols = bcline.split()
            in_sample = bccols[0]
            out_sample = bccols[1]
            output = wd + "/" + out_sample
            so_s = seqprep_output + "/" + in_sample
            bwa_output = wd + "/" + out_sample + "/BWA_" + key
            bo_s = bwa_output + "/" + out_sample

            dmnd_output = wd + "/DMND_output"
            do_s = dmnd_output + "/" + out_sample
            newdir = dmnd_output
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
            krkn_output = wd + "/Kraken_output"
            newdir = krkn_output
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
            ko_s = krkn_output + "/" + out_sample

            logfile.write("Performing DIAMOND analysis" + out_sample + "\n")
            bash_command(
                "diamond blastx -p " + threads + " -q " + so_s + ".all.SP.cf.rd.fasta -d /data/BLAST/diamond_nr.dmnd -o " + do_s + ".all.SP.cf.rd.dmnd.matches.txt")
            logfile.write("Done with DIAMOND analysis" + "\n")
            logfile.write("----------------------------------------------------" + "\n")

            # 	### turn diamond results into rma files viewable in MEGAN
            bash_command(
                "/home/hpglab/install/megan-ce/tools/blast2rma -i " + do_s + ".all.SP.cf.rd.dmnd.matches.txt -f BlastTab -r " + so_s + ".all.SP.cf.rd.fasta -o " + do_s + ".all.SP.cf.rd.dmnd.matches.rma -g2t /data/BLAST/gi2tax-July2016.bin")

            # 	### turn diamond output into Krona html file for visualization
            bash_command(
                "perl /data/install/KronaTools-2.6.1/scripts/ImportBLAST.pl " + do_s + ".all.SP.cf.rd.dmnd.matches.txt -o " + do_s + ".all.SP.cf.rd.dmnd.krona.html")
            # 	### runs kraken
            bash_command(
                "kraken --threads 16 --db /data/krakenDB/ " + so_s + ".all.SP.cf.rd.fasta | kraken-translate --db /data/krakenDB/ > " + ko_s + ".kraken.labels | cut -f2")

            bash_command("python /data/scripts/make_counts.py > " + ko_s + ".kraken4krona")

            # 	# output in krona html visual format
            bash_command("ktImportText " + ko_s + ".kraken4krona -o " + ko_s + ".kraken.krona.html")

    print "\nCompressing..."

    # Rezip files
    bar = progressbar.ProgressBar()
    for i in bar(range(len(filestorezip))):
        ftrz = filestorezip[i]
        if os.path.isfile(ftrz) and os.path.isfile(ftrz + ".gz"):
            os.remove(ftrz + ".gz")
        if os.path.isfile(ftrz):
            with open(ftrz, 'rb') as f_in, gzip.open(ftrz + ".gz", 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(ftrz)

    logfile.close()
    outfile.close()
    print "1_map.py complete."
    exit(0)













