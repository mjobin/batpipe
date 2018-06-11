#!/usr/bin/python


#####################################
#####        HPG Lab            #####
#####    updated April 2018     #####
#####       MJJ                 #####
#####################################

# Written and maintained by Matthew Jobin, UCSC Human Paleogenomics Lab
# This runs the 0-4 scripts in succession

import argparse
from argparse import RawTextHelpFormatter
import subprocess
from subprocess import Popen, PIPE


def bash_command(cmd):
    subp = subprocess.Popen(['/bin/bash', '-c', cmd])
    subp.wait()
    return subp.returncode


def adddict(thedict):
    cl = ""
    for key in thedict:
        value = thedict[key]
        if type(value) is bool:
            if value:
                cl += " -"
                cl += str(key)
        else:
            cl += " -"
            cl += str(key)
            cl += " "
            cl += str(value)
    return cl


if __name__ == "__main__":
    print ("\n\n")
    print(
        "        __.--'\     \.__./     /'--.__\n    _.-'       '.__.'    '.__.'       '-._\n  .'                                      '.\n /                  BATPIPE                 \\\n|                                            |\n|                                            |\n \         .---.              .---.         /\n  '._    .'     '.''.    .''.'     '.    _.'\n     '-./            \  /            \.-'\n                      ''")

    parser = argparse.ArgumentParser(description="Script runs all 0-4 scripts:\n\t"
                                                 "1. You can choose whether to run 2b\n"
                                                 "2. Anything you put on the command line here will go to the appropriate scirpt\n"
                                                 "", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-twob', dest='twob', help='Run 2b_contammix.',
                        action='store_true')
    parser.set_defaults(fourb=False)
    parser.add_argument('-skipzero', dest='skipzero', help='Skip zero script.',
                        action='store_true')
    parser.set_defaults(skipzero=False)
    parser.add_argument('-bc_file', metavar='<bc_file>', help='location of barcode files, Must have a newline at end.',
                        required=True)
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
                        default="23")
    parser.add_argument('-bc_trim', metavar='<bc_trim>', help='Location of barcode trimmed files',
                        default='/data/barcodetrimmed')
    parser.add_argument('-seqprep_output', metavar='<seqprep_output>', help='seqprep_output',
                        default='/data/adaptertrimmed')
    parser.add_argument('-seqprep_output_in_output', metavar='<seqprep_output_in_output>',
                        help='Prepend output directory to seqprep_output',
                        default=False)
    parser.add_argument('-bformat', dest='bformat', help='Barcode format type new?.',
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

    # ONE
    parser.add_argument('-mtdna', metavar='<mtdna>', help='mtDNA mapping',
                        default='/data/genomes/rCRS.fas')
    parser.add_argument('-mia_ref', metavar='<mia_ref>', help='mia_ref',
                        default='rCRS')
    parser.add_argument('-index_algorithm', metavar='<index_algorithm>',
                        help='If reference is <2Gb use is, if >2Gb use bwtsw',
                        default='bwtsw')
    parser.add_argument('-seed_disable', metavar='<seed_disable>',
                        help='Following ancient DNA data processing protocols',
                        default="1024")
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
    parser.add_argument('-frag_length_r', metavar='<frag_length_r>', help='frag_length_r',
                        default='/data/scripts/frag_len_hist.R')
    parser.add_argument('-bwaindex', dest='bwaindex', help='Need to index if never used the reference genome before.',
                        action='store_true')
    parser.set_defaults(bwaindex=False)
    parser.add_argument('-nomia', dest='nomia', help='Do not run MIA.',
                        action='store_true')
    parser.set_defaults(nomia=False)
    parser.add_argument('-malt', dest='malt', help='Run MALT.',
                        action='store_true')
    parser.set_defaults(malt=False)
    parser.add_argument('-nosexest', dest='nosexest', help='Do not run sex estimation.',
                        action='store_true')
    parser.set_defaults(nosexest=False)
    parser.add_argument('-kraken', dest='kraken', help='Release the kraken.',
                        action='store_true')
    parser.set_defaults(kraken=False)
    parser.add_argument('-nodedup', dest='nodedup', help='Use samtools instead of DeDup to remove duplicates.',
                        action='store_true')
    parser.set_defaults(nodedup=False)
    parser.add_argument('-skipprinseq', dest='skipprinseq', help='Skip the prinseq part of script.',
                        action='store_true')
    parser.set_defaults(overwrite=False)
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
                        default="/data/db/malt/nr/")
    parser.add_argument('-maltblast', metavar='<maltblast>', help='MALT BLAST type.',
                        default="BlastX")
    parser.add_argument('-megandir', metavar='<megandir>', help='MEGAN6 directory',
                        default='/opt/megan')
    parser.add_argument('-krakendb', metavar='<krakendb>', help='KrakenDB location',
                        default='/data/db/krakenDB')
    parser.add_argument('-blastdir', metavar='<blastdir>', help='BLAST db directory',
                        default='/data/db/BLAST')
    parser.add_argument('-scriptsdir', metavar='<scriptsdir>', help='Default scripts directory',
                        default='/data/scripts')

    # TWO A
    parser.add_argument('-mt311', metavar='<mt311>', help='mt311',
                        default='/data/genomes/mt311.fna')

    # THREE
    parser.add_argument('-mpa_pkl', metavar='<mpa_pkl>', help='mpa_pkl',
                        default='/data/install/metaphlan2/db_v20/mpa_v20_m200.pkl')
    parser.add_argument('-bowtie2db', metavar='<bowtie2db>',
                        help='I thought this was called Bowie2DB and was momentarily happy',
                        default='/data/install/metaphlan2/db_v20/mpa_v20_m200')

    # FOUR
    parser.add_argument('-pmdref', metavar='<pmdref>', help='',
                        default="rCRS")
    parser.add_argument('-pmd_threshold', metavar='<pmd_threshold>', help='PMDtools threshold',
                        default="3")

    # RESULTS
    parser.add_argument('-raw', metavar='<raw>', help='Location of raw reads',
                        default='/data/raw')
    parser.add_argument('-mito_ref', metavar='<mito_ref>', help='mito_ref',
                        default='rCRS')
    parser.add_argument('-wobc_output', metavar='<wobc_output>', help='wobc_output',
                        default='/data/barcodetrimmed')

    thisdict = {}
    alldict = {}
    zerodict = {}
    onedict = {}
    resultsdict = {}
    twoadict = {}
    twobdict = {}
    threedict = {}
    fourdict = {}

    args = parser.parse_args()
    twob = bool(args.twob)
    skipzero = bool(args.skipzero)

    wd = args.wd
    alldict["wd"] = wd
    bcfile = args.bc_file
    alldict["bc_file"] = bcfile
    verbose = bool(args.verbose)
    alldict["verbose"] = verbose
    bcbinclip = args.bc_bin_clip

    rawreads = args.rawreads
    zerodict["rawreads"] = rawreads

    univ_illumina = args.univ_illumina
    zerodict["univ_illumina"] = univ_illumina

    ref_barcodes = args.ref_barcodes
    zerodict["ref_barcodes"] = ref_barcodes

    seqprep_min_length = args.seqprep_min_length
    zerodict["seqprep_min_length"] = seqprep_min_length

    seqprep_overlap = int(args.seqprep_overlap)
    zerodict["seqprep_overlap"] = seqprep_overlap

    mismatch = int(args.mismatch)
    zerodict["mismatch"] = mismatch

    bc_trim = args.bc_trim
    zerodict["bc_trim"] = bc_trim

    seqprep_output_in_output = bool(args.seqprep_output_in_output)
    zerodict["seqprep_output_in_output"] = seqprep_output_in_output

    seqprep_output = args.seqprep_output
    zerodict["seqprep_output"] = seqprep_output

    bformat = bool(args.bformat)
    zerodict["bformat"] = bformat

    fqloc = args.fqloc
    zerodict["fqloc"] = fqloc

    overwrite = bool(args.overwrite)
    zerodict["overwrite"] = overwrite

    threads = args.threads
    zerodict["threads"] = threads

    # 1
    mtdna = args.mtdna
    mia_ref = args.mia_ref
    index_algorithm = args.index_algorithm
    seed_disable = args.seed_disable
    q = args.q
    adna_matrix = args.adna_matrix
    max_misincorp_frequency = args.max_misincorp_frequency
    read_plot_length = args.read_plot_length
    max_read_length = args.max_read_length
    frag_length_r = args.frag_length_r
    bwaindex = bool(args.bwaindex)
    nomia = bool(args.nomia)
    malt = bool(args.malt)
    nosexest = bool(args.nosexest)
    kraken = bool(args.kraken)
    nodedup = bool(args.nodedup)
    skipprinseq = bool(args.skipprinseq)
    refs = args.refs

    arefs = ["/data/genomes/hg19.fa", "/data/genomes/rCRS.fas"]
    for ref in refs:
        arefs.append(ref)

    refs = ""
    for ref in arefs:
        refs = refs + ref
        refs = refs + " "

    sexref = args.sexref
    bwamaxedit = args.bwamaxedit
    clipleft = args.clipleft
    clipright = args.clipright
    maltdb = args.maltdb
    maltblast = args.maltblast
    megandir = args.megandir
    krakendb = args.krakendb
    blastdir = args.blastdir
    scriptsdir = args.scriptsdir

    onedict["mtdna"] = mtdna
    onedict["mia_ref"] = mia_ref
    onedict["index_algorithm"] = index_algorithm
    onedict["seed_disable"] = seed_disable
    onedict["threads"] = threads
    onedict["q"] = q
    onedict["adna_matrix"] = adna_matrix
    onedict["max_misincorp_frequency"] = max_misincorp_frequency
    onedict["read_plot_length"] = read_plot_length
    onedict["max_read_length"] = max_read_length
    onedict["frag_length_r"] = frag_length_r
    onedict["bwaindex"] = bwaindex
    onedict["nomia"] = nomia
    onedict["malt"] = malt
    onedict["nosexest"] = nosexest
    onedict["kraken"] = kraken
    onedict["nodedup"] = nodedup
    onedict["skipprinseq"] = skipprinseq
    onedict["overwrite"] = overwrite
    onedict["seqprep_output_in_output"] = seqprep_output_in_output
    onedict["seqprep_output"] = seqprep_output
    onedict["refs"] = refs
    onedict["sexref"] = sexref
    onedict["bwamaxedit"] = bwamaxedit
    onedict["clipleft"] = clipleft
    onedict["clipright"] = clipright
    onedict["maltdb"] = maltdb
    onedict["maltblast"] = maltblast
    onedict["bformat"] = bformat
    onedict["megandir"] = megandir
    onedict["krakendb"] = krakendb
    onedict["blastdir"] = blastdir
    onedict["scriptsdir"] = scriptsdir

    # TWO_A
    mt311 = args.mt311
    twoadict["mt311"] = mt311
    twoadict["overwrite"] = overwrite

    # TWO B
    twobdict["mt311"] = mt311
    twobdict["mia_ref"] = mia_ref
    twobdict["seqprep_output"] = seqprep_output
    twobdict["seqprep_output_in_output"] = seqprep_output_in_output
    twobdict["threads"] = threads
    twobdict["overwrite"] = overwrite

    # THREE
    mpa_pkl = args.mpa_pkl
    bowtie2db = args.bowtie2db
    threedict["mpa_pkl"] = mpa_pkl
    threedict["bowtie2db"] = bowtie2db
    threedict["seqprep_output"] = seqprep_output
    threedict["seqprep_output_in_output"] = seqprep_output_in_output
    threedict["overwrite"] = overwrite

    # FOUR
    pmd_threshold = args.pmd_threshold
    pmdref = args.pmdref
    fourdict["pmdref"] = pmdref
    fourdict["q"] = q
    fourdict["pmd_threshold"] = pmd_threshold
    fourdict["overwrite"] = overwrite

    # RESULTS
    raw = args.raw
    mito_ref = args.mito_ref
    wobc_output = args.wobc_output

    resultsdict["raw"] = raw
    resultsdict["seqprep_output"] = seqprep_output
    resultsdict["seqprep_output_in_output"] = seqprep_output_in_output
    resultsdict["wobc_output"] = wobc_output
    resultsdict["nomia"] = nomia
    resultsdict["refs"] = refs
    resultsdict["bformat"] = bformat

    runline = "0_trim_barcodes.py"
    runline += adddict(alldict)
    runline += adddict(zerodict)
    if not skipzero:
        rc = bash_command(runline)
        if rc != 0:
            print "Error in script 0! Halting. Check logs!"
            exit()

    runline = "1_map.py"
    runline += adddict(alldict)
    runline += adddict(onedict)
    rc = bash_command(runline)
    if rc != 0:
        print "Error in script 1! Halting. Check logs!"
        exit()

    runline = "2a_contam.py"
    runline += adddict(alldict)
    runline += adddict(twoadict)
    rc = bash_command(runline)
    if rc != 0:
        print "Error in script 2a! Halting. Check logs!"
        exit()

    if twob:
        runline = "2b_contammix.py"
        runline += adddict(alldict)
        runline += adddict(twobdict)
        rc = bash_command(runline)
        if rc != 0:
            print "Error in script 2b! Halting. Check logs!"
            exit()

    runline = "3_metaphlan.py"
    runline += adddict(alldict)
    runline += adddict(threedict)
    rc = bash_command(runline)
    if rc != 0:
        print "Error in script 3! Halting. Check logs!"
        exit()

    runline = "4_pmdtools.py"
    runline += adddict(alldict)
    runline += adddict(fourdict)
    rc = bash_command(runline)
    if rc != 0:
        print "Error in script 4! Halting. Check logs!"
        exit()

    runline = "results.py"
    runline += adddict(alldict)
    runline += adddict(resultsdict)
    rc = bash_command(runline)
    if rc != 0:
        print "Error in script results! Halting. Check logs!"
        exit()

    print "batpipe.py complete"

    exit(0)







