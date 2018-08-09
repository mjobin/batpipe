# batpipe
NGS/aDNA data analysis pipeline for the Human Paleogenomics Lab at UCSC.
Batpipe is the ancient DNA-focused Next-Generation Sequencing data analysis pipeline developed by and used at Lars Fehren-Schmitz’s Human Paleogenomics Lab at the University of California Santa Cruz. The pipeline was originally developed as a series of bash scripts by Kelly Harkins and Pete Heintzman, and has been converted to Python, extended and is currently maintained by Matthew Jobin. 
The pipeline can be run in separate stages by invoking each of the scripts 0 through 4 individually. Alternatively, the batpipe.py script can be invoked to run the desired scripts in turn. Note that scripts 2a through 4 normally require scripts 0 through 1 to already have been run in the working directory.
The pipeline proceeds through trimming barcodes to mapping the resultant sequences to desired references, after which two tests for contamination on the mitochondria may be run, followed by metagenome analysis and an assessment of post-mortem damage. Many of the functions of the pipeline involve the invocation of other, separate software programs, for which a list is provided in the wiki link below. It is important for the new user to compile and/or install these programs and verify that they are running before using Batpipe. Where possible, Batpipe assumes that external software is located within the search path of the local machine, so it is thus required for the user to make sure this is true for this software on the machine he or she is using. Where external software cannot be placed in a path, a header variable is put in place in order to ensure proper invocation. The user may either install that software in the default location specified (i.e. /data/scripts) or change the location using the command line.
The user may set the working directory and locations of the input files. Input is assumed to be in FASTQ format. There are also default locations for intermediate files such as files with barcodes trimmed. The user may change these on the command line, or edit the code to change the default.

For full details regarding the operation of the pipeline, please consult: https://sites.google.com/a/ucsc.edu/hpglab-wiki/

Please note! Batpipe integrates numerous software programs written neither by me nor by the UCSC Human Paleogenomics Lab. These programs must be installed and working on your system for batpipe to operate. Here follows a list of the needed software:

bc_bin_clip    by Sam Vohr
SeqPrep2    https://github.com/jeizenga/SeqPrep2
