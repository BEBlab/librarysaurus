# Librarysaurus

A python script for a quick analysis of plasmidsaurus results on your library miniprep.

This script uses a control sequence for (roughly) analyzing the sequencing error rate.
Then it analyzes your sequence(s) of interest in the region where you cloned it.

## Usage:

Place the script, the input file (plasmidsaurus locate this file in the "per-base-data" folder).
Run through the terminal. Recommended to be run in conda environment. 

If you run 'python Plasmidsaurus_stats_analysis.py -h' in your terminal you should get:

usage: Librarysaurus.py [-h] [-m {DMS,pool}] [-i INFILE] [-p STATS_FILE]

Analyze Plasmidsaurus Per Base Stats files in the context of library cloning

optional arguments:
  -h, --help            show this help message and exit
  -m {DMS,pool}, --mode {DMS,pool}
                        Choose what you want to analyze
  -i INFILE, --input INFILE
                        Input File (see example_file.txt)
  -p STATS_FILE, --plasmidsaurus STATS_FILE
                        Plasmidsaurus Stats file

mode: it depends on the type of your library. Can be a NNK library (single/doubles mutations) or a pool (with different sequences).
INFILE: is '.txt' tab separated file where you put the control sequence + adaptor sequence 5' + adaptor sequence 3', and your sequence(s) of interest + adaptor sequence 5' + adaptor sequence 3'. (see example.txt).
STATS_FILE: stats .tsv file from Plasmidsaurus located the "per-base-data" folder.

## Dependencies:

Python 3 and following packages needed: Pandas, Matplotlib, Seaborn, Biopython, Argparse, Os, Warnings
Recommended to be run in conda environment.

## Out:

Two bar plots showing the mutation rate per nucleotide position:
-one for the control sequences
-one for your sequence of interest

