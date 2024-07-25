# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 10:27:20 2023

@author: mmartin
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 18 16:38:51 2022

@author: mmartin
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import os
import warnings

warnings.filterwarnings("ignore")


def stats_df(stats_file):
    df = pd.read_csv(stats_file, sep="\t")
    df['MutRate'] = df.mismatches/df.reads_all
    return df

def reverse_sequence(df, ctt_5, ctt_3):
    ctt_5_rev = Seq(ctt_3).reverse_complement()
    start_p = ('').join(df.ref.tolist()).find(str(ctt_5_rev))+len(ctt_5_rev)
    ctt_3_rev = Seq(ctt_5).reverse_complement()
    end_p = ('').join(df.ref.tolist()).find(str(ctt_3_rev))
    return start_p, end_p

def start_end_p(df, wt_dna, ctt_5, ctt_3):
    start_p = ('').join(df.ref.tolist()).find(ctt_5)+len(ctt_5)
    end_p = ('').join(df.ref.tolist()).find(ctt_3)
    if end_p == -1 or start_p == -1:
        start_p, end_p = reverse_sequence(df, ctt_5, ctt_3)
    return start_p, end_p

def regioninterest(mode, df, wt_dna, ctt_5, ctt_3):
    start_p, end_p = start_end_p(df, wt_dna, ctt_5, ctt_3)
    region_interest = df.iloc[start_p:end_p]
    region_interest['comp'] = [x for x in str(Seq(('').join(region_interest.ref.tolist())).complement())]
    region_interest = region_interest.sort_values(by='pos')
    region_interest.reset_index(inplace=True)
    region_interest['pos'] = region_interest.index.tolist()
    consensus = ('').join(region_interest.ref.tolist())
    rev_consensus = ('').join(region_interest.comp.tolist())
    if wt_dna == consensus:
        print('Consensus sequence is the WT sequence of interest')
    elif wt_dna == rev_consensus:
        print('Consensus sequence is the WT sequence of interest')   
    elif mode == 'pool':
        print('Analyzing pools...')
        if len(wt_dna) == len(region_interest):
            print(f'Inserts length is the same as the designed sequences:', len(region_interest), 'nucleotides')
        else:
            print("There's a difference between the length of designed inserts and what we have in plasmidsaurus seqs...")
            print("However this is not necessaily bad and could be related to the sequencing errors...")
    else:
        print('The consensus sequence is not the WT sequence of interest!!')
        print(wt_dna)
        print(consensus)
        region_interest = pd.DataFrame()
        start_p = 0
        end_p = 0
    return region_interest, start_p, end_p
    

def plot_hist(mode, ID, library_name, region_interest, df, start_p, wt_dna):
    plt.figure(figsize=(20,5))
    sns.barplot(data=region_interest, x='pos', y='MutRate')
    plt.axhline(y=df.iloc[:start_p].MutRate.mean(), linestyle='--', label='Mean Error Rate')
    if mode == 'DMS':
        plt.xticks(ticks=range(len(region_interest.pos)), labels = wt_dna, rotation=0)
        plt.legend()
        plt.savefig(ID+'_'+library_name+'_plasmidsaurus_seq_error_rate.png')
    elif mode == 'pool':
        plt.xticks(ticks=range(0, len(region_interest.pos), 10), labels = range(0, len(region_interest.pos), 10), rotation=45)
        plt.legend()
        plt.savefig(library_name+'_plasmidsaurus_seq_error_rate.png')
    plt.show()
    

def run_analysis(mode, library_name, in_file, stats_file):
    in_file_df = pd.read_csv(in_file, sep="\t")
    df = stats_df(stats_file)
    if mode == 'pool':
        ID = in_file_df.loc[0, 'ID']
        wt_dna = in_file_df.loc[0, 'Sequence']
        ctt_5 = in_file_df.loc[0, 'ctt_5']
        ctt_3 = in_file_df.loc[0, 'ctt_3']   
        region_interest, start_p, end_p = regioninterest(mode, df, wt_dna, ctt_5, ctt_3)
        if region_interest.empty:
            print('')
        else:
            print('Showing plot of Mut Rate vs Pos of insert')
            plot_hist(mode, ID, library_name, region_interest, df, start_p, wt_dna)
    else:
        for i in in_file_df.index:
            ID = in_file_df.loc[i, 'ID']
            wt_dna = in_file_df.loc[i, 'Sequence'].upper()
            ctt_5 = in_file_df.loc[i, 'ctt_5'].upper()
            ctt_3 = in_file_df.loc[i, 'ctt_3'].upper() 
            region_interest, start_p, end_p = regioninterest(mode, df, wt_dna, ctt_5, ctt_3)
            if region_interest.empty:
                print('')
            else:
                print('Showing plot of Mut Rate vs Pos of insert')
                plot_hist(mode, ID, library_name, region_interest, df, start_p, wt_dna)
    print("That's all folks!")
    return region_interest, df

if __name__ == "__main__":
    width = os.get_terminal_size().columns
    stars = '*'*85
    print('')
    print(stars.center(width))
    print("**           LIBRARYSAURUS           **".center(width))  
    print("**           Hi! this script will do the Plasmidsaurus sequence analysis           **".center(width))
    print("** Let's see if the region of interest is cloned and how's the mutation rate there **".center(width))
    print(stars.center(width))
    print("(by Mariano Martín - last update July 2024 - with tons of BEB lab input)".center(width))
    print('')
    print("Check in the running folder the plots generated!".center(width))
    parser = argparse.ArgumentParser(exit_on_error=False, 
                                     description='Analyze Plasmidsaurus Stats files in the context of library cloning')
    parser.add_argument('-m', '--mode', choices=['DMS', 'pool'], help='Choose what you want to analyze')
    parser.add_argument('-n', '--name', dest='library_name', action = "store", help = "Library Name")         
    parser.add_argument('-i', '--input',dest='infile', action = "store", help = "Input File (see example_file.txt)")
    parser.add_argument('-p', '--plasmidsaurus', dest='stats_file', action = "store", help = "Plasmidsaurus stats txt file")         
    args = parser.parse_args()  
    if args.mode != None and args.infile != None and args.stats_file != None:
        if args.mode == 'DMS':
            print('You designed a DMS library and want to analyze the plasmidsaurus results')
            print('')          
        elif args.mode == 'pool':
            print('You designed an oPool with several sequences and want to analyze the plasmidsaurus results...')
            print('')
        mode = args.mode
        library_name = args.library_name
        in_file = args.infile  
        stats_file = args.stats_file
        run_analysis(mode, library_name, in_file, stats_file)
    else:
        print('There is an error on the arguments passed!')
        parser.print_help()
    