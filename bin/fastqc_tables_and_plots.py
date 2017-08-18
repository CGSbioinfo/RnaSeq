#!/usr/bin/env python

import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import multiprocessing
import subprocess
#sys.path.insert(0,'/usr/local/bin')
import functions
import argparse

def tables(i):
    outdir = re.sub('fastqc_data.txt', '', i)
    os.system("python bin/create_fastqcTables.py " + i + " all " + outdir)

def plots(i):
    #outdir = re.sub('fastqc_data.txt', '', i)
    os.system('/usr/bin/Rscript bin/create_fastqcPlots_perSample.R ' + in_dir + ' ' + i + ' ' + readType + ' ' + out_dir_report + ' ' + suffix_name + ' ' + args.plot_device)


##############################################################

__version__ = 'v02'

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='fastqc_tables_and_plots.py', description = 'Generates tables used for FastQC plots')
    parser.add_argument('-v','--version', action='version', version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads')
    #parser.add_argument('--out_dir', help='Path to out put folder. Default=rawReads/', default='rawReads')
    parser.add_argument('--readType', help='Read Type: pairedEnd or singleEnd. Default=pairedEnd', default='pairedEnd')
    parser.add_argument('--out_dir_report', help='Path to out put folder. Default=Report/figure', default='Report/figure')
    parser.add_argument('--suffix_name', help='Suffix to optionally put to the output name. Default=', default='_plot')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names.txt', default='sample_names.txt')
    parser.add_argument('--plot_device', action='store', help='Specify the format of the plot output. Default=png', default='png')
    #parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    # Read analysis info file and cd
    ai=functions.read_analysis_info_file(args.analysis_info_file)
    os.chdir(ai['project_location'])    

    #Ncores
    ncores=ai['ncores']

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not 'rawReads/'
    in_dir= args.in_dir
    out_dir_report= args.out_dir_report
    readType=ai['readType']
    suffix_name=args.suffix_name
    
    # Create tables
    files=functions.get_filepaths(in_dir)
    files = [files[y] for y, x in enumerate(files) if re.findall("fastqc_data.txt", x)] 
    Parallel(n_jobs=8)(delayed(tables)(i) for i in files)
    print "Got data from fastqc output... \n"    

    # Create plots
    functions.make_sure_path_exists(out_dir_report)
    Parallel(n_jobs=8)(delayed(plots)(i) for i in sampleNames)
    print "Made plots per sample... \n"
    os.system('/usr/bin/Rscript bin/create_fastqcPlots_allSamples.R ' + in_dir + ' ' + sample_names_file + ' ' + readType + ' ' + out_dir_report + ' ' + suffix_name + ' ' + args.plot_device)
    print "Made plots all samples... \n"
