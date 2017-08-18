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
#sys.path.insert(0,'./')
import functions
import argparse

__version__='v01'
# created 17/08/2016

if __name__ == '__main__':

    # Parser
    parser = argparse.ArgumentParser(prog='organizeWorkingDirectory.py',description = 'Organize working directory of the analysis')
    parser.add_argument('-v','--version', action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    parser.add_argument('--in_dir', help='directory with fastq files. Default= corresponding to bcl2fastq_output', default='bcl2fastq_output')
    args=parser.parse_args()

    # Read analysis info file
    ai=functions.read_analysis_info_file(args.analysis_info_file)
    
    # Change dir
    os.chdir(ai['project_location'])

    # Read sample names
    sample_names_file = args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)
        
    # Create rawReads folder
    # Check if rawReads exists
    project_location=ai['project_location']
    folders = os.listdir(project_location)
    readsFiles = [folders[i] for i, x in enumerate(folders) if re.findall('rawReads',x)]
    # print readsFiles    

    # Collect fastq files analysis_info_file
    if args.in_dir == 'bcl2fastq_output':
        allFiles=functions.get_filepaths(ai['project_location'] + '/' + ai[args.in_dir])
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[fastq[y] for y,x in enumerate(fastq) if not re.findall('Undetermined', x)]
    elif args.in_dir != 'bcl2fastq_output':
        allFiles=os.listdir(args.in_dir)
        fastq=[allFiles[y] for y, x in enumerate(allFiles) if re.findall("fastq.gz", x)]
        fastq=[args.in_dir + x for x in fastq]
    # print fastq
    
    # Move reads
    if not readsFiles:
        functions.make_sure_path_exists('rawReads')
        sampleDir = []
        for sample in sampleNames:
            reads = [fastq[i] for i,x in enumerate(fastq) if re.findall(sample,x)]
            if sample not in sampleDir:
                functions.make_sure_path_exists('rawReads/'+sample)
            for r in reads:
               os.system('mv ' + '"' + r + '"' + ' rawReads/' + sample)
            sampleDir.append(sample)
    else:
        print "rawReads/ already folder exists"

