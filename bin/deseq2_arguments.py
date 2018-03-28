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
import functions

##################################################################
# Added in 'plot_device =' as an extra line in the txt file output
##################################################################

if __name__ == '__main__':
    try: 
        outfile_name=sys.argv[1]
    except:
        outfile_name='deseq2_arguments.txt'
    lines=['indir = ', 'outdir = ', 'sample_info = ', 'comparisons = ', 'design = ', 'gtfFile = ', 'plot_device =']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()

    try: 
        outfile_name=sys.argv[2]
    except:
        outfile_name='deseq2_arguments.txt'
    lines=['indir = ', 'outdir = ', 'sample_info = ', 'comparisons = ', 'design = ', 'gtfFile = ', 'plot_device =']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()

    try: 
        outfile_name=sys.argv[3]
    except:
        outfile_name='sample_info.csv'
    lines=['SampleID,Group']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()

    try: 
        outfile_name=sys.argv[4]
    except:
        outfile_name='comparisons.csv'
    lines=['baselineGroup,comparisonGroup']
    outfile=open(outfile_name,'w')
    for l in lines:
        outfile.write(l + '\n')
    outfile.close()


    
