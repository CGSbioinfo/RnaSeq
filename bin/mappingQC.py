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
#sys.path.insert(0,'/usr/local/bin/')
import functions
import argparse

__version__ = 'v02'

def junctions(i):
    os.system("junction_annotation.py -i " + 
        in_dir + "/" + i +"Aligned.sortedByCoord.out.bam -o " + 
        out_dir + '/' + i + " -r " + bedFile + " & ")
    os.system("junction_saturation.py -i " + 
        in_dir + '/' + i +"Aligned.sortedByCoord.out.bam  -o " + 
        out_dir + '/' + i + " -r " + bedFile)


# Collect Metrics
def picard_collect_metrics(i):
    bamfiles = os.listdir(in_dir)
    bamfiles = [bamfiles[y] for y, x in enumerate(bamfiles) if re.findall(i, x)]
    bamfiles = [bamfiles[y] for y, x in enumerate(bamfiles) if re.findall(r'sortedByCoord\.out\.bam$', x)][0]
    os.system("java -jar /usr/local/picard/picard.jar CollectRnaSeqMetrics " + 
        " REF_FLAT=" + refFlat + 
        " RIBOSOMAL_INTERVALS=" + rRNA_interval_list + 
        " STRAND_SPECIFICITY=" + strand_piccard + 
        " INPUT=" + in_dir + "/" + bamfiles +  
        " OUTPUT=" + out_dir + "/" + i + "_metrics.txt")

def pct(i):
    os.system('/usr/bin/Rscript bin/mapping_distribution.R ' + out_dir + '/ ' + i)


#########################

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='mappingQC.py',description = 'Quality control of mapped data')
    parser.add_argument('-v','--version',action='version',version='%(prog)s-'+__version__)
    parser.add_argument('--analysis_info_file', help='Text file with details of the analysis. Default=analysis_info.txt', default='analysis_info.txt')
    parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=alignedReads/', default='alignedReads/')
    parser.add_argument('--out_dir', help='Path to out put folder. Default=alignedReads/QC/', default='alignedReads/QC/')
    parser.add_argument('--out_dir_report', help='Path to out put folder. Default=Report/figure/mappingQC/', default='Report/figure/mappingQC/')
    parser.add_argument('--sample_names_file', help='Text file with sample names. Default=sample_names_info.txt', default='sample_names.txt')
    #parser.add_argument('--suffix_name', help='Suffix to optionally put to the output name. Default=', default='_plot')
    parser.add_argument('--run', help='Choose a section of the pipeline to run. Possible options: mapping_summary; gene_body_coverage; junctions; picard_tools; all. Default = all')
    parser.add_argument('--plot_device', action='store', help='Specify the format of the plot output. Default=png', default='png')
    #parser.add_argument('--ncores', help='Number of cores to use. Default=8', default='8')
    args=parser.parse_args()

    ai=functions.read_analysis_info_file(args.analysis_info_file)
    path=os.getcwd()
    refGenome=ai['reference_genome']
    os.chdir(path)
    bedFile_10k=ai['BedFile10K']
    bedFile=ai['BedFile']
    refFlat=ai['refFlat']
    rRNA_interval_list=ai['rRNA_interval_list']
    strand=ai['strand']
    strand_piccard, strand_htseq = functions.get_strand(strand)
    plot_device = args.plot_device
   # suffix_name = args.suffix_name
    
    #Ncores
    ncores=int(ai['ncores'])

    # Read sample names text file
    sample_names_file=args.sample_names_file
    sampleNames = functions.read_sample_names(sample_names_file)

    # Set input and output directories if not '/'
    in_dir=path + '/' + args.in_dir
    out_dir=path + '/' + args.out_dir
    functions.make_sure_path_exists(out_dir)
    out_dir_report=path + '/' + args.out_dir_report
    functions.make_sure_path_exists(out_dir_report)

    # Detect if files are gz
    #gz = functions.check_gz(in_dir)
    
    # Mapping summary
    if args.run == 'all' or args.run == 'mapping_summary':
        os.system("/usr/bin/Rscript bin/mapping_summary.R " + in_dir + '/ ' + out_dir + '/' )
        os.system("/usr/bin/Rscript bin/mapping_summary.R " + in_dir + '/ ' + out_dir_report + '/' )

    # Gene body coverage
    if args.run == 'all' or args.run == 'gene_body_coverage':
        os.system("ls " + in_dir + "/*Aligned.sortedByCoord.out.bam > tempbamfiles.txt")
        os.system("/usr/local/minicondaexport/bin/geneBody_coverage.py -r " + bedFile_10k + " -i tempbamfiles.txt -o " + out_dir + "/10KGenes")
        os.system("rm tempbamfiles.txt")
        os.system('cp ' + out_dir + '/10KGenes.geneBodyCoverage.curves.pdf ' + out_dir_report + '/10KGenes_geneBodyCoverage_curves.pdf')

    # Junction QC
    if args.run == 'all' or args.run == 'junctions':
        Parallel(n_jobs=ncores)(delayed(junctions)(i) for i in sampleNames)
        os.system('grep "y=c(" ' + out_dir + '/*junctionSaturation*  | sed \'s/:y=c(/,/g\' | sed \'s/.junctionSaturation_plot.r//g\' | sed \'s/)//g\' | sed \"s/.*\///g\"  > ' + out_dir + '/junctionSat_all.csv')
        os.system('Rscript ~/bin/junctionPlotAll.R ' + out_dir + ' ' + out_dir)
        os.system('cp ' + out_dir + '/junctionSaturationAll.pdf ' + out_dir_report)

    # Picard tools
    if args.run == 'all' or args.run == 'picard_tools':
       # Parallel(n_jobs=ncores)(delayed(picard_collect_metrics)(i) for i in sampleNames)
       # Parallel(n_jobs=ncores)(delayed(pct)(i) for i in sampleNames)
        os.system('/usr/bin/Rscript bin/read_distribution_genomic_context.R ' + out_dir + ' ' + out_dir_report+' '+plot_device )
