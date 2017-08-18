#!/usr/bin/env python
import sys

# USAGE: python fastq_data.txt per_base_sequence_content.py
fastq_data = sys.argv[1].replace('\\','')
plot = sys.argv[2] # all, per_base_sequence_content, per_base_sequence_quality, per_sequence_gc_content, per_sequence_quality_scores, seq_dup_levels, kmer_content
outdir = sys.argv[3] + '/'

dictionary = {'seq_length':'>>Sequence Length Distribution','per_base_sequence_content':'>>Per base sequence content', 'per_base_sequence_quality':'>>Per base sequence quality', 'kmer_content':'>>Kmer Content', 'seq_dup_levels':'>>Sequence Duplication Levels', 'per_sequence_quality_scores':'>>Per sequence quality scores', 'per_sequence_gc_content': '>>Per sequence GC content'}

if plot == 'all':
    for d in dictionary.keys():
        f = open(fastq_data,'r')
        i=0
        for line in f:
            line = line.strip()
            if line.startswith(dictionary[d]):
                name = d + '.txt'
                line_number = i
            i += 1
        f.close()
        output = open(outdir + name,'w')
        f = open(fastq_data,'r')
        i=0
        for line in f:
            line = line.strip()
            if i > line_number:
                if line != '>>END_MODULE':
                    output.write(line + '\n')
                else:
                    break
            i += 1
        f.close()
else:
    f = open(fastq_data,'r')
    i=0
    for line in f:
        line = line.strip()
        if line.startswith(dictionary[plot]):
            name = plot + '.txt'
            line_number = i
        i += 1
    f.close()
    output = open(outdir + name,'w')
    f = open(fastq_data,'r')
    i=0
    for line in f:
        line = line.strip()
        if i > line_number:
            if line != '>>END_MODULE':
                output.write(line + '\n')
            else:
                break
        i += 1

