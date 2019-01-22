RnaSeq pipeline
================

-   [Getting Started](#getting-started)
-   [Running the pipeline](#running-the-pipeline)
    -   [Step 1: Analysis info file](#step-1-analysis-info-file)
        -   [Format of the analysis info file](#format-of-the-analysis-info-file)
        -   [Links to tools](#links-to-tools)
        -   [Parameter choices for small RNA seq projects and for Lexogen projects](#parameter-choices-for-small-rna-seq-projects-and-for-lexogen-projects)
        -   [Strand](#strand)
        -   [Reference genome](#reference-genome)
        -   [How to create the analysis info file](#how-to-create-the-analysis-info-file)
    -   [Step 2: Obtain FASTQ files](#step-2-obtain-fastq-files)
        -   [Using bcl2fastq](#using-bcl2fastq)
        -   [Downloading data from basespace (Currently used for Lexogen projects)](#downloading-data-from-basespace-currently-used-for-lexogen-projects)
        -   [Move the reads to a new folder named rawReads](#move-the-reads-to-a-new-folder-named-rawreads)
        -   [Create sample names file](#create-sample-names-file)
    -   [Step 3: Quality control of Fastq files](#step-3-quality-control-of-fastq-files)
    -   [Step 4: Table and plot of number of reads per sample](#step-4-table-and-plot-of-number-of-reads-per-sample)
    -   [Step 5: FastQC plots](#step-5-fastqc-plots)
    -   [Step 6: Trim low quality bases and adapters](#step-6-trim-low-quality-bases-and-adapters)
        -   [Step 6 (Optional for Lexogen): Trim polyA sequences](#step-6-optional-for-lexogen-trim-polya-sequences)
    -   [Step 7: Trimming summary](#step-7-trimming-summary)
    -   [Step 8: Trimming QC plots](#step-8-trimming-qc-plots)
    -   [Step 9: Mapping](#step-9-mapping)
    -   [Step 10: Mapping QC](#step-10-mapping-qc)
        -   [10a: Mapping summary](#a-mapping-summary)
        -   [10b: Gene body coverage](#b-gene-body-coverage)
        -   [10c: Genomic context of mapped reads.](#c-genomic-context-of-mapped-reads.)
    -   [Step 11: Counting of reads](#step-11-counting-of-reads)
    -   [Step 12: Counting QC](#step-12-counting-qc)
        -   [12a: Counting QC part1](#a-counting-qc-part1)
        -   [12b: Counting QC part2](#b-counting-qc-part2)
    -   [Step 13: EdgeR and DESeq2](#step-13-edger-and-deseq2)
        -   [13a: Create requiered files:](#a-create-requiered-files)
        -   [13b: DESeq2 PCA plot](#b-deseq2-pca-plot)
        -   [13c: Differential gene expression.](#c-differential-gene-expression.)

Getting Started
---------------

1.  Create a new folder in fs3 or any other location with the project name - called the "project folder" in this README.
2.  Create a bin/ folder in the project folder.
3.  Download and copy the scripts to the bin/ folder.

Running the pipeline
--------------------

### Step 1: Analysis info file

A central part of the pipeline is the **analysis info** file. It has information about the project location, the original run folder, the reference fasta, gtf and bed files, and the parameters used throughout the analysis.

#### Format of the analysis info file

The analysis info file is a simple .txt file, with each line providing information. Parameters are separated by the semicolon (i.e ";") character.

The script that creates the analysis info file chooses some parameters by default. Below there is more information and links to the tools used if you want to change the parameters.

The following is an example of the analysis info file:

|                                                                                                         |
|:--------------------------------------------------------------------------------------------------------|
| **working\_directory =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/      |
| **run\_folder =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/xxx          |
| **run\_samplesheet =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/xxx     |
| **bcl2fastq\_output =** /mnt/cgs-fs3/Sequencing/Pipelines/RNASeq\_pipeline/example\_small\_files/fastq/ |
| **readType =** pairedEnd                                                                                |
| **reference\_genome =**                                                                                 |
| **trimgalore\_params =** --gzip; --paired; --fastqc; --fastqc\_args '--nogroup --extract'               |
| **mapping\_params =** --runThreadN 4; --outSAMtype BAM SortedByCoordinate; readFilesCommand zcat        |
| **ncores =** 8                                                                                          |
| **htseq\_params =** -a 10; -m union; -s reverse; -t exon                                                |
| **strand =** reverse                                                                                    |

The following is the explanation of the analysis info file:

|                                                                                                                                     |
|:------------------------------------------------------------------------------------------------------------------------------------|
| **working\_directory =** *&lt;path to directory of the analysis&gt;*                                                                |
| **run\_folder =** *&lt;path to the run folder&gt;*                                                                                  |
| **run\_samplesheet =** *&lt;sample sheet to be used to generate fastq files. This is created using the Illumina Expert Manager&gt;* |
| **bcl2fastq\_output =** *&lt;path to output of bcl2fastq. The defaults is fastq/ and the folder will be created automatically&gt;*  |
| **readType =** *&lt;either pairedEnd or singleEnd&gt;*                                                                              |
| **reference\_genome =** *&lt;path to the STAR reference genome folder that will be used at the mapping step&gt;*                    |
| **trimgalore\_params =** *&lt;parameters to be passed to trim galore&gt;*                                                           |
| **mapping\_params =** *&lt;parameters to be passed to star&gt;*                                                                     |
| **ncores =** *&lt;Number of cores to use to pararellize analysis&gt;*                                                               |
| **htseq\_params =** *&lt;parameters to be passed to htseq-count&gt;*                                                                |
| **strand =** *&lt; expected mapping strand. &gt;*                                                                                   |

<br>

#### Links to tools

-   Trim galore: <https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md>
-   STAR: <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>
-   HTSeq: <http://htseq.readthedocs.io/en/release_0.9.0/count.html>

#### Parameter choices for small RNA seq projects and for Lexogen projects

**For small RNA seq projects:**
- consider changing the mapping parameters in analysis info to --runThreadN 4; --outSAMtype BAM SortedByCoordinate; --readFilesCommand zcat; --outFilterMultimapNmax 20; --outFilterMismatchNoverLmax 0.03; --outFilterScoreMinOverLread 0; --outFilterMatchNminOverLread 0; --outFilterMatchNmin 16; --alignSJDBoverhangMin 1000; --alignIntronMax 1.
- consider using the gtf file from mirBase.
- **If the library is single end and maps to the forward strand, and you are using the gtf file from mirbase,** change the counting parameters to -a 0; -m union; -s yes; -t miRNA; -i Name.

**For lexogen projects:**
- you need to add the argument *--clip\_R1 12* in the trimgalore parameters in the analysis info file to remove the first 12 bases as recommended.

<br>

#### Strand

This depends on the format of the data: whether is strand specific or not, and on which strand are reads expected to map.

Strand is specified for the mapping QC section where picard tools is used. Picard options are: NONE, FIRST\_READ\_TRANSCRIPTION\_STRAND, and SECOND\_READ\_TRANSCIRPTION\_STRAND. Here is a link to more explanation, look in the CollectRnaSeqMetrics section: <http://broadinstitute.github.io/picard/command-line-overview.html>. The following is some brief explanation:

|                                                                                                                                                                                                                                                                      |
|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| For strand-specific library prep. For unpaired reads, use FIRST\_READ\_TRANSCRIPTION\_STRAND if the reads are expected to be on the transcription strand. Required. Possible values: {NONE, FIRST\_READ\_TRANSCRIPTION\_STRAND, SECOND\_READ\_TRANSCRIPTION\_STRAND} |

-   **When we do paired end analysis** we normally use SECOND\_READ\_TRANSCRIPTION\_STRAND.
-   **When we do single end analysis** we normally use FIRST\_READ\_TRANSCRIPTION\_STRAND.

The value chosen should also match the strand specified in the htseq\_params field. The specific argument fot htseq is *-s*, which could be: *yes*, *no* or *reverse*. <http://htseq.readthedocs.io/en/release_0.9.0/count.html>

-   **When we do paired end analysis** we normally use reverse.
-   **When we do single end analysis** we normally use yes.

<br>

#### Reference genome

The reference genome should be the output of STAR's *genomeGenerate* command. We currently download genomes from [ensembl](ftp://ftp.ensembl.org/pub/).

<br>

#### How to create the analysis info file

**Script:** bin/analysis\_info.txt.

**Command example:**

``` bash
$ python bin/analysis_info.py --outfile analysis_info.txt
```

This will create a file analysis\_info.txt, which you can open in a text editor and fill.

<br>

### Step 2: Obtain FASTQ files

#### Using bcl2fastq

**Script:** bin/run\_bcl2fastq.py

**Input:** the analysis\_info.txt file, which has information about the **run folder**, the **run samplesheet**, and the **bcl2fastq\_output** folder. If you have more than one run, you can run this command for each run by changing the run\_folder and run\_samplesheet fields in the analysis\_info.txt file).

**Command:**

``` bash
$ python bin/run_bcl2fastq.py --analysis_info_file analysis_info.py   
```

**Output:** fastq/ folder, which has a project subfolder with the fastq.gz files, among other files and subfolders.

#### Downloading data from basespace (Currently used for Lexogen projects)

For lexogen projects, we currently download the data from Basespace. Normally, the downloaded data will have the following structure:

-   A main folder with the project name
-   Subfolders, which correspond to each of the sequenced samples
-   Within each subfolder, four fastq.gz files from the corresponding sample

The tasks done are:

-   Unzip the fastq files:
    -   Move to the downloaded folder.
    -   Use the following commands:

``` bash
$ ls -d * | parallel -j 4 --no-notice "cd {} ; gzip -d *"      
```

-   Change the folder name of each sample to a string consisting of **the sample name, an underscore, a capital S, and one or more numbers (eg. sample1\_S1)**. You can use the following commands, however **make sure that the sample names consists of the original sample name, an underscore, a capital S, and one or more numbers (eg. sample1\_S1)**:
    -   Check what the outcome would be. This command will print the original folder name, followed by a space and the new folder name:

``` bash
$ for i in $(ls -d *); do sample_name=$(ls ${i} | sed 's/.*\///g' | sed 's/_L00[[:digit:]]\+.*//g'  | sort | uniq); echo ${i} ${sample_name}; done         
```

    -   If it looks correct, use mv to change the folder names:    

``` bash
$ for i in $(ls -d *); do sample_name=$(ls ${i} | sed 's/.*\///g' | sed 's/_L00[[:digit:]]\+.*//g'  | sort | uniq); mv ${i} ${sample_name}; done     
```

-   For each sample, concatenate the four files into one:
    -   The following command tests whether the outcome of the command is right. It prints the sample name, it's four fasta files, and the name of the output concatenated file.

``` bash
$ ls -d * | parallel -j 4 --no-notice "echo {}; ls {}/{}*L001_R1*  {}/{}*L002_R1* {}/{}*L003_R1*  {}/{}*L004_R1*; echo {}/{}_R1_001.fastq"     
```

    -   If that looks correct, concatenate the files for R1 reads:   

``` bash
$ ls -d * | parallel -j 4 --no-notice "cat {}/{}*L001_R1*  {}/{}*L002_R1* {}/{}*L003_R1*  {}/{}*L004_R1* > {}/{}_R1_001.fastq"     
```

    -   And for R2 reads for paired end reads:    

``` bash
$ ls -d * | parallel -j 4 --no-notice "cat {}/{}*L001_R2*  {}/{}*L002_R2* {}/{}*L003_R2*  {}/{}*L004_R2* > {}/{}_R2_001.fastq"    
```

<br>

#### Move the reads to a new folder named rawReads

Use mv command in bash.

<br>

#### Create sample names file

**Bash command:**

``` bash
$ ls rawReads | sed 's/_R1_001.fastq//g' | sed 's/_R2_001.fastq//g' | sort | uniq > sample_names.txt   
```

<br>

### Step 3: Quality control of Fastq files

To run fastqc in all the samples, use the script bin/qcReads.py.

**Script:** bin/qcReads.py

**Arguments:**

|                                                                                                                 |
|:----------------------------------------------------------------------------------------------------------------|
| -h, --help show this help message and exit                                                                      |
| -v, --version show program's version number and exit                                                            |
| --analysis\_info\_file ANALYSIS\_INFO\_FILE. Text file with details of the analysis. Default=analysis\_info.txt |
| --in\_dir IN\_DIR. Path to folder containing fastq files. Default=rawReads/                                     |
| --out\_dir OUT\_DIR. Path to out put folder. Default=rawReads/                                                  |
| --sample\_names\_file SAMPLE\_NAMES\_FILE. Text file with sample names to run. Default=sample\_names.txt        |

**Output:**
- Fastqc files and folders will be created for each sample in the rawReads/ folder.
- Log file with version of the software used in rawReads/ folder.

**Command example:**

``` bash
$ python bin/qcReads.py --analysis_info_file analysis_info.txt --in_dir rawReads/ --out_dir rawReads/ --sample_names_file sample_names.txt
```

<br>

### Step 4: Table and plot of number of reads per sample

**Script:** bin/indexQC.R

**Requires:** ggplot2 and reshape R libraries

**Input:** The input directory should have the *fastqc\_data.txt* files of all samples of the project. The 7th line in the fastqc\_data.txt files have the total number of sequences. This is the information that the script collects.

**Command:**

``` bash
$ /usr/bin/Rscript bin/indexQX.R <input dir> <output dir>
```

**Command example:**

``` bash
$ /usr/bin/Rscript bin/indexQX.R rawReads/ Report/figure/data/
```

**Output:** The output directory will have a csv file with the Total number of Sequences and a bar plot with the proportion of reads represented by each sample from the total number of reads from the run.

<br>

### Step 5: FastQC plots

**main script:** bin/fastqc\_tables\_and\_plots.py

**sub scripts:** bin/create\_fastqcPlots\_allSamples.R; bin/create\_fastqcPlots\_perSample.R; bin/create\_fastqcTables.py

**Arguments main script:**

|                                                                                                  |
|:-------------------------------------------------------------------------------------------------|
| -h, --help show this help message and exit                                                       |
| -v, --version show program's version number and exit                                             |
| --in\_dir IN\_DIR Path to folder containing fastq files. Default=rawReads/                       |
| --out\_dir\_report OUT\_DIR\_REPORT Path to out put folder. Default=Report/figure                |
| --suffix\_name SUFFIX\_NAME Suffix to optionally put to the output name. Default=''              |
| --sample\_names\_file SAMPLE\_NAMES\_FILE Text file with sample names. Default=sample\_names.txt |
| --plot\_device PLOT\_DEVICE Specify the format of the plot output. Default=png                   |

**Command example:**

``` bash
$ python bin/fastqc_tables_and_plots.py --in_dir rawReads/ --out_dir_report Report/figure/rawQC/ --suffix_name _raw --sample_names_file sample_names.txt 
```

<br>

### Step 6: Trim low quality bases and adapters

**Check the trimgalore\_params in the analysis info file. Most are default, but paired end or single end will need to be changed depending on the project.** Manual: <https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md>.
- add in '--paired' into trimgalore params in analysis info file, as well as changing readType to pairedEnd.

**main script:** bin/trimmingReads.py

**other scripts:** bin/trimming\_summary.R

**analysis infofile:** Define the trimgalore parameters that you want to pass to trim galore and the number of cores.

**Lexogen projects analysis infofile:** For lexogen projects, the argument *--clip\_R1 12* needs to be added in the analysis\_info.txt file (trimgalore\_params line) to remove the first 12 bases as recommended.

**Arguments:**

|                                                                                                                 |
|:----------------------------------------------------------------------------------------------------------------|
| -h, --help show this help message and exit                                                                      |
| -v, --version show program's version number and exit                                                            |
| --analysis\_info\_file ANALYSIS\_INFO\_FILE. Text file with details of the analysis. Default=analysis\_info.txt |
| --in\_dir IN\_DIR. Path to folder containing fastq files. Default=rawReads/                                     |
| --out\_dir OUT\_DIR. Path to out put folder. Default=trimmedReads/                                              |
| --sample\_names\_file SAMPLE\_NAMES\_FILE. Text file with sample names to run. Default=sample\_names.txt        |

**Command example:**

``` bash
$ python bin/trimmingReads.py
```

<br>

#### Step 6 (Optional for Lexogen): Trim polyA sequences

PolyA sequences may also affect the mapping of reads if the alignment tool used involves alignment of reads from end-to-end. The tool used for mapping in this analysis is STAR, as explained in the following section.

STAR does local alignment and clipping of 3’-end of reads in case that such an alignment gives a better score. This means that the majority of polyA sequences are clipped from reads during the alignment step. For this reason it is not strictly necessary to trim polyA sequences from the reads.

To remove polyA sequences I have used prinseq (argument -trim\_right 8, but may change).

<br>

### Step 7: Trimming summary

The script bin/trimming\_summary creates a table with the number of raw sequences, the number of sequences after trimming, and the percentage of sequences removed after trimming.

It expects an input directory of **raw data** with a fastqc files for each sample (fastqc\_data.txt; these should be in subdirectories).

It also expects an input directory of **trimmed data** with a fastqc files for each sample (fastqc\_data.txt; these should be in subdirectories).

**Arguments:**

``` bash
/usr/bin/Rscript bin/trimming_summary.R <raw_data_indir> <trimmed_data_indir> <outdir>   
```

**Command example:**

``` bash
/usr/bin/Rscript bin/trimming_summary.R rawReads/ trimmedReads/ Report/figure/data/
```

**Output:** A csv table in the output directory.

<br>

### Step 8: Trimming QC plots

**Make sure the arguments, specially *in\_dir* and *out\_dir* arguments, correspond to trimmed reads folders**

**main script:** bin/fastqc\_tables\_and\_plots.py

**other scripts:** bin/create\_fastqcPlots\_allSamples.R; bin/create\_fastqcPlots\_perSample.R; bin/create\_fastqcTables.py

**Arguments main script:**

|                                                                                                  |
|:-------------------------------------------------------------------------------------------------|
| -h, --help show this help message and exit                                                       |
| -v, --version show program's version number and exit                                             |
| --in\_dir IN\_DIR Path to folder containing fastq files. Default=rawReads/                       |
| --out\_dir\_report OUT\_DIR\_REPORT Path to out put folder. Default=Report/figure                |
| --suffix\_name SUFFIX\_NAME Suffix to optionally put to the output name. Default=''              |
| --sample\_names\_file SAMPLE\_NAMES\_FILE Text file with sample names. Default=sample\_names.txt |
| --plot\_device PLOT\_DEVICE Specify the format of the plot output. Default=png                   |

**Command example:**

``` bash
$ python bin/fastqc_tables_and_plots.py --in_dir trimmedReads --out_dir_report Report/figure/trimmedQC/ --suffix_name _trimmed --sample_names_file sample_names.txt
```

<br>

### Step 9: Mapping

**main script:** bin/mappingReads.py

**software:**
1. samtools <http://www.htslib.org/doc/samtools.html>
2. STAR <https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf>

**analysis infofile:** The following fields need to be defined in the analysis info file:
1. Number of cores to be used. Normally 2 works fine.
2. Reference genome folder, which corresponds to the STAR genomeGenerate folder output.

**For small RNA seq projects:**
- consider changing the parameters in analysis info to --runThreadN 4; --outSAMtype BAM SortedByCoordinate; --readFilesCommand zcat; --outFilterMultimapNmax 20; --outFilterMismatchNoverLmax 0.03; --outFilterScoreMinOverLread 0; --outFilterMatchNminOverLread 0; --outFilterMatchNmin 16; --alignSJDBoverhangMin 1000; --alignIntronMax 1

**Comments:** a temp dir can be a local directory. Transfer of data between servers (i.e fs3 and cluster) is often very slow.

**Arguments main script:**

|                                                                                                                |
|:---------------------------------------------------------------------------------------------------------------|
| -h, --help show this help message and exit                                                                     |
| -v, --version show program's version number and exit                                                           |
| --in\_dir IN\_DIR Path to folder containing fastq files. Default=trimmedReads/                                 |
| --out\_dir OUT\_DIR Path to out put folder. Default=alignedReads/                                              |
| --temp\_dir TEMP\_DIR Path to a temp dir folder. Must be within /home/tmp/. Default=none                                               
| --sample\_names\_file SAMPLE\_NAMES\_FILE. Text file with sample names to map. Default=sample\_names\_info.txt |

**Command example:**

``` bash
$ python bin/mappingReads.py --temp_dir /home/tmp/temp_dir
```

<br>

### Step 10: Mapping QC

#### 10a: Mapping summary

**main script:** bin/mappingQC.py

**other scripts:** bin/mapping\_summary.R

**Input:** folder of mapped reads with STAR output files \*Log.final.out

**Command example:**

``` bash
python bin/mappingQC.py --run mapping_summary --in_dir alignedReads --out_dir_report Report/figure/mappingQC/
```

**Ouput:** csv file in the out\_dir\_report with percentage of uniquely mapped reads, reads mapped to multiple loci, unmapped reads, starting number of reads, and number of mapped and usable reads (uniquely mapped reads)

<br>

#### 10b: Gene body coverage

**main script:** bin/mappingQC.py

**main tools:** /usr/local/minicondaexport/bin/geneBody\_coverage.py

**Input:** folder of mapped reads with bam files.

**Bedfile:** bedFile\_10k specified in the analyis info file.

**Command example:**

``` bash
python bin/mappingQC.py --run gene_body_coverage --in_dir alignedReads --out_dir alignedReads/QC --out_dir_report Report/figure/mappingQC/
```

**Output:** pdf file with gene body coverage plot in the out\_dir\_report. *I sometimes re-do this plot using the 10KGenes.geneBodyCoverage.r, which is also an output from the script, and is stored in the out\_dir folder*

<br>

#### 10c: Genomic context of mapped reads.

**main script:** bin/mappingQC.py

**main tools:** /usr/local/picard/picard.jar CollectRnaSeqMetrics

**other scripts:** bin/read\_distribution\_genomic\_context.R, bin/mapping\_distribution.R

**Input:** folder of mapped reads with bam files.

**Other files:** ref flat file and ribosomal rna interval list, both specified in the analysis info file.

**Other options:** *strand*, specified in the analysis info file.

**Command example:**

``` bash
python bin/mappingQC.py --run picard_tools --in_dir alignedReads --out_dir alignedReads/QC --out_dir_report Report/figure/mappingQC/
```

**Output:** csv and txt files in out\_dir. read\_distribution\_genomic\_context png file and strand\_mappingQC csv file in out\_dir\_report.

<br>

### Step 11: Counting of reads

**main script:** bin/countingQC.py

**software:**
- samtools <http://www.htslib.org/doc/samtools.html>
- htseq-count <http://htseq.readthedocs.io/en/release_0.9.0/count.html>

**gtfFile:** gtf file specified in the analysis info file

**Other options:** *htseq\_params*, specified in the analysis info file

**For small RNA seq projects:**
- consider using the gtf file from mirBase. Make sure the chromosome names are in the same format. [Mice](ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3). [Human](ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3)
- **If the library is single end and maps to the forward strand, and you are using the gtf file from mirbase,** change the counting parameters to -a 0; -m union; -s yes; -t miRNA; -i Name.

**Command example:**

``` bash
python bin/countingReads.py --in_dir alignedReads --out_dir countedReads 
```

**Output:** \*\_count.txt\* files for each sample in the outdir with counts for each gene.

### Step 12: Counting QC

#### 12a: Counting QC part1

**Script:** bin/countingQC\_part1.R

**R packages required:** RColorBrewer, matrixStats, gplots, ggplot2, reshape, GenomicRanges, rtracklayer, Rsamtools, grid.

**Input:** folder with \_count.txt files

**Other files:** mapping\_summary.csv, output of mappinQC

**Command:**

``` bash
/usr/bin/Rscript bin/countingQC_part1.R <indir> <outdir> <mapping summary csv>  
```

**Command example:**

``` bash
/usr/bin/Rscript bin/countingQC_part1.R countedReads/ Report/figure/countingQC Report/figure/mappingQC/mapping_summary.csv
```

**Output:** files and plots in outdir

#### 12b: Counting QC part2

**Script:** bin/countingQC\_part2.R

**R packages required:** gplots, rtracklayer, reshape, ggplot2.

**Input:** folder with \_count.txt files

**Other files:** GTF file

**Command:**

``` bash
/usr/bin/Rscript bin/countingQC_part2.R <indir> <outdir> <gtf file>
```

**Command example:**

``` bash
/usr/bin/Rscript bin/countingQC_part2.R countedReads/ Report/figure/countingQC/ /mnt/cgs-fs3/Sequencing/Genome/Mouse/gtf/ensembl/grcm38/release-84/Mus_musculus.GRCm38.84.gtf
```

**Output:** files and plots in outdir

### Step 13: EdgeR and DESeq2

At the moment, edgeR is used to do the differential expression analysis. DESeq2 is only used to do PCA plots with **rlog** transformed data.

#### 13a: Create requiered files:

**Scripts:** bin/deseq2\_arguments.py, bin/edger\_arguments.py

**Commands:**

``` bash
python bin/edger_arguments.py   
python bin/deseq2_arguments.py
```

**Output:**
- **sample\_info.csv.** This file will have columns **SampleID** and **Group**. You can add a **Sibship** column for paired analysis. Fill this file manually.
- **comparisons.csv.** This file will have columns **baselineGroup** and **comparisonGroup**. The baselineGroup is normally the WT or control group. List the pairwise comparisons you want to make, using the corresponding group names.
- **deseq2\_arguments.txt.** Fill this file manually. Explanation and example below.
- **edger\_arguments.** Fill this file manually. Explanation and example below.

**deseq2 arguments txt file explanation:**

|                                                                                                            |
|:-----------------------------------------------------------------------------------------------------------|
| indir = \*&lt;folder with \_count.txt files&gt;\*                                                          |
| outdir = *&lt;output dir for results&gt;*                                                                  |
| sample\_info = *&lt;csv file with columns SampleID and Group (and Sibship for paired analysis)&gt;*        |
| comparisons = *&lt;csv file with columns baselineGroup and comparisonGroup with a list of comparisons&gt;* |
| design = *&lt;either pairedSamples OR non-pairedSamples&gt;*                                               |
| gtfFile = *&lt;gtf file with gene info. Can be copied from the analysis\_info file&gt;*                    |

**deseq2 arguments txt file example:**

|                                                                                                          |
|:---------------------------------------------------------------------------------------------------------|
| indir = countedReads/                                                                                    |
| outdir = deseq2/                                                                                         |
| sample\_info = sample\_info.csv                                                                          |
| comparisons = comparisons.csv                                                                            |
| design = non-pairedSamples                                                                               |
| gtfFile = /mnt/cgs-fs3/Sequencing/Genome/Mouse/gtf/ensembl/grcm38/release-84/Mus\_musculus.GRCm38.84.gtf |

**edger arguments txt file explanation:**

|                                                                                                                                |
|:-------------------------------------------------------------------------------------------------------------------------------|
| indir = \*&lt;folder with \_count.txt files&gt;\*                                                                              |
| outdir = *&lt;output dir for results&gt;*                                                                                      |
| sample\_info = *&lt;csv file with columns SampleID and Group (and Sibship for paired analysis)&gt;*                            |
| comparisons = *&lt;csv file with columns baselineGroup and comparisonGroup with a list of comparisons&gt;*                     |
| min.cpm = *&lt;Filtering threshold: miminum count per million value a gene should have to be kept&gt;*                         |
| min.nsamples = *&lt;Filtering threshold: miminum number of samples with cpm specified above a gene should have to be kept&gt;* |
| design = *&lt;either pairedSamples OR non-pairedSamples&gt;*                                                                   |
| gtfFile = *&lt;gtf file with gene info. Can be copied from the analysis\_info file&gt;*                                        |

**edger arguments txt file example:**

|                                                                                                          |
|:---------------------------------------------------------------------------------------------------------|
| indir = countedReads/                                                                                    |
| outdir = &lt;edger/                                                                                      |
| sample\_info = sample\_info.csv                                                                          |
| comparisons = comparisons.csv                                                                            |
| min.cpm = 1                                                                                              |
| min.nsamples = 1                                                                                         |
| design = non-pairedSamples                                                                               |
| gtfFile = /mnt/cgs-fs3/Sequencing/Genome/Mouse/gtf/ensembl/grcm38/release-84/Mus\_musculus.GRCm38.84.gtf |

#### 13b: DESeq2 PCA plot

**main script:** bin/deseq2.R

**other scripts:** bin/deseq2\_functions.R

**R packages:** DESeq2, gplots, rtracklayer, ggplot2, gridExtra, matrixStats, dplyr.

The following script uses DESeq2 package to transform raw counts to rlog values.

The input directory, output directory, and other arguments need to be specified in the deseq2\_arguments.txt.
**Command:**

``` bash
/usr/bin/Rscript bin/deseq2.R deseq2_arguments.txt   
```

Output: PCA and heatmap plots in the outdir specified.

#### 13c: Differential gene expression.

**main script:** bin/edger.R

**Other scripts:** bin/edger\_functions.R

**R packages:** edgeR, gplots, rtracklayer.

**Command:**

``` bash
/usr/bin/Rscript bin/edger.R edger_arguments.txt
```

**Output:** results in the output dir.
