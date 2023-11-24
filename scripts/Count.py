#!/usr/bin/env python3
# Author: zhanglele

import argparse
import os
import subprocess
import gzip
import multiprocessing as mp
import contextlib
import shutil
from datetime import datetime

@contextlib.contextmanager
def cd(cd_path):
    saved_path = os.getcwd()
    os.chdir(cd_path)
    yield
    os.chdir(saved_path)



def create_STAR_index(args):
    cmd = 'STAR'+' --runMode genomeGenerate' \
        +' --runThreadN '+str(args.nthread) \
        +' --genomeDir '+args.STAR_index \
        +' --genomeFastaFiles '+args.ref_genome \
        +' --sjdbGTFfile '+args.annotation \
        +' --sjdbOverhang 149' \
    
    return cmd

def STAR(args):
    cmd = 'STAR'+' --runThreadN '+str(args.nthread) \
        +' --genomeDir '+args.STAR_index \
        +' --outFileNamePrefix '+align_dir+'/'+args.prefix+'_' \
        +' --readFilesIn '+args.fastq1+' '+args.fastq2 \
        +' --readFilesCommand zcat' \
        +' --outSAMtype BAM SortedByCoordinate' \
        +' --outSAMattributes NH HI NM MD XS AS' \
        +' --outFilterMultimapNmax 500' 


def STARsolo(args):
    cmd = 'STAR'+' --soloType CB_UMI_Simple' \
        +' --runThreadN '+str(args.nthread) \
		+' --genomeDir '+args.STAR_index \
        +' --readFilesIn '+args.fastq2+args.fastq1 \
		+' --outFileNamePrefix'+align_dir+'/'+args.prefix+'_' \
		+' --outFilterMultimapNmax'+args.multimap \
		+' --outSAMattributes NH HI NM MD XS AS' \
        +' --readFilesCommand zcat'\
        +' --clipAdapterType CellRanger4'\
        +' --outFilterScoreMin 30'\
        +' --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts' \
        +' --soloUMIfiltering MultiGeneUMI_CR'\
        +' --soloUMIdedup 1MM_CR '\
        +' --soloMultiMappers EM'\
        +' --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM'
    if args.nthreadsort:
        cmd += ' --outBAMsortingThreadN '+str(args.nthreadsort)
    if args.nRAMsort:
        cmd += ' --limitBAMsortRAM '+str(args.nRAMsort)
    if args.version=="V3":
        cmd += ' --soloUMIlen 12  --soloCBwhitelist '+ script_dir+'/data/3M-february-2018.txt.gz'
    if args.version=="V2":
        cmd += ' --soloCBwhitelist '+ script_dir+'/data/737K-august-2016.txt'
    return cmd

parser = argparse.ArgumentParser(description='SERVE: pipeline for detecting expressed ERVs')
parser.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser.add_argument('-v', '--version',default="V2", help='specific version of the 10X chemistry (default:V3)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format')
parser.add_argument('-S', '--STAR_index', default='./STAR_index', help='Path to the directory where STAR index generated (default: STAR_index)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run SERVE (default: 1)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')


args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))
align_dir='./align'


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running Counting on {0:d} threads.'.format(args.nthread), flush=True)


with cd(args.output_dir):
    if not os.path.exists(args.STAR_index):
        os.makedirs(args.STAR_index)
        if not args.annotation:
            print('ERROR: Lack annotation file (--annotation)')
        if not args.ref_genome:
            print('ERROR: Lack reference genome (--ref_genome)')

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to create STAR index.', flush=True)

        cmd = create_STAR_index(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./STAR_index/).', flush=True)


    if not os.path.exists(align_dir):    
        os.makedirs(align_dir)

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to Count Singlecell RNASEQ.', flush=True)
    
    if args.seqtype=="Smart"|"bulk":
        cmd =STAR(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] STAR is done , now FeatureCounts is starting', flush=True)
        
        cmd = 'featureCounts -O --fracOverlap 0.1  -p -B -M --fraction -T'+ args.nthread +' -a '+ args.annotation + ' -o  sample_featurecount -R BAM '+align_dir+'/'+args.prefix+'_Aligned.sortedByCoord.out.bam' 
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Reformating bam files', flush=True)
        
        cmd = 'samtools view -x BX -x NH -x AS -x nM -x HI -x IH -x NM -x uT -x MD -x jM -x jI -x XN -x XS -x XS'\
            + align_dir+'/'+args.prefix+'_Aligned.sortedByCoord.out.bam.featureCounts.bam | grep -e "XT:Z" > sample_featurecount.sam.txt'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        
        cmd = 'Rscript '+script_dir+'/redistribute_multiple_aligned_reads.r  -f sample_featurecount.sam.txt -r' + args.annotations + ' -n '+ args.prefix + ' -s 50 -m 1 -p '+args.nthread
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 
        
        cmd = 'rm sample_featurecount.sam.txt'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 

        cmd = 'rm'+ args.prefix+'_Aligned.sortedByCoord.out.bam.featureCounts.bam'
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')     
        
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Counting is done', flush=True)

    else:
        cmd = STARsolo(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')


   


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Counting is done', flush=True)