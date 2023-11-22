#!/usr/bin/env python3
# Author: Janky

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

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to align RNA-seq reads to reference genome.', flush=True)

    cmd = STARsolo(args)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./align/).', flush=True)
    

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to extract ERV reads.', flush=True)

    cmd = create_BAM_index(args)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    cmd = Extract(args)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./1_align/).', flush=True)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] ERV identification is done. Now you can do ERV merge with SERVE_merge (Recommond).', flush=True)