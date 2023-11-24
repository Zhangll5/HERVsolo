#!/usr/bin/env python3
# Author: Zhanglele

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
        + ' --runThreadN '+str(args.nthread) \
        + ' --genomeDir '+args.STAR_index \
        + ' --genomeFastaFiles '+args.ref_genome \
        + ' --sjdbGTFfile '+args.annotation \
        + ' --sjdbOverhang 149' \
        + ' --genomeSAindexNbases '+str(args.genomeSAindexNbases)
    
    return cmd


def STAR(args):
    cmd = 'STAR'+' --runThreadN '+str(args.nthread) \
        + ' --genomeDir '+args.STAR_index \
        + ' --outFileNamePrefix '+align_dir+'/'+args.prefix+'_' \
        + ' --readFilesIn '+args.fastq1+' '+args.fastq2 \
        + ' --readFilesCommand zcat' \
        + ' --outSAMtype BAM SortedByCoordinate' \
        + ' --outSAMattributes NH HI AS nM NM' \
        + ' --outFilterMultimapNmax 500' 
    if args.nthreadsort:
        cmd += ' --outBAMsortingThreadN '+str(args.nthreadsort)
    if args.nRAMsort:
        cmd += ' --limitBAMsortRAM '+str(args.nRAMsort)
    
    return cmd


def create_BAM_index(args):
    cmd = 'sambamba index'+' -t '+str(args.nthread) \
                          +' '+align_bam
    
    return cmd

def Extract(args):
    cmd = 'sambamba view'+' -t '+str(args.nthread) \
                         +' -L '+args.erv_bed \
                         +' -f bam' \
                         +' -o '+ERV_bam \
                         +' '+align_bam

    


def Assemble(args):
    cmd = 'stringtie'+ERV_bam \
                     + ' -j ' +args.junction     \
                     + ' -s '+args.single \
                     + ' -c '+args.multi \
                     + ' -f '+args.isoform \
                     + ' -p '+str(args.nthread) \
                     + ' -o '+ERV_gtf
    if args.stranded_type:
        cmd += ' --'+args.stranded_type
    
    return cmd

def AssembleFulllength(args):
    cmd = 'stringtie'+align_bam \
                     + ' -j ' +args.junction     \
                     + ' -s '+args.single \
                     + ' -c '+args.multi \
                     + ' -f '+args.isoform \
                     + ' -p '+str(args.nthread) \
                     + ' -o '+ALL_gtf
    
    return cmd




parser = argparse.ArgumentParser(description='SERVE: pipeline for detecting expressed ERVs')
parser.add_argument('-f', '--Full_length',default=False, help='Read1 in FASTQ format (required)')
parser.add_argument('-fq1', '--fastq1', help='Read1 in FASTQ format (required)')
parser.add_argument('-fq2', '--fastq2', help='Read1 in FASTQ format (required)')
parser.add_argument('-e', '--erv_bed', help='ERV position in BED format (required)')
parser.add_argument('-p', '--prefix', default='SERVE', help='Prefix for output file name (default: SERVE)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format')
parser.add_argument('-S', '--STAR_index', default='./STAR_index', help='Path to the directory where STAR index generated (default: STAR_index)')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run SERVE (default: 1)')
parser.add_argument('--nthreadsort', type=int, help='Number of threads for BAM sorting')
parser.add_argument('--nRAMsort', type=int, help='Maximum available RAM (bytes) for sorting BAM.')
parser.add_argument('-j', '--junction', default='2',help='minimum junction coverage (default: 2)')
parser.add_argument('--single-exon', default='5',help='minimum reads per bp coverage to consider for single-exon transcript (default: 5)')
parser.add_argument('--multi-exon', default='2',help='minimum reads per bp coverage to consider for single-exon transcript (default: 2)')
parser.add_argument('--filter', default='0.05',help='filter transcripts with a low proportion of expression  (default: 0.05)')
parser.add_argument('--stranded_type', default=None, help='Strand-specific RNA-seq read orientation: rf or fr (default: None)')
parser.add_argument('-m', '--nRAMassem', default='10G', help='Maximum available RAM (Gb) for assembly (default: 10G)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')

args = parser.parse_args()
script_dir = os.path.abspath(os.path.dirname(__file__))

align_dir = './align'
assem_dir = './assem'


align_bam = align_dir+'/'+args.prefix+'_Aligned.sortedByCoord.out.bam'
ERV_bam = align_dir+'/'+args.prefix+'_ERV.bam'
ERV_gtf = assem_dir+'/'+args.prefix+'_ERV.gtf'
ALL_gtf = assem_dir+'/'+args.prefix+'_ALL.gtf'


if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S") +
      '] Running SERVE on {0:d} threads.'.format(args.nthread), flush=True)


with cd(args.output_dir):
    if not os.path.exists(args.STAR_index):
        os.makedirs(args.STAR_index)
        if not args.annotation:
            print('ERROR: Lack annotation file (--annotation)')
        if not args.ref_genome:
            print('ERROR: Lack reference genome (--ref_genome)')

        print('['+datetime.now().strftime("%b %d %H:%M:%S") +
              '] Start to create STAR index.', flush=True)

        cmd = create_STAR_index(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

        print('['+datetime.now().strftime("%b %d %H:%M:%S") +
              '] Finish (In ./STAR_index/).', flush=True)

    if not os.path.exists(align_dir):
        os.makedirs(align_dir)

    print('['+datetime.now().strftime("%b %d %H:%M:%S") +
          '] Start to align RNA-seq reads to reference genome.', flush=True)

    cmd = STAR(args)
    subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    print('['+datetime.now().strftime("%b %d %H:%M:%S") +
          '] Finish (In ./align/).', flush=True)
    
    if args.Full_length==False:

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to extract ERV reads.', flush=True)
        cmd = create_BAM_index(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

        cmd = Extract(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish (In ./1_align/).', flush=True)



        if not os.path.exists(assem_dir):
            os.makedirs(assem_dir)

        
        print('['+datetime.now().strftime("%b %d %H:%M:%S") +
            '] Start to de novo assemble ERVs.', flush=True)

        cmd = Assemble(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

    if args.Full_length==True:   
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to de novo assemble Full_length_transcriptome .', flush=True)
        
        cmd = AssembleFulllength(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')

print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Assembling is down. Now you can do GTF merge HERV_RNA_merge', flush=True)