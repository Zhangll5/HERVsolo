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


def Taco_merge(args):
    cmd = 'taco_run'+' '+args.input_gtf_list \
		+' -p '+str(args.nthread) \
		+' --gtf-expr-attr '+'TPM' \
		+' -o '+merge_dir
    return cmd

def Intersect(args):
    cmd = 'bedtools intersect'\
        + ' -a '+args.erv_bed \
        + ' -b '+ERV_gtf \
        + ' -f 1.00 ' \
        + ' -wb '\
        + ' > ' + merge_bed;\
        +' bedtools intersect '\
        + ' -a '+args.erv_bed \
        + ' -b '+ERV_gtf \
        + " -wo |awk '$NF > 200' >>" + merge_bed
    return cmd

def Selectgenes(args):
    cmd = "awk '{print $28}'"+merge_bed+ '| uniq >' +select_genes
    return cmd
def Selecttrans(args):
    cmd = "awk '{print $26}'"+merge_bed+ '| uniq >' +select_trans
    return cmd

def GtfReturnERV(args):
    cmd = "grep -v " + select_genes +ERV_gtf +'>' +Final_ERV_gtf
    return cmd
def GtfReturnFulllength(args):
    cmd = "grep -v " + select_trans +ERV_gtf +'>' +Final_ERV_gtf    
    return cmd
def MergeGTF(args):
    cmd = 'cat '\
        + Final_ERV_gtf \
        + args.annotation\
        + '>'\
        + Merge_gtf
      


parser = argparse.ArgumentParser(description='SERVE_merge: merge expressed ERVs')
parser.add_argument('-i', '--input_gtf_list', help='A text file with a list of SERVE GTF files (required)')
parser.add_argument('-p', '--prefix', default='SERVE', help='Prefix for output file name (default: SERVE)')
parser.add_argument('-e', '--erv_bed', help='ERV position in BED format (required)')
parser.add_argument('-r', '--ref_genome', help='Reference genome in FASTA format (required)')
parser.add_argument('-a', '--annotation', help='Genome annotation in GTF format')
parser.add_argument('-t', '--nthread', type=int, default=1, help='Number of threads to run SERVE_merge (default: 1)')
parser.add_argument('-l', '--length', default=200, help='Minimum ERV overlaplength (bp) (default: 200)')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory (default: .)')

args = parser.parse_args()

merge_dir= args.output_dir+'/ERV_merge'

###assem_gtf = args.output_dir+'/'+args.prefix+'_assem.gtf'
ERV_gtf = args.output_dir+'/'+args.prefix+'_ERV_merge.gtf'
Full_gtf = args.output_dir+'/'+args.prefix+'_FullLength_merge.gtf'
merge_bed = args.output_dir+'/'+args.prefix+'_merge.bed'
select_trans = args.output_dir+'/'+args.prefix+'_select_trans'
select_genes = args.output_dir+'/'+args.prefix+'_select_genes'
Final_ERV_gtf = args.output_dir+'/'+args.prefix+'_Final_ERV.gtf'
Merge_gtf = args.output_dir+'/'+args.prefix+'_merge.gtf'

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Running SERVE_merge on {0:d} threads.'.format(args.nthread), flush=True)


with cd(args.output_dir):
    if args.Full_length==False:
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to merge ERV gtfs.', flush=True)
        
        cmd = Taco_merge(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        shutil.copyfile('./ERV_merge/assembly.gtf',ERV_gtf)
        shutil.rmtree(merge_dir)
        
        cmd = Intersect(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 
    
        cmd = Selectgenes(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 

        cmd = GtfReturnERV(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 
        
        cmd=MergeGTF(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 

        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)

    if args.Full_length==True:
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Start to merge Full_length_transcriptome gtfs.', flush=True)
        cmd = Taco_merge(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash')
        shutil.copyfile('./ERV_merge/assembly.gtf',Full_gtf)
        shutil.rmtree(merge_dir)
        print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] Finish.', flush=True)
   
        cmd = Intersect(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 
    
        cmd = Selecttrans(args)
        subprocess.check_call(cmd, shell=True, executable='/bin/bash') 

        cmd = GtfReturnFulllength(args)


print('['+datetime.now().strftime("%b %d %H:%M:%S")+'] ERV merge is done. Now you can quantify expressed ERVs with SERVE_quant.', flush=True)
