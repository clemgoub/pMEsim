#!/usr/bin/python3

import argparse
import sys
import subprocess

# arguments for the wrapper (this script)
## common arguments
parser = argparse.ArgumentParser(description='Create a simulated VCF file with random TE insertions or deletions with TSD')
parser.add_argument('-v', '--variants', type = int, metavar='N', help='the total number of variants to simulate', dest = 'nb_var', default = 100)
parser.add_argument('-r', '--ratio', type = float, metavar='[0-1]', help='insertion/deletion ratio', dest = 'ratio', default = .7)
parser.add_argument('-g', '--genome', type = str, metavar='STR', help='reference genomes (fa/fa.gz)', dest = 'ref_genome', default = '../simData/hg38.UCSC.1to22XY.fa.gz')
parser.add_argument('-t', '--target_chr', type = str, metavar='STR', help='name of target chromosome where insertion/deletion will be performed', dest = 'target_chrom', default = 'chr22')
parser.add_argument('-b', '--bed', type=str, metavar='STR', help='bed file with reference insertions to use', dest = 'bed_in', default = '../simData/hg38.AluY.L1HSPA2.SVA_EF.bed')
parser.add_argument('-R', '--tsdrange', type = str, metavar='min,max', help='TSD length range', dest = 'tsdrange', default = '4,22' )
parser.add_argument('-V', '--verbose', action="store_true", help="increase output verbosity")
## 01-specific arguments
parser.add_argument('-f', '--target_fasta1', type = str, metavar='STR', help='fasta file for the target chromosome', dest = 'target_fasta1', default = '../simData/hg38.chr22.fa') # first step the target is a the ref genome chromosome 22
parser.add_argument('-o', '--out1', type=str, metavar='STR', help='output file prefix (will be <prefix>.vcf)', dest = 'out_prefix1', default = 'simRef')
## 02-specific arguments
parser.add_argument('-F', '--target_fasta2', type = str, metavar='STR', help='fasta file for the target chromosome', dest = 'target_fasta2', default = 'simRef.simseq.genome.fa') # now the target is the refSim genome (made from chr22)
parser.add_argument('-B', '--bed2', type=str, metavar='STR', help='bed file with simRef TEs to delete', dest = 'bed2_in', default = 'ins2del.bed') # new bed to take the DEL coordinates
parser.add_argument('-O', '--out2', type=str, metavar='STR', help='output file prefix (will be <prefix>.vcf)', dest = 'out_prefix2', default = 'simAlt')
# parse the arguments
args = parser.parse_args()

### according to verbose or not, we have two ways per command:
if args.verbose:
    command01 = str("python3 " + "01_makeIns2Del.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta1) + " -o " + str(args.out_prefix1) + " -V")
    command02 = str("python3 " + "02_makeAltGenome.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta2) + " -B " + str(args.bed2_in) + " -o " + str(args.out_prefix2) + " -V")
else:
    command01 = str("python3 " + "01_makeIns2Del.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta1) + " -o " + str(args.out_prefix1))
    command02 = str("python3 " + "02_makeAltGenome.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta2) + " -B " + str(args.bed2_in) + " -o " + str(args.out_prefix2))

### to change
if args.verbose:
    print("running script 1...")
    print(command01)
subprocess.call(str(command01), 
        shell = True, 
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)
if args.verbose:
    print("running script 2...")
    print(command02)
subprocess.call(str(command02), 
        shell = True, 
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)