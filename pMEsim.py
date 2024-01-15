#!/usr/bin/python3

import argparse
import sys
import subprocess
import os
import glob
import datetime
from pathlib import Path

# Create a timestamp
def stamp():
    t=str(datetime.datetime.now().strftime("[%Y-%M-%d %H:%M:%S]"))
    return(t)

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
parser.add_argument('-d', '--out-dir', type = str, metavar='STR', help='directory to store the simulation\'s output -- will create if doesn\'t exist', dest = 'out_dir', default = 'out')
## 01-specific arguments
parser.add_argument('-f', '--target_fasta1', type = str, metavar='STR', help='fasta file for the target chromosome', dest = 'target_fasta1', default = '../simData/hg38.chr22.fa') # first step the target is a the ref genome chromosome 22
parser.add_argument('-o', '--out1', type=str, metavar='STR', help='output file prefix (will be <prefix>.vcf)', dest = 'out_prefix1', default = 'simRef')
## 02-specific arguments
parser.add_argument('-F', '--target_fasta2', type = str, metavar='STR', help='fasta file for the target chromosome', dest = 'target_fasta2', default = 'simRef.simseq.genome.fa') # now the target is the refSim genome (made from chr22)
parser.add_argument('-B', '--bed2', type=str, metavar='STR', help='bed file with simRef TEs to delete', dest = 'bed2_in', default = 'ins2del.bed') # new bed to take the DEL coordinates
parser.add_argument('-O', '--out2', type=str, metavar='STR', help='output file prefix (will be <prefix>.vcf)', dest = 'out_prefix2', default = 'simAlt')
# parse the arguments
args = parser.parse_args()

# we need to clean the working directory
if args.verbose:
    print("cleaning up last session...")
for filename in glob.glob('./simRef*'):
    try:
        os.remove(filename)
    except OSError:
        pass
for filename in glob.glob('./simAlt*'):
    try:
        os.remove(filename)
    except OSError:
        pass
try:
    os.remove('ins2del.bed')
except OSError:
    pass

### according to verbose or not, we have two ways per command:
if args.verbose:
    command01 = str("python3 " + "01_makeIns2Del.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta1) + " -o " + str(args.out_prefix1) + " -V")
    command02 = str("python3 " + "02_makeAltGenome.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta2) + " -B " + str(args.bed2_in) + " -o " + str(args.out_prefix2) + " -V")
else:
    command01 = str("python3 " + "01_makeIns2Del.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta1) + " -o " + str(args.out_prefix1))
    command02 = str("python3 " + "02_makeAltGenome.py " + "-v " + str(args.nb_var) + " -r " + str(args.ratio) + " -g " + str(args.ref_genome) + " -t " + str(args.target_chrom) + " -b " + str(args.bed_in) + " -R " + str(args.tsdrange) + " -f " + str(args.target_fasta2) + " -B " + str(args.bed2_in) + " -o " + str(args.out_prefix2))

### actually run the scripts
if args.verbose:
    print("running script 1...")
    print(command01)
script1 = subprocess.Popen(str(command01), 
    shell = True#, 
    #stdout=subprocess.DEVNULL,
    #stderr=subprocess.STDOUT
    )
script1.wait()
if args.verbose:
    print("running script 2...")
    print(command02)
script2 = subprocess.Popen(str(command02), 
    shell = True#, 
    #stdout=subprocess.DEVNULL,
    #stderr=subprocess.STDOUT
    )
script2.wait()

### organize the outputs
out_d = args.out_dir
# create dir, equivalent to bash `mkdir -P`
try:
    Path(out_d).mkdir(parents=True, exist_ok=True)
except OSError:
    print('[ERROR!]' + stamp() + 'output dir: ' + out_d + ' cannot be created')
    exit()
# move the files
for filename in glob.glob('simRef*'):
    os.replace(filename, out_d + "/" + filename)
for filename in glob.glob('simAlt*'):
    os.replace(filename, out_d + "/" + filename)
os.replace('ins2del.bed', out_d + '/ins2del.bed')

