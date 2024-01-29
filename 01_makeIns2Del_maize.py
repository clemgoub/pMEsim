#!/usr/bin/python3

import vcfpy # must be >= 0.13.7 or .8; not working with 0.12.x
import argparse
import pandas
import pysam
import random
import sys
import statistics
import subprocess
from alive_progress import alive_bar
import time
import datetime

# Define a custom argument type for a list of integers
def list_of_ints(arg):
    return list(map(int, arg.split(',')))
# Create a timestamp
def stamp():
    t=str(datetime.datetime.now().strftime("[%Y-%M-%d %H:%M:%S]"))
    return(t)

# create the parser for the user input
parser = argparse.ArgumentParser(description='Create a simulated VCF file with random TE insertions or deletions -- maize version')
parser.add_argument('-v', '--variants', type = int, metavar='N', help='the total number of variants to simulate', dest = 'nb_var', default = 100)
# the ratio ins/del is unknown a priori, so we go with 50/50
parser.add_argument('-r', '--ratio', type = float, metavar='[0-1]', help='insertion/deletion ratio', dest = 'ratio', default = .5) 
# use Mo17 T2T refernce
parser.add_argument('-g', '--genome', type = str, metavar='STR', help='reference genomes (fa/fa.gz)', dest = 'ref_genome', default = '../simData/Zm-Mo17-REFERENCE-CAU-2.0.fa')
# target is chromosome 10
parser.add_argument('-t', '--target_chr', type = str, metavar='STR', help='name of target chromosome where insertion/deletion will be performed', dest = 'target_chrom', default = 'chr10')
parser.add_argument('-f', '--target_fasta', type = str, metavar='STR', help='fasta file for the target chromosome', dest = 'target_fasta', default = '../simData/Zm-Mo17-REFERENCE-CAU-2.0__chr10__.fa')
# this bed is parsed so each entry is a potential TE to ins2del
parser.add_argument('-b', '--bed', type=str, metavar='STR', help='bed file with reference insertions to use', dest = 'bed_in', default = '../simData/Zm-Mo17-REFERENCE-CAU-2.0.fa.RMout_CopGypDTA4sim.bed')
# remove TSD range as it will be custom for each super-family
#parser.add_argument('-R', '--tsdrange', type = list_of_ints, metavar='min,max', help='TSD length range', dest = 'tsdrange', default = '4,22' )
parser.add_argument('-o', '--out', type=str, metavar='STR', help='output file prefix (will be <prefix>.vcf)', dest = 'out_prefix', default = 'simRef')
parser.add_argument('-V', '--verbose', action="store_true", help="increase output verbosity")
# parse the arguments
args = parser.parse_args()

# store the input arguments
ratio = args.ratio
# make the number of variants even
if (args.nb_var % 2) == 0:
    nb_var = args.nb_var
else:
    nb_var = args.nb_var + 1
    print("[WARNING]" +  stamp() + " you have asked for an odd number of variants  " + str(nb_var-1) + ",  " + str(nb_var) + " will be used instead.")
ref_genome = args.ref_genome
target_fasta = args.target_fasta
target_chrom = args.target_chrom
out_prefix = args.out_prefix
# # force the Alu, L1, SVA ratio to realistic values
# te_props = [0.7,0.2,0.1]
# calculate ins/del number based on input
ins_nb = round(ratio*nb_var)
del_nb = nb_var - ins_nb #round(nb_var-(ratio*nb_var))
# load the reference genome
fasta = pysam.FastaFile(ref_genome)
# read and store the input bed file; there is one more column that humans for "super-family"
bed_in_pre = pandas.read_csv(args.bed_in, sep = '\t',
                          names = ['chrom', 'start', 'end', 'TE', 'score', 'strand', 'super'])
# make sure the chromosomes of the bed file are in the reference genome
bed_in = bed_in_pre[bed_in_pre['chrom'].isin(fasta.references)]

# THIS IS UNNECESSARY AS LONG AS YOU KEEP VARIANT NUMBER REASONABLE
# # check input length to make sure we don't sample more than there is; here we only want to make insertions that will be later deleted.
# #in_del_nb = len(bed_in[bed_in['chrom'].isin([target_chrom])])
# in_ins_nb = len(bed_in[-bed_in['chrom'].isin([target_chrom])])
# # here we want to use the del number from the ratio, even though we make "insertions" -- this is because we need to insert the right number of TE to be later deleted.
# if in_ins_nb < ins_nb:
#     print('*******')
#     print('[ERROR!]' + stamp() + ' not enough TEs in the bed file to insert in the target chromosome!')
#     print('*******')
#     parser.print_help(sys.stderr)
#     sys.exit(1)

# sample the insertions (future deletions) according to predefined proportions (50/50 LTR/DTA)
i_LTR = bed_in[-bed_in['chrom'].isin([target_chrom]) & (bed_in.super.str.match('Copia') | bed_in.super.str.match('Gypsy'))].sample(round(max(del_nb*0.5,1)))
i_DTA = bed_in[-bed_in['chrom'].isin([target_chrom]) & bed_in.super.str.match('DTA')].sample(round(max(del_nb*0.5,1)))
# concatenate
repmask_subset = pandas.concat([i_LTR, i_DTA], axis=0).reset_index().sort_values(by=['chrom', 'start'])
if args.verbose:
    print(repmask_subset)
# read and store the header
reader = vcfpy.Reader.from_path('header_maize.vcf')
# set the name of the sample in the genotype column for the header we just stored
reader.header.samples = vcfpy.SamplesInfos(['sim'])
# initialize the output file. For now it only includes the header with the modified sample name
writer = vcfpy.Writer.from_path(out_prefix + '.vcf', reader.header)


#########################
# Simulate TE ins / del #
#########################
# initialize the current chromosome and its sequence
current_chrom = None
current_chrom_seq = None
# define which chromosome will receive the annotations
target_chrom_seq = fasta[target_chrom]
# create an empty list to store the positions of the simulated ins/del (1-based)
sim_pos = []
# loop over each line of the bed file
print('[info]' + stamp() + ' creating new reference pMEI to delete in the next round...')
if args.verbose:
    print(repmask_subset)
with alive_bar(len(repmask_subset.index), bar = 'circles', spinner = 'classic') as bar:
    for index, repeat in repmask_subset.iterrows():
        chrom = repeat['chrom']
        start = repeat['start']
        end = repeat['end']
        te = repeat['TE']
        sup = repeat['super']
        # define TSD range according to superfamily
        # tsdrange = []
        if sup == "DTA":
            tsdrange = [8,8]
        else: # LTR
            tsdrange = [4,6]

        # update in memory chromosome sequence
        if current_chrom != chrom:
            #print("Switching to " + chrom)
            current_chrom = chrom
            current_chrom_seq = fasta[current_chrom]
        if args.verbose:
            print("creating insertions from: " + chrom + '...')
        # pick a random position within the target chromosome, avoid N
        ref_sequence = 'N'
        while ref_sequence == 'N':
            rnd_pos = random.randint(1,fasta.get_reference_length(target_chrom))
            # get the 1 base at the site, will be the reference sequence
            # we assume rnd_pos is 1-based (VCF). Thus we -1 it to get it in python
            ref_sequence = target_chrom_seq[rnd_pos - 1] # this is a single base-pair
        # now we need to add the TSD. We pick a random number within the tsd range
        # tsdrange = args.tsdrange
        tsdR = range(tsdrange[0], tsdrange[1], 1) # we first make a range based on user input
        tsdL = random.sample(tsdR, 1) # then we sample 1 value in the range as our TSD length
        tsdSeq = target_chrom_seq[rnd_pos : rnd_pos + int(tsdL[0])]
        if args.verbose:
            print("TE = " + str(te) + " TSD length = " + str(tsdL))
            print("TSD = " + tsdSeq)
        # get the sequence of the TE from the bed file, add the ref bp to if in 5', this is our alternative sequence
        # we don't offset start because it is from a 0-based bed and we don't need the previous base, we insert only the TE
        alt_sequence = ref_sequence + tsdSeq + current_chrom_seq[start : end - 1]
        alt_len = len(alt_sequence)
        # get TE name
        rep_class = repeat['TE'] + "-_-" + sup
        # if args.verbose:
        #     # DEBUG
        #     print(target_chrom)
        #     print(rnd_pos)
        #     print(chrom)
        #     print(start)
        #     print(end)
        #     print(ref_sequence)
        #     print(alt_sequence)
        #     print(rep_class)
        #     print(alt_len)
        #     print(tsdSeq)
        # create the VCF line
        rec = vcfpy.Record(CHROM = target_chrom, POS = rnd_pos, ID = ['pMEI_INS_' + chrom + "_" + str(start) + "_" + str(end)],
                           REF = ref_sequence, ALT = [vcfpy.Substitution("INS", alt_sequence)],
                           QUAL = 999, FILTER = ["PASS"], INFO = {"repClass" : rep_class, "SVLEN" : alt_len, "TSD" : tsdSeq},
                           FORMAT = ["GT"],
                        calls = [
                            vcfpy.Call(sample = "sim",
                                        data = vcfpy.OrderedDict(GT = 1))]
                       )
        if args.verbose:
            print("went through!")
        # write the line
        writer.write_record(rec)
        # store the position used
        sim_pos.append(rnd_pos)
        # increase the progress bar
        bar()
    writer.close()


#######################
# simulate with simuG #
#######################
# simuG is a general purpose genome simulator written by Jia-Xing Yue (GitHub ID: yjx1217)
# Github https://github.com/yjx1217/simuG (MIT license)
simug_cmd = str("perl simuG/simuG.pl -refseq " + str(target_fasta) + " -indel_vcf " + str(out_prefix) + ".vcf -prefix " + str(out_prefix))
print('[info]' + stamp() + ' simulating genome [simuG]...')
if not args.verbose:
    simug_process = subprocess.Popen(str(simug_cmd), 
        shell = True, 
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)
else:
    simug_process = subprocess.Popen(str(simug_cmd), 
        shell = True)
simug_process.wait()

# report the new coordinates in bed format to delete in the next step
if args.verbose:
    print("[info] saving ins2del coordinates in .bed")
# wrapping a bash 1-liner in sooooo many python lines! =)
command = "paste <(awk 'NR > 1' simRef.refseq2simseq.map.txt | sort -k12,12 | cut -f 6-8) <(grep -v '#' simRef.vcf | sort -k3,3 | cut -f 8) | sed 's/repClass=//g;s/;SVLEN=/\\t/g;s/;TSD=/\\t/g' | awk '{print $1\"\\t\"$2\"\\t\"$3\"\\t\"$4\"\\t1000\\t.\\t\"$6}' > ins2del.bed"

# Execute the commands and capture the output
subprocess.run(command, shell=True, text=True, executable='/bin/bash')
#output.wait()