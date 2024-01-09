#!/usr/bin/python3

import vcfpy
import argparse
import pandas
import pysam
import random
import sys
import statistics
import subprocess

# Define a custom argument type for a list of integers
def list_of_ints(arg):
    return list(map(int, arg.split(',')))

# create the parser for the user input
parser = argparse.ArgumentParser(description='Create a simulated VCF file with random TE insertions or deletions')
parser.add_argument('-v', '--variants', type = int, metavar='N', help='the total number of variants to simulate', dest = 'nb_var', default = 100)
parser.add_argument('-r', '--ratio', type = float, metavar='[0-1]', help='insertion/deletion ratio', dest = 'ratio', default = .7)
parser.add_argument('-g', '--genome', type = str, metavar='STR', help='reference genomes (fa/fa.gz)', dest = 'ref_genome', default = '../simData/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.gz') # same ref genome to get the INS TE from
parser.add_argument('-t', '--target_chr', type = str, metavar='STR', help='name of target chromosome where insertion/deletion will be performed', dest = 'target_chrom', default = 'chr22')
parser.add_argument('-f', '--target_fasta', type = str, metavar='STR', help='fasta file for the target chromosome', dest = 'target_fasta', default = 'simRef.simseq.genome.fa') # now the target is the refSim genome (chr22)
parser.add_argument('-b', '--bed', type=str, metavar='STR', help='bed file with reference insertions to use', dest = 'bed_in', default = '../simData/hg38.AluY.L1HSPA2.SVA_EF.bed') # same bed to take the INS coordinates
parser.add_argument('-B', '--bed2', type=str, metavar='STR', help='bed file with simRef TEs to delete', dest = 'bed2_in', default = 'ins2del.bed') # new bed to take the DEL coordinates
parser.add_argument('-R', '--tsdrange', type = list_of_ints, metavar='min,max', help='TSD length range', dest = 'tsdrange', default = '4,22' )
parser.add_argument('-o', '--out', type=str, metavar='STR', help='output file prefix (will be <prefix>.vcf)', dest = 'out_prefix', default = 'simAlt')
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
    print("warning: you have asked for an odd number of variants  " + str(nb_var-1) + ",  " + str(nb_var) + " will be used instead.")
ref_genome = args.ref_genome
target_fasta = args.target_fasta
target_chrom = args.target_chrom
out_prefix = args.out_prefix
# force the Alu, L1, SVA ratio to realistic values
te_props = [0.7,0.2,0.1]
# calculate ins/del number based on input
ins_nb = round(ratio*nb_var)
# here we don't need to sample the deletions are they are all set in the ins2del.bed file produced by step 01
# del_nb = nb_var - ins_nb #round(nb_var-(ratio*nb_var))
# load the reference genome
fasta = pysam.FastaFile(ref_genome)
# read and store the input bed file
bed_in_pre = pandas.read_csv(args.bed_in, sep = '\t',
                          names = ['chrom', 'start', 'end', 'TE', 'score', 'strand'])
# make sure the chromosomes of the bed file are in the reference genome
bed_in = bed_in_pre[bed_in_pre['chrom'].isin(fasta.references)]
# check input length to make sure we don't sample more than there is; we don't need to check the deletion as their number is already set
in_ins_nb = len(bed_in[-bed_in['chrom'].isin([target_chrom])])
if in_ins_nb < ins_nb:
    print('*******')
    print('ERROR!: not enough TEs in the bed file to insert in the target chromosome!')
    print('*******')
    parser.print_help(sys.stderr)
    sys.exit(1)

# read and store the loci to delete from ins2del.bed
bed_del = pandas.read_csv(args.bed2_in, sep = '\t',
                          names = ['chrom', 'start', 'end', 'TE', 'score', 'strand', 'TSD'])

# sample the insertions (future deletions) according to predefined proportions
i_alu = bed_in[-bed_in['chrom'].isin([target_chrom]) & bed_in.TE.str.match('Alu') & (bed_in['end']-bed_in['start'] >= 250)].sample(round(max(ins_nb*te_props[0],1)))
i_line = bed_in[-bed_in['chrom'].isin([target_chrom]) & bed_in.TE.str.match('L1') & (bed_in['end']-bed_in['start'] >= 900)].sample(round(max(ins_nb*te_props[1],1)))
i_sva = bed_in[-bed_in['chrom'].isin([target_chrom]) & bed_in.TE.str.match('SVA') & (bed_in['end']-bed_in['start'] >= 900)].sample(round(max(ins_nb*te_props[2],1)))
# here I replaced "bed_ins" with "repmask_subset", as it is the file used later and I don't want to change names
repmask_subset = pandas.concat([i_alu, i_line, i_sva], axis=0).reset_index().sort_values(by=['chrom', 'start'])
# concatenate
# repmask_subset = pandas.concat([bed_del, bed_ins], axis=0).reset_index().sort_values(by=['chrom', 'start'])
if args.verbose:
    print(bed_ins)
# read and store the header
reader = vcfpy.Reader.from_path('header.vcf')
# set the name of the sample in the genotype column for the header we just stored
reader.header.samples = vcfpy.SamplesInfos(['sim'])
# initialize the output file. For now it only includes the header with the modified sample name
writer = vcfpy.Writer.from_path(out_prefix + '.vcf', reader.header)


#########################
# Simulate TE ins / del #
#########################
# subset the bed file (repmask) for the chromosomes present in the reference genome
# repmask_subset = repmask[repmask['chrom'].isin(fasta.references) ]
# initialize the current chromosome and its sequence
current_chrom = None
current_chrom_seq = None
# define which chromosome will receive the annotations
target_chrom_seq = fasta[target_chrom]
# create an empty list to store the positions of the simulated ins/del (1-based)
sim_pos = []
# loop over each line of the bed file
print('creating pMEI insertions/deletions VCF entries...')
if args.verbose:
    print(repmask_subset)
for index, repeat in repmask_subset.iterrows():
    chrom = repeat['chrom']
    start = repeat['start']
    end = repeat['end']

    # update in memory chromosome sequence
    if current_chrom != chrom:
        #print("Switching to " + chrom)
        current_chrom = chrom
        current_chrom_seq = fasta[current_chrom]
    
    # we skip the next part because we don't want to make deletion at this point

    # # if the bed chromosome is the same as the target chromosome, we will do a deletion
    # if current_chrom == target_chrom:    
    #     # for each line, update the reference and alternative sequence
    #     if args.verbose:
    #         print("creating deletions in: " + target_chrom + '...')
    #     # start is 0-based (bed), so no need to offset, but we need to pick 1 base before start to get the ref with the alt allele in 5'
    #     # end is 1-based, so we -1 it to get the right position in the string
    #     ref_sequence = current_chrom_seq[start - 1 : end - 1]
    #     alt_sequence = ref_sequence[0]
    #     rep_class = repeat['TE']

    #     # create the VCF line
    #     rec = vcfpy.Record(CHROM = chrom, POS = start, ID = ['pMEI_DEL_' + chrom + "_" + str(start) + "_" + str(end)],
    #                        REF = ref_sequence, ALT = [vcfpy.Substitution("DEL", alt_sequence)],
    #                        QUAL = 999, FILTER = ["PASS"], INFO = {"repClass" : rep_class, "SVLEN" : - (end - start)},
    #                        FORMAT = ["GT"],
    #                     calls = [
    #                         vcfpy.Call(sample = "sim",
    #                                     data = vcfpy.OrderedDict(GT = 1))]
    #                    )
    #     # write the line
    #     writer.write_record(rec)
    #     # store the position used
    #     sim_pos.append(start)
    # # else, we will do a random insertion
    # else:
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
    tsdR = range(tsdrange[0], tsdrange[1], 1) # we first make a range based on user input
    tsdL = random.sample(tsdR, 1) # then we sample 1 value in the range as our TSD length
    tsdSeq = current_chrom_seq[end : end + int(tsdL)]
    if args.verbose:
        print("TSD = " + tsdSeq)
    # get the sequence of the TE from the bed file, add the ref bp to if in 5', this is our alternative sequence
    # we don't offset start because it is from a 0-based bed and we don't need the previous base, we insert only the TE
    alt_sequence = ref_sequence + tsdSeq + current_chrom_seq[start : end - 1]
    alt_len = len(alt_sequence)
    # get TE name
    rep_class = repeat['TE']
    # create the VCF line
    rec = vcfpy.Record(CHROM = target_chrom, POS = rnd_pos, ID = ['pMEI_INS_' + chrom + "_" + str(start) + "_" + str(end)],
                       REF = ref_sequence, ALT = [vcfpy.Substitution("INS", alt_sequence)],
                       QUAL = 999, FILTER = ["PASS"], INFO = {"repClass" : rep_class, "SVLEN" : alt_len, "TSD" : tsdSeq},
                       FORMAT = ["GT"],
                    calls = [
                        vcfpy.Call(sample = "sim",
                                    data = vcfpy.OrderedDict(GT = 1))]
                   )
    # write the line
    writer.write_record(rec)
    # store the position used
    sim_pos.append(rnd_pos)

#############################
# simulate random intervals #
#############################

print('creating background insertions/deletions VCF entries...')
# create a list of random position on the target chrom
# they must not be already occupied by a simulated pMEI
randit_pos = []
for i in range(nb_var):
    j = random.randint(1,fasta.get_reference_length(target_chrom))
    ref_sequence = target_chrom_seq[j - 1]
    while ref_sequence == 'N' or j in sim_pos:
        j = random.randint(1,fasta.get_reference_length(target_chrom))
        ref_sequence = target_chrom_seq[j - 1]
    randit_pos.append(j)

# load the distribution of SV length for HG002
sv_len_in = pandas.read_csv('HG002_SVs_Tier1_v0.6_filtered_SV_length', sep = '\t', names = ['SVchr', 'SVlen'])
# filter size to 100bp 
sv_len = sv_len_in[(sv_len_in['SVlen'] < -100) | (sv_len_in['SVlen'] > 100)]
# sample them
randit_len = sv_len['SVlen'].sample(len(randit_pos))
# sample some chromosomes
randit_ins_chr = sv_len[sv_len['SVchr'] != 'chr22'].sample(round(.5*nb_var))['SVchr']
randit_del_chr = sv_len[sv_len['SVchr'] == 'chr22'].sample(round(.5*nb_var))['SVchr']
#randit_chr = randit_ins_chr.append(randit_del_chr)
randit_chr = pandas.concat([randit_ins_chr, randit_del_chr], axis = 0, ignore_index = True) #sv_len['SVchr'].sample(len(randit_pos))
# combine in a table to iterate over
randit_table = pandas.DataFrame(zip(randit_chr, randit_len), columns = ['SVchr','SVlen'])
# create an empty table to store the validated random intervals
val_rnd = pandas.DataFrame(columns=['chrom', 'start', 'end', 'name', 'type', 'ins_site'])
# gather the list of simulated deletions
dels = repmask_subset[repmask_subset['chrom'].isin(list(target_chrom))]
# iterate the random intervals table to validate them (i.e.: check they don't intersect a simulated deletion)
for index, rnd in randit_table.iterrows():
    chrom = rnd['SVchr']
    rndlen = abs(rnd['SVlen'])
    # we take a random base on the chromosome (being careful not to go overboard)
    if chrom == target_chrom:
        start = randit_pos[index]
        # check if the interval we will create is not intersecting a TE deletion on the target chromosome!
        end = start + rndlen
        if args.verbose:
            print('random del...' + str(start) + ' ' + str(end))
        test = 'false'
        while test == 'false':
            if args.verbose:
                print('enters while loop...')
            for i in range(len(dels)):
                if max(start, dels['start'].values[i-1]) < min(end, dels['end'].values[i-1]):
                    # we retry
                    if args.verbose:
                        print('intersect simulated del ' + str(dels['start'].values[i-1]) + ' ' + str(dels['end'].values[i-1]) + ' Regenerate breakpoint...')
                    start = random.randint(rndlen + 1,fasta.get_reference_length(chrom) - rndlen)
                    test = 'false'
                    break
                else:
                    if args.verbose:
                        print('pass simulated del ' + str(dels['start'].values[i-1]) + ' ' + str(dels['end'].values[i-1]) + ' continuing...')
            test = 'pass'
        if args.verbose:
            print('random del...' + str(start) + ' ' + str(end) + ' passed!')
        # increment the list
        ins_site = start
        row = {'chrom':chrom, 'start':start, 'end':end, 'name': 'sim' + str(index), 'type':'DEL', 'ins_site':ins_site}
        newline = pandas.DataFrame([row])
        val_rnd = pandas.concat([val_rnd, newline], axis = 0, ignore_index = True)
    else:
        start = random.randint(rndlen + 1,fasta.get_reference_length(chrom) - rndlen)
        # check if the interval we will create is not intersecting a TE deletion on the target chromosome!
        end = start + rndlen
        # check the sequence we insert is not made of 'N'
        current_chrom_seq = fasta[chrom]
        ins_sequence = current_chrom_seq[start : end - 1]
        while ins_sequence.count('N')/len(ins_sequence) > 0.10:
            start = random.randint(rndlen + 1,fasta.get_reference_length(chrom) - rndlen)
            end = start + rndlen
            ins_sequence = current_chrom_seq[start : end - 1]
        # pick the insertion site on the target chromosome and increment the list
        ins_site = randit_pos[index]
        row = {'chrom':chrom, 'start':start, 'end':end, 'name': 'sim' + str(index), 'type':'INS', 'ins_site':ins_site}
        newline = pandas.DataFrame([row])
        val_rnd = pandas.concat([val_rnd, newline], axis = 0, ignore_index = True)

if args.verbose:
    print('all bgSV passed checks!')
    print('the background SV will be:')
    print(val_rnd)

##### Now write the bgSV in the VCF:
for index, repeat in val_rnd.iterrows():
    chrom = repeat['chrom']
    start = repeat['start']
    end = repeat['end']
    ins_site = repeat['ins_site']
    # update in memory chromosome sequence
    if current_chrom != chrom:
        #print("Switching to " + chrom)
        current_chrom = chrom
        current_chrom_seq = fasta[current_chrom]
    # if the bed chromosome is the same as the target chromosome, we will do a deletion
    if current_chrom == target_chrom:    
        # for each line, update the reference and alternative sequence
        if args.verbose:
            print("creating background deletions in: " + target_chrom + '...')
        # start is 0-based (bed), so no need to offset, but we need to pick 1 base before start to get the ref with the alt allele in 5'
        # end is 1-based, so we -1 it to get the right position in the string
        ref_sequence = current_chrom_seq[start - 1 : end - 1]
        alt_sequence = ref_sequence[0]
        rep_class = repeat['type']

        # create the VCF line
        rec = vcfpy.Record(CHROM = chrom, POS = start, ID = ['bgSV_DEL_' + chrom + "_" + str(start) + "_" + str(end)],
                           REF = ref_sequence, ALT = [vcfpy.Substitution("DEL", alt_sequence)],
                           QUAL = 999, FILTER = ["PASS"], INFO = {"repClass" : rep_class, "SVLEN" : - (end - start)},
                           FORMAT = ["GT"],
                        calls = [
                            vcfpy.Call(sample = "sim",
                                        data = vcfpy.OrderedDict(GT = 1))]
                       )
        # write the line
        writer.write_record(rec)
        # store the position used
        sim_pos.append(start)
    # else, we will do a random insertion
    else:
        if args.verbose:
            print("creating background insertions from: " + chrom + '...')
        rnd_pos = start
        # get the 1 base at the site, will be the reference sequence
        # we assume rnd_pos is 1-based. Thus we -1 it to get it in python
        ref_sequence = target_chrom_seq[ins_site - 1]
        # get the sequence of the TE from the bed file, add the ref bp to if in 5', this is our alternative sequence
        # we don't offset start because it is from a 0-based bed and we don't need the previous base, we insert only the TE
        alt_sequence = ref_sequence + current_chrom_seq[start : end - 1]
        # get SV name
        rep_class = repeat['type']
        # create the VCF line
        rec = vcfpy.Record(CHROM = target_chrom, POS = ins_site, ID = ['bgSV_INS_' + chrom + "_" + str(start) + "_" + str(end)],
                           REF = ref_sequence, ALT = [vcfpy.Substitution("INS", alt_sequence)],
                           QUAL = 999, FILTER = ["PASS"], INFO = {"repClass" : rep_class, "SVLEN" : end - start},
                           FORMAT = ["GT"],
                        calls = [
                            vcfpy.Call(sample = "sim",
                                        data = vcfpy.OrderedDict(GT = 1))]
                       )
        # write the line
        writer.write_record(rec)
        # store the position used
        sim_pos.append(rnd_pos)

# close the VCF file (write it)
writer.close()

#######################
# simulate with simuG #
#######################
# simuG is a general purpose genome simulator written by Jia-Xing Yue (GitHub ID: yjx1217)
# Github https://github.com/yjx1217/simuG (MIT license)
simug_cmd = str("perl simuG/simuG.pl -refseq " + str(target_fasta) + " -indel_vcf " + str(out_prefix) + ".vcf -prefix " + str(out_prefix))
print('simulating genome [simuG]...')
if not args.verbose:
    simug_process = subprocess.Popen(str(simug_cmd), 
        shell = True, 
        stdout=subprocess.DEVNULL,
        stderr=subprocess.STDOUT)
else:
    simug_process = subprocess.Popen(str(simug_cmd), 
        shell = True)
simug_process.wait()