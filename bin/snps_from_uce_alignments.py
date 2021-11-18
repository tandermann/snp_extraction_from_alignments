#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 12:09:12 2021

@author: Tobias Andermann (tobiasandermann88@gmail.com)

requirements:
- python 3.9.7
- numpy 1.21.4
- Biopython 1.79
"""

import os, glob
from argparse import ArgumentParser
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def add_arguments(parser):
    parser.add_argument(
        '--input',
        required=True,
        metavar='dir',
        help="The directory containing all fasta-alignment files."
    )
    parser.add_argument(
        '--config',
        required=True,
        metavar='txt-file',
        help="A configuration file containing all sequence IDs that you want to extract SNPs for (one sequence ID per line)."
    )
    parser.add_argument(
        '--output',
        required=True,
        metavar='dir',
        help="The output directory where results will be saved."
    )
    parser.add_argument(
        '--include_missing',
        action='store_true',
        default=False,
        help="Use this flag if you want to include sites with missing data into the SNP alignment (leads to a higher SNP yield)."
    )
    parser.add_argument(
        '--export_nucleotides',
        action='store_true',
        default=False,
        help="Use this flag if you want to export variable positions as nucleotides (A/C/T/G) rather than binary SNPs (default, 0/1)."
    )
    parser.add_argument(
        '--extract_all',
        action='store_true',
        default=False,
        help="By default the program extracts only one SNP per alignment. Use this flag if you instead want to extract all variable sites."
    )
    parser.add_argument(
        '--phased',
        action='store_true',
        default=False,
        help="Use this flag if alignments contain phased sequences, i.e. multiple allele sequences per individual.",
    )
    parser.add_argument(
        '--delimiter',
        type=str,
        default='_',
        metavar='str',
        help="When using the '--phased' setting, specify here the delimiter that separates the different alleles from the sample name in the alignment. Example: If your alignment contains two alleles for sample1 which are named sample1_allele1 and sample1_allele2, the delimiter would be '_' (provide delimiter without quotation marks)."
    )
    parser.add_argument(
        '--seed',
        type=int,
        default=None,
        metavar='int',
        help="Set seed for reproducibility."
    )


def variable_positions(alignment,valid_chars='ACTG-',include_missing=False):
    valid_chars_list = list(valid_chars)
    alignment_T = alignment.T
    variable_chars_per_col = [list(set(i)) for i in alignment_T]
    valid_variable_chars_per_col = np.array([[j for j in i if j in valid_chars_list] for i in variable_chars_per_col],dtype=object)# reduced to only the ones we are actually counting
    all_col_ids = np.arange(alignment_T.shape[0])
    if not include_missing:
        # don't even consider columns that contain anything but information
        allowed_chars = valid_chars_list
        valid_col_ids = all_col_ids[[not any(j not in allowed_chars for j in i) for i in variable_chars_per_col]]
    else:
        valid_col_ids = all_col_ids
    # extract only those with exactly two variants
    selected_indices = valid_col_ids[[len(i)==2 for i in valid_variable_chars_per_col[valid_col_ids]]]
    # only keep positions where both variants exist in more than 1 sequence (avoiding phylogenetically uninformative sites)
    count_lists = [[list(alignment_T[i]).count(j) for j in valid_variable_chars_per_col[i]] for i in selected_indices]
    final_selected_indices = [selected_indices[id] for id,i in enumerate(count_lists) if not 1 in i]
    return(final_selected_indices)


# load arguments
parser = ArgumentParser()
add_arguments(parser)
args = parser.parse_args()
# # for trouble shooting activate the code below
# arguments = ['--input',
#              'personal_data/PHASED-DATA_all9-taxa-incomplete-mafft-nexus-edge-trimmed-fasta',
#              '--config',
#              'personal_data/snp-conf.txt',
#              '--output',
#              'personal_data/test',
#              '--snps_per_locus',
#              'one',
#              '--phased',
#              '--delimiter',
#              '_',
#              '--seed',
#              '1234']
# args = parser.parse_args(arguments)


# get input data
work_dir = args.input
out_dir = args.output
# Create the output directory
if not os.path.exists(out_dir):
    os.makedirs(out_dir)
target_taxa = np.loadtxt(args.config,dtype=str)
valid_chars = 'ACTG-'
if args.seed is None:
    seed = np.random.choice(np.arange(0,9999999))
else:
    seed = args.seed
np.random.seed(seed)

# screen output
print( "Written by Tobias Andermann")
print( "Version 2.0, November 2021")
print( "Launching SNP extraction script with seed %i ..."%seed)


if not args.phased:
    print( "\nScript is treating data as unphased alignment (add flag --phased if your data contains multiple allele sequences per individual).")
else:
    delimiter = args.delimiter
if not args.include_missing:
    print( "\nSequences with missing data are not considered during SNP extraction (add flag --include_missing for higher SNP yield, including sites with missing data).\n")
else:
    print( "\nMissing data option is activated. The final SNP alignment will contain sites with missing nucleotides, which are coded as \'?\'.\n")

fasta_files = glob.glob(os.path.join(work_dir,'*.fasta'))
if len(fasta_files) == 0:
    fasta_files = glob.glob(os.path.join(work_dir, '*.fa'))

snp_alignment = []
for fasta in fasta_files:
    print( "_" * 50)
    print( "Processing file:", fasta)
    fasta_content_tmp = SeqIO.parse(fasta, "fasta")
    fasta_content = np.array([[record.name,str(record.seq)] for record in fasta_content_tmp])
    list_sequence_names = list(fasta_content[:,0])
    length_alignment = len(fasta_content[0,1])
    # check which of the specified taxa are not present in alignment
    missing_elements_boolean = [not i in list_sequence_names for i in target_taxa]
    missing_taxa = target_taxa[missing_elements_boolean]
    if len(missing_taxa)>0:
        if args.include_missing:
            print( "%i taxa missing in alignment. These are substituted with dummy sequences:" %len(missing_taxa))
            # add dummy sequences for missing taxa
            dummy_seqs = ["N" * length_alignment for i in missing_taxa]
            fasta_content = np.vstack([fasta_content,np.array([missing_taxa,dummy_seqs]).T])
            list_sequence_names = list(fasta_content[:, 0])
        else:
            print("Alignment discarded as it does not contain all target samples. Add --include_missing flag to surpress discarding of sites with missing information.")
            continue
    # get rows that are in target_taxa in same order as provided in conf file
    row_indices = [list_sequence_names.index(i) for i in target_taxa]
    selected_fasta_content = fasta_content[row_indices]

    # find variable positions in alignment
    sequence_names = selected_fasta_content[:,0]
    alignment = np.array([list(i) for i in selected_fasta_content[:,1]])
    var_pos_list = variable_positions(alignment,valid_chars=valid_chars,include_missing=args.include_missing)

    # extract the SNP output
    if len(var_pos_list) == 0:
        print("No SNP extraction performed for '%s' due too a lack of polymorphic sites."%os.path.basename(fasta))
    else:
        if args.extract_all:
            selected_index = var_pos_list
            print("Sampling positions %s" %selected_index)
            temp_snp_align = alignment[:,selected_index]
        else:
            # chooses randomly one snp position and saves position-coordinate
            selected_index = np.random.choice(var_pos_list, 1)[0]
            print("Sampling position %s"%selected_index)
            # creates an alignment with only the extracted position
            temp_snp_align = alignment[:,selected_index].reshape(len(alignment),1)
    for i in temp_snp_align.T:
        i = i.copy()
        counts = np.array([list(i).count(j) for j in valid_chars])
        # get the most common variant as ancestral state
        ancestral_state_nucleotide = list(valid_chars)[np.argsort(counts)[-1]]
        # second common is derived state
        derived_state_nucleotide = list(valid_chars)[np.argsort(counts)[-2]]
        i[[j not in list(valid_chars) for j in i]] = '?'
        if not args.export_nucleotides:
            if args.phased:
                # find pairs of allele sequences for each sample
                reduced_names = np.array([args.delimiter.join(i.split(args.delimiter)[:-1]) for i in sequence_names])
                allele_pair_ids = [list(np.where(reduced_names==i)[0]) for i in np.unique(reduced_names)]
                # see whether sample is homozygous for ancestral, state, derived state, or heterozygous
                nucleotide_allele_pairs = np.array(['_'.join(i[j]) for j in allele_pair_ids])
                homo_ancestral = '_'.join([ancestral_state_nucleotide,ancestral_state_nucleotide])
                homo_derived = '_'.join([derived_state_nucleotide,derived_state_nucleotide])
                heterozygous = ['_'.join([ancestral_state_nucleotide,derived_state_nucleotide]),'_'.join([derived_state_nucleotide,ancestral_state_nucleotide])]
                id_ancestral = np.where(nucleotide_allele_pairs == homo_ancestral)[0]
                id_heterozygous = np.concatenate([np.where(nucleotide_allele_pairs == i)[0] for i in heterozygous])
                id_derived = np.where(nucleotide_allele_pairs == homo_derived)[0]
                nucleotide_allele_pairs[id_ancestral] = 0
                nucleotide_allele_pairs[id_heterozygous] = 1
                nucleotide_allele_pairs[id_derived] = 2
                nucleotide_allele_pairs[[j not in ['0','1','2'] for j in nucleotide_allele_pairs]] = '?'
                # export single sequence per sample
                i = nucleotide_allele_pairs
            else:
                # replace with the new values
                i[i==ancestral_state_nucleotide] = 0
                i[i==derived_state_nucleotide] = 1
        # if '?' in i:
        #     print(fasta)
        #     exit()
        snp_alignment.append(i)

# convert into numpy array
final_snp_alignment = np.array(snp_alignment).T
if args.phased and not args.export_nucleotides:
    try:
        reduced_names = np.array([args.delimiter.join(i.split(args.delimiter)[:-1]) for i in sequence_names])
    except NameError:
        exit('No suitable sites for SNP extraction were found in the alignment data. Consider using --include_missing flag to allow for sites containing missing data. Alternatively reduce the set of taxa in the config file to only samples with decent data throughout the alignments.')
    final_sequence_names = np.unique(reduced_names)
else:
    final_sequence_names = sequence_names

sequence_collection = []
for i,seq in enumerate(final_snp_alignment):
    seq_name = final_sequence_names[i]
    sequence = ''.join(seq)
    sequence_collection.append(SeqRecord(seq=Seq(sequence), id=seq_name, name=seq_name,description=''))
if args.export_nucleotides:
    outfile = os.path.join(out_dir,'snps_nucleotides.fasta')
else:
    outfile = os.path.join(out_dir,'snps_binary.fasta')


if args.export_nucleotides:
    outfile = outfile.replace('.fasta','_nucleotides.fasta')
if args.phased:
    outfile = outfile.replace('.fasta','_phased.fasta')
if args.include_missing:
    outfile = outfile.replace('.fasta','_incl_missing.fasta')
if args.extract_all:
    outfile = outfile.replace('.fasta','_all_snps.fasta')
SeqIO.write(sequence_collection, outfile, 'fasta-2line')

