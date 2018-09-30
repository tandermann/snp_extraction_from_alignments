#!/usr/local/opt/python/bin/python
#author: Tobias Hofmann, tobiashofmann@gmx.net
#workflow inspired by Yann Bertrand

#_____________________________________________________________________________________
#%%% Imports %%%
import os
import re
import sys
import glob
import shutil
import random
import argparse
import fileinput
from cogent import LoadSeqs, DNA
from cogent.core.alignment import Alignment
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped

# Complete path function
class CompletePath(argparse.Action):
	"""give the full path of an input file/folder"""
	def __call__(self, parser, namespace, values, option_string=None):
		setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Clean and trim raw Illumina read files",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--input',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing all fasta-alignment files'
	)
	parser.add_argument(
		'--config',
		required=True,
		action=CompletePath,
		help='A configuration file containing all sequence IDs that you want to extract SNPs for (one sequence ID per line)'
	)
	parser.add_argument(
		'--phased',
		action='store_true',
		default=False,
		help='Use flag if alignments contain phased sequences'
	)
	parser.add_argument(
		'--missing',
		action='store_true',
		default=False,
		help='Use flag if you want to include missing data into the SNP alignment for a higher SNP yield. WARNING: all uncertainties/ambiguities have to be coded as "N"'
	)
	parser.add_argument(
		'--base_export',
		action='store_true',
		default=False,
		help='If you want to extract variable positions as bases (rather than binary SNPs) from the alignments'
	)
	parser.add_argument(
		'--snps',
		choices=["one", "all"],
		default="one",
		help="If 'one', then only one SNP is extracted per locus (recommended for unlinked SNPs e.g. in SNAPP). If 'all', then all variable positions are extracted."
	)
	parser.add_argument(
		'--delimiter',
		type=str,
		default='_',		
		help='What is the delimiter that separates the different alleles from the sample name in the alignment? Example: If your alignment contains two alleles for sample1 which are named sample1_allele1 and sample1_allele2, the delimiter would be _'
	)
	parser.add_argument(
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The output directory where results will be saved'
	)
	return parser.parse_args()

# Get arguments
args = get_args()
# Set working directory
work_dir = args.input
out_dir = args.output
# Create the output directory
if not os.path.exists(out_dir):
	os.makedirs(out_dir)
config = args.config
name_file = open (config, 'r')
taxa_list = name_file.readlines()
taxa_names = []
for element in taxa_list:
	taxa_names.append(element.replace("\n", "").replace(" ", ""))
snp_mode = args.snps

print "\n\n"
print " ______________________________________________________"
print "|Launching SNP extraction script...                    |"
print "|                                                      |"
print "|Written by Tobias Andermann, inspired by Yann Bertrand|"
print "|Version 1.3, September 2018                           |"
print "|______________________________________________________|"

if not args.phased:
	print "\n\nScript is treating data as unphased alignment (add flag --phased to command if your data is phased)\n\n"
else:
	delimiter = args.delimiter 

if not args.missing:
	print "\nSequences with missing data are not considered during SNP extraction (add flag --missing to command for higher SNP yield, including missing data)\n\n"
else:
	print "\nThe missing data option is activated. The final SNP alignment will contain missing sites (\'?\')\n\n"
#_______________________________________________________________________________
#%%% Functions %%%
def find_names(sequence_names):
	names = []
	for element in sequence_names:
		name, allele = element.rsplit(delimiter, 1)
		names.append(name)
	return names


def variable_positions(alignment):
	var_col = []
	for x in range(len(alignment)):
		col = alignment[x].iterPositions()
		column = list(col)
		nucleotides = set(column[0])
		# if there is a 'N' in the set, just ignore it and export position if more exactly 2 states present at this position
		if not "N" in nucleotides:
			#if not args.base_export:
			count_list = []
			if len(nucleotides) == 2:
				for element in nucleotides:
					base_count = column[0].count(element)
					count_list.append(base_count)
				# this additional filter avoids positions that are only present in a single sequence (phylogenetically uninformative)
				if 1 not in count_list:
					var_col.append(x)
			# in case the base_export flag is activated, we want all variable positions, no matter how many states are present
			#else:
			#	if len(nucleotides) > 1:
			#		var_col.append(x)
	return var_col


def variable_positions_incl_missing(alignment):
	var_col = []
	# iterate through columns of alignment
	for x in range(len(alignment)):
		col = alignment[x].iterPositions()
		column = list(col)
		nucleotides = set(column[0])
		# if there is a 'N' in the set, just ignore it and export position if more exactly 2 states present at this position
		if "N" in nucleotides:
			nucleotides.remove('N')
		#if not args.base_export:
		count_list = []
		if len(nucleotides) == 2:
			for element in nucleotides:
				base_count = column[0].count(element)
				count_list.append(base_count)
			# this additional filter avoids positions that are only present in a single sequence (phylogenetically uninformative)
			if 1 not in count_list:
				var_col.append(x)
		# in case the base_export flag is activated, we want all variable positions, no matter how many states are present
		#else:
		#	if len(nucleotides) > 1:
		#		var_col.append(x)
	return var_col


def unphased_snps(list_var):

	if len(list_var)>0:
		list_positive = list_var
	else:
		print("no snp extraction performed due too a lack of polymorphic sites")
		return None

	if snp_mode == 'one':
		#chooses randomly one snp position and saves position-coordinate
		snp = random.sample(list_positive, 1)[0]
		print "sampling position", snp
		#creates an alignment with only the extracted position
		temp_snp_align = edited_alignment[snp]
	elif snp_mode == 'all':
		snp_list = list_positive
		print "sampling positions", snp_list
		temp_snp_align = edited_alignment[snp_list]
		
	#creates dictionary from the extracted snp position
	seq_dict = temp_snp_align.todict()
	new_dict = seq_dict.copy()	


	if not args.base_export:
		if snp_mode == 'one':
			set_values = []
			list_chars = list(seq_dict.values())
			for entry in set(list_chars):
				if not entry in ['N','n']:
					set_values.append((entry, list_chars.count(entry)))
			set_values.sort(key = lambda x: -x[1])
			snp_dict = {}
			# Code the base-letters into 0 and 1
			if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
				original = set_values[0][0]
				derived = set_values[1][0]
				snp_dict = {}
				for name_seq, nucleotide in seq_dict.items():
					if nucleotide == original:
						snp_dict[name_seq] = "0"
					elif nucleotide == derived:
						snp_dict[name_seq] = "1"
					else:
						print "Programming error: There are ambiguities in the final SNP dict, even though it should only contain SNPs without ambiguities. Either choose --missing flag to also include SNPs with missing data or contact programmer of this script (tobias.hofmann@bioenv.gu.se)"
		elif snp_mode == 'all':
			temp_dict = {}
			for entry in new_dict:
				temp_dict.setdefault(entry,list(new_dict[entry]))
			name = random.choice(seq_dict.keys())
			positions = range(len(list(seq_dict[name])))
			snp_dict = {}
			for pos in positions:
				# a list of all values from this position accross all samples 
				list_chars = []
				for name_seq, nucleotide in temp_dict.items():
					nucleotide_0 = nucleotide[pos]
					nucleotide_1 = nucleotide[pos]
					list_chars.append(nucleotide_0)
					list_chars.append(nucleotide_1)
				set_values = []
				for entry in set(list_chars):
					if not entry in ['N','n']:
						set_values.append((entry, list_chars.count(entry)))
				set_values.sort(key = lambda x: -x[1])

				if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
					original = set_values[0][0]
					derived = set_values[1][0]
					for name_seq, nucleotide in temp_dict.items():
						nucleotide_pos = nucleotide[pos]
						if nucleotide_pos == original:
							snp_value = "0"
						elif nucleotide_pos == derived:
							snp_value = "1"
						else:
							print "Programming error: There are ambiguities in the final SNP dict, even though it should only contain SNPs without ambiguities. Either choose --missing flag to also include SNPs with missing data or contact programmer of this script (tobias.hofmann@bioenv.gu.se)"
						snp_dict.setdefault(name_seq,[])
						snp_dict[name_seq].append(snp_value)
	else:
		snp_dict = seq_dict.copy()
	return snp_dict


def phased_snps(list_var):
	if len(list_var)>0:
		list_positive = list_var
	else:
		print("no SNP extraction performed due too a lack of polymorphic sites")
		return None

	if snp_mode == 'one':
		#chooses randomly one snp position and saves position-coordinate
		snp = random.sample(list_positive, 1)[0]
		print "sampling position", snp
		#creates an alignment with only the extracted position
		temp_snp_align = edited_alignment[snp]
	elif snp_mode == 'all':
		snp_list = list_positive
		print "sampling positions", snp_list
		temp_snp_align = edited_alignment[snp_list]		

	#creates dictionary from the extracted snp position
	seq_dict = temp_snp_align.todict()

	if not args.base_export:	
		#finds elements in the dictionary that belong to the same sample (in case of multiple alleles)
		clean_names = sorted(find_names(taxa_names))
		no_duplicates = list(set(clean_names))
		new_dict = {taxon : [value for key, value in seq_dict.items() if key.startswith(taxon)] for taxon in no_duplicates}
	else:
		new_dict = seq_dict.copy()		



	if not args.base_export:
		if snp_mode == 'one':
			#returns which two bases are present in the dictionary (ordered, to replace them properly)
			set_values = []
			list_chars = list(seq_dict.values())
			for entry in set(list_chars):
				if not entry in ['N','n']:
					set_values.append((entry, list_chars.count(entry)))
			set_values.sort(key = lambda x: -x[1])
			snp_dict = {}
			# Code the base-letters into 0, 1 or 2
			if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
				original = set_values[0][0]
				derived = set_values[1][0]
				for name_seq, nucleotide in new_dict.items():
					if nucleotide[0] == nucleotide[1] == original:
						snp_dict[name_seq] = "0"
					elif nucleotide[0] == nucleotide[1] == derived:
						snp_dict[name_seq] = "2"
					elif nucleotide[0] != nucleotide[1] and nucleotide[0] != "N" and nucleotide[1] != "N":
						snp_dict[name_seq] = "1"
					else:
						print "Programming error: There are ambiguities in the final SNP dict, even though it should only contain SNPs without ambiguities. Either choose --missing flag to also include SNPs with missing data or contact programmer of this script (tobias.hofmann@bioenv.gu.se)"


		elif snp_mode == 'all':
			temp_dict = {}
			for entry in new_dict:
				temp_dict.setdefault(entry,list(new_dict[entry]))
			name = random.choice(seq_dict.keys())
			positions = range(len(list(seq_dict[name])))
			snp_dict = {}
			for pos in positions:
				# a list of all values from this position accross all samples 
				list_chars = []
				for name_seq, nucleotide in temp_dict.items():
						list_bases = []
						for element in nucleotide:
							list_bases.append(list(element))
						nucleotide_0 = list_bases[0][pos]
						nucleotide_1 = list_bases[1][pos]
						list_chars.append(nucleotide_0)
						list_chars.append(nucleotide_1)
				set_values = []
				for entry in set(list_chars):
					if not entry in ['N','n']:
						set_values.append((entry, list_chars.count(entry)))
				set_values.sort(key = lambda x: -x[1])
				
				# Code the base-letters into 0, 1 or 2
				if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
					original = set_values[0][0]
					derived = set_values[1][0]
					for name_seq, nucleotide in temp_dict.items():
						# make a list with all SNPs of the sample for the locus
						list_bases = []
						for element in nucleotide:
							list_bases.append(list(element))
						nucleotide_0 = list_bases[0][pos]
						nucleotide_1 = list_bases[1][pos]
						snp_value = ""
						if nucleotide_0 == nucleotide_1 == original:
							snp_value = "0"
						elif nucleotide_0 == nucleotide_1 == derived:
							snp_value = "2"
						elif nucleotide_0 != nucleotide_1 and nucleotide_0 != "N" and nucleotide_1 != "N":
							snp_value = "1"
						else:
							print "Programming error: There are ambiguities in the final SNP dict, even though it should only contain SNPs without ambiguities. Either choose --missing flag to also include SNPs with missing data or contact programmer of this script (tobias.hofmann@bioenv.gu.se)"
						
						snp_dict.setdefault(name_seq,[])
						snp_dict[name_seq].append(snp_value)
	else:
		snp_dict = new_dict.copy()
	return snp_dict


def unphased_snps_missing(list_var):
	if len(list_var)>0:
		list_positive = list_var
	elif len(variable_positions_incl_missing(edited_alignment))>0:
		list_positive = variable_positions_incl_missing(edited_alignment)
		print "Variable positions including ambiguities (N):", list_positive
	else:
		print("no snp extraction performed due too a lack of polymorphic sites")
		return None

	if snp_mode == 'one':
		#chooses randomly one snp position and saves position-coordinate
		snp = random.sample(list_positive, 1)[0]
		print "sampling position", snp
		#creates an alignment with only the extracted position
		temp_snp_align = edited_alignment[snp]
	elif snp_mode == 'all':
		snp_list = list_positive
		print "sampling positions", snp_list
		temp_snp_align = edited_alignment[snp_list]
		
	#creates dictionary from the extracted snp position
	seq_dict = temp_snp_align.todict()
	new_dict = seq_dict.copy()	


	if not args.base_export:
		if snp_mode == 'one':
			set_values = []
			list_chars = list(seq_dict.values())
			for entry in set(list_chars):
				if not entry in ['N','n']:
					set_values.append((entry, list_chars.count(entry)))
			set_values.sort(key = lambda x: -x[1])
			snp_dict = {}
			# Code the base-letters into 0 and 1
			if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
				original = set_values[0][0]
				derived = set_values[1][0]
				snp_dict = {}
				for name_seq, nucleotide in seq_dict.items():
					if nucleotide == original:
						snp_dict[name_seq] = "0"
					elif nucleotide == derived:
						snp_dict[name_seq] = "1"
					else:
						snp_dict[name_seq] = "?"
		elif snp_mode == 'all':
			temp_dict = {}
			for entry in new_dict:
				temp_dict.setdefault(entry,list(new_dict[entry]))
			name = random.choice(seq_dict.keys())
			positions = range(len(list(seq_dict[name])))
			snp_dict = {}
			for pos in positions:
				# a list of all values from this position accross all samples 
				list_chars = []
				for name_seq, nucleotide in temp_dict.items():
					nucleotide_0 = nucleotide[pos]
					nucleotide_1 = nucleotide[pos]
					list_chars.append(nucleotide_0)
					list_chars.append(nucleotide_1)
				set_values = []
				for entry in set(list_chars):
					if not entry in ['N','n']:
						set_values.append((entry, list_chars.count(entry)))
				set_values.sort(key = lambda x: -x[1])

				if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
					original = set_values[0][0]
					derived = set_values[1][0]
					for name_seq, nucleotide in temp_dict.items():
						nucleotide_pos = nucleotide[pos]
						snp_value = ""
						if nucleotide_pos == original:
							snp_value = "0"
						elif nucleotide_pos == derived:
							snp_value = "1"
						else:
							snp_value = "?"
						snp_dict.setdefault(name_seq,[])
						snp_dict[name_seq].append(snp_value)

	else:
		snp_dict = seq_dict.copy()
	return snp_dict


def phased_snps_missing(list_var):
	if len(list_var)>0:
		list_positive = list_var
	elif len(variable_positions_incl_missing(edited_alignment))>0:
		list_positive = variable_positions_incl_missing(edited_alignment)
		print "Variable positions including ambiguities (N):", list_positive
	else:
		print("no SNP extraction performed due too a lack of polymorphic sites")
		return None

	if snp_mode == 'one':
		#chooses randomly one snp position and saves position-coordinate
		snp = random.sample(list_positive, 1)[0]
		print "sampling position", snp
		#creates an alignment with only the extracted position
		temp_snp_align = edited_alignment[snp]
	elif snp_mode == 'all':
		snp_list = list_positive
		print "sampling positions", snp_list
		temp_snp_align = edited_alignment[snp_list]		

	#creates dictionary from the extracted snp position
	seq_dict = temp_snp_align.todict()

	if not args.base_export:	
		#finds elements in the dictionary that belong to the same sample (in case of multiple alleles)
		clean_names = sorted(find_names(taxa_names))
		no_duplicates = list(set(clean_names))
		new_dict = {taxon : [value for key, value in seq_dict.items() if key.startswith(taxon)] for taxon in no_duplicates}
	else:
		new_dict = seq_dict.copy()		



	if not args.base_export:
		if snp_mode == 'one':
			#returns which two bases are present in the dictionary (ordered, to replace them properly)
			set_values = []
			list_chars = list(seq_dict.values())
			for entry in set(list_chars):
				if not entry in ['N','n']:
					set_values.append((entry, list_chars.count(entry)))
			set_values.sort(key = lambda x: -x[1])
			snp_dict = {}
			# Code the base-letters into 0, 1 or 2
			if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
				original = set_values[0][0]
				derived = set_values[1][0]
				for name_seq, nucleotide in new_dict.items():
					if nucleotide[0] == nucleotide[1] == original:
						snp_dict[name_seq] = "0"
					elif nucleotide[0] == nucleotide[1] == derived:
						snp_dict[name_seq] = "2"
					elif nucleotide[0] != nucleotide[1] and nucleotide[0] != "N" and nucleotide[1] != "N":
						snp_dict[name_seq] = "1"
					else:
						snp_dict[name_seq] = "?"

		elif snp_mode == 'all':
			temp_dict = {}
			for entry in new_dict:
				temp_dict.setdefault(entry,list(new_dict[entry]))
			name = random.choice(seq_dict.keys())
			positions = range(len(list(seq_dict[name])))
			snp_dict = {}
			for pos in positions:
				# a list of all values from this position accross all samples 
				list_chars = []
				for name_seq, nucleotide in temp_dict.items():
						list_bases = []
						for element in nucleotide:
							list_bases.append(list(element))
						nucleotide_0 = list_bases[0][pos]
						nucleotide_1 = list_bases[1][pos]
						list_chars.append(nucleotide_0)
						list_chars.append(nucleotide_1)
				set_values = []
				for entry in set(list_chars):
					if not entry in ['N','n']:
						set_values.append((entry, list_chars.count(entry)))
				set_values.sort(key = lambda x: -x[1])
				
				# Code the base-letters into 0, 1 or 2
				if "A" or "C" or "T" or "G" or "a" or "c" or "t" or "g" in set_values:
					original = set_values[0][0]
					derived = set_values[1][0]
					for name_seq, nucleotide in temp_dict.items():
						# make a list with all SNPs of the sample for the locus
						list_bases = []
						for element in nucleotide:
							list_bases.append(list(element))
						nucleotide_0 = list_bases[0][pos]
						nucleotide_1 = list_bases[1][pos]
						snp_value = ""
						if nucleotide_0 == nucleotide_1 == original:
							snp_value = "0"
						elif nucleotide_0 == nucleotide_1 == derived:
							snp_value = "2"
						elif nucleotide_0 != nucleotide_1 and nucleotide_0 != "N" and nucleotide_1 != "N":
							snp_value = "1"
						else:
							snp_value = "?"
						
						snp_dict.setdefault(name_seq,[])
						snp_dict[name_seq].append(snp_value)
	else:
		snp_dict = new_dict.copy()
	return snp_dict


#_____________________________________________________________________________________
#%%% Workflow %%%
#get all fasta files in working directory
os.chdir(work_dir)
fasta_files = [x for x in os.listdir(os.getcwd()) if ".fasta" in x and "snp.fasta" not in x]
final_dict = {}
list_of_genes=[]
for fasta in fasta_files[:]:
	print "_" * 50
	print "processing file:", fasta
	aln = LoadSeqs(fasta, moltype=DNA)
	list_sequences = aln.Names
	# Check if all taxa that are specified in the control file exist in the alignment
	list_check = set(taxa_names).issubset(set(list_sequences))
	# If there are some taxa not present in alignment the following code will simulate a sequence full of "N" of the correct length for those missing taxa and add it to the alignment
	if list_check == False:
		missing_elements = []
		for element in taxa_names:
			#print element,list_sequences,element in list_sequences
			if element not in list_sequences:
				missing_elements.append(element)
		print "These taxa are missing in alignment:", missing_elements, "\nSequences for missing taxa will be generated only containing \"N\"."
		seq = Alignment(aln)
		string_list = seq.todict().values()
		length_alignment = ""
		for element in string_list:
			length_alignment = len(element)
		simulated_seq = []
		for element in missing_elements:
			fake_aln = "N" * length_alignment
			simulated_seq.append((element,fake_aln))
		fake_seqs = LoadSeqs(data = simulated_seq)
		aln = aln.addSeqs(fake_seqs)
	# Apply filter of user-set taxa names to be used for snp-extraction
	edited_alignment = aln.takeSeqs(sorted(taxa_names))
	# Get the variable positions for each fasta file
	if not args.missing:
		var_pos_list = variable_positions(edited_alignment)
	else:
		var_pos_list = variable_positions_incl_missing(edited_alignment)
	print var_pos_list
	D = ""
	if not args.phased:
		if not args.missing:
			D = unphased_snps(var_pos_list)
		else:
			D = unphased_snps_missing(var_pos_list)
	else:
		if not args.missing:
			D = phased_snps(var_pos_list)
		else:
			D = phased_snps_missing(var_pos_list)
	if D != None:
		list_of_genes.append(fasta)
		for key, values in D.items():
			if key not in final_dict.keys():
				final_dict[key] = [D[key]]
			else:
				final_dict[key].append(D[key])

#join all snps into one dictionary
final_snp_alignment = {}
if snp_mode == 'one':
	for key, value in final_dict.items():
		final_snp_alignment[key] = "".join(value)
elif snp_mode == 'all':
	for key, values in final_dict.items():
		value = sum(values, [])
		final_snp_alignment[key] = "".join(value)

# Create the output file in output directory
output_file_fasta = os.path.join(out_dir,'snp.fasta')
#print the snp dictionary into a fasta-file
with open(output_file_fasta, "wb") as f:
	for k, v in final_snp_alignment.items():
		f.write(">" + k+ "\n")
		f.write(v+ "\n")

# Create output file for SNAPP
output_file_nexus = os.path.join(out_dir,'snp.nexus')
aln = AlignIO.read(open(output_file_fasta), "fasta", alphabet=Gapped(IUPAC.ambiguous_dna))
with open(output_file_nexus, "wb") as n:
	n.write(aln.format("nexus"))
if not args.phased:
	for line in fileinput.input(output_file_nexus, inplace = 1):
		print line.replace("format datatype=dna missing=? gap=-;", "format datatype=binary symbols=01 missing=?;").rstrip()
else:
	for line in fileinput.input(output_file_nexus, inplace = 1):
		print line.replace("format datatype=dna missing=? gap=-;", "format datatype=integerdata symbols=\"012\" missing=?;").rstrip()
