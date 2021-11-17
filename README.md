# snp_extraction_from_alignments
*Written by Tobias Andermann (tobiasandermann88@gmail.com)*

This repository contains a freely accessible program for the extraction of Single Nucleotide Polymorphisms (SNPs) from multiple sequence alignments. To run the program, you can [download this GitHub repo](https://github.com/tobiashofmann88/snp_extraction_from_alignments/archive/refs/heads/master.zip). This includes the program as well as some example data and this readme file.

The program is designed to extract binary positions, i.e. sites with exactly two nucleotides variants, which each are shared by at least two samples. The latter ensures that the extracted sites are phylogenetically informative. The output SNP alignment is in fasta format and can be used for downstream applications, such as [SNAPP](https://doi.org/10.1093/molbev/mss086) or [Admixture](https://doi.org/10.1101/gr.094052.109).


The program is written in Python and can be executed from your **bash command line**. I recommend running it with [**Python3.9**](https://www.python.org/downloads/). Older Python3 versions may work as well but are not tested. Other requirements are the biopython and numpy libraries, which can be installed with pip:

```bash
pip install numpy
pip install Biopython
```

To run the program you have to call python, followed by the path pointing to the python script:

```bash
python bin/snps_from_uce_alignments.py
```

Use the help function (`-h`) for more detailed explanation of individual settings. This should give you this output:

```
python bin/snps_from_uce_alignments.py -h

usage: snps_from_uce_alignments.py [-h] --input dir --config txt-file --output dir [--include_missing] [--export_nucleotides] [--snps_per_locus str] [--phased]
                                   [--delimiter str] [--seed int]

optional arguments:
  -h, --help            show this help message and exit
  --input dir           The directory containing all fasta-alignment files.
  --config txt-file     A configuration file containing all sequence IDs that you want to extract SNPs for (one sequence ID per line).
  --output dir          The output directory where results will be saved.
  --include_missing     Use this flag if you want to include sites with missing data into the SNP alignment (leads to a higher SNP yield).
  --export_nucleotides  Use this flag if you want to export variable positions as nucleotides (A/C/T/G) rather than binary SNPs (default, 0/1).
  --snps_per_locus str  If 'one' (default), then only one SNP is extracted per locus (recommended for unlinked SNPs, e.g. for use in SNAPP). If 'all', then all variable
                        positions are extracted.
  --phased              Use this flag if alignments contain phased sequences, i.e. multiple allele sequences per individual.
  --delimiter str       When using the '--phased' setting, specify here the delimiter that separates the different alleles from the sample name in the alignment. Example:
                        If your alignment contains two alleles for sample1 which are named sample1_allele1 and sample1_allele2, the delimiter would be '_' (provide
                        delimiter without quotation marks).
  --seed int            Set seed for reproducibility.

```

The script requires a folder with fasta-formatted alignment files of different loci. The naming of samples across alignments should be consistent. The program is designed to work for incomplete alignments, i.e. it is no problem if some taxa are missing in some of the alignments. In this case it is recommended to use the `--include_missing ` settings (see more information below).

Besides the input alignments, the program also requires a `.txt` file with a list of the sample names to extract SNPs for. This allows more flexibility to extract SNPs for different subsets of the samples. Examples of input files can be found in the `example_files/` folder.

Below you can find some different use cases demonstrated on example data available from the [GitHub repo](https://github.com/tobiashofmann88/snp_extraction_from_alignments).

## Extracting SNPs from alignments

### Simulated, clean alignments
These alignments were simulated with the simulation program BPP. They are therefore free of missing sites or ambiguities. This is usually not the case for empirical data, as we see in the second example below.

The alignments contain two allele sequences for each sample per locus. These sequences are named in the following manner in the fasta file, indicating the sequences for allele 0 and allele 1 for sampleA:

```
>sampleA_0
>sampleA_1
```

We can run the SNP extraction on these data as shown below. The command includes the `--phased` flag, which tells the program to look for multiple allele sequences for each sample, sharing the same name stem. To tell the program how the different allele sequence headers are delimited, you can use the `--delimiter` flag followed by the delimiter that sets apart the sequence headers belonging to the same sample. In this example they are delimited from the sample name stem by `_`, ending in `_0` and `_1`. We're also telling the program to only extract a single SNP per locus (`--snps_per_locus one`) to ensure unlinked SNPs. Alternatively this setting can be set to `--snps_per_locus all`.

Now lets run the command:

`python bin/snps_from_uce_alignments.py --input example_files/alignments/simulated_alignments --config example_files/config_file_simulated.txt --output example_files/snps/simulated_data --snps_per_locus one --phased --delimiter _`

By default the program extracts SNPs in binary format, where 0 = ancestral state, and 1 = derived state. In case of phased data the extracted SNPs are coded as 0 = hommozygous for ancestral state, 1 = heterozygous, 2 = homozygous for derived state (following convention of SNAPP). Instead we can also extract the actual nucleotides at the variable positions by adding `--export_nucleotides` to the command. This may be useful for other applications (e.g. for use in BEAST):

`python bin/snps_from_uce_alignments.py --input example_files/alignments/simulated_alignments --config example_files/config_file_simulated.txt --output example_files/snps/simulated_data --snps_per_locus one --phased --delimiter _ --export_nucleotides`


### Real, messy alignments

Empirical data often have additional complications, which were not present in the simulated alignments above. These are mostly missing sequences or parts of sequences for several samples at a given locus, as well as ambiguity codes at sites where there is uncertainty about the nucleotide call.

The program only considers the base DNA nucleotides (`A,C,T.G`) as well as gaps (`-`) as valid characters. When finding suitable positions for SNP extraction, these positions have to show exactly two variants among these 5 allowed characters, each of which has to be present in more than two sequences to avoid phylogeneticlaly uninformative sites. By default positions that have ambiguities of any form for any of the samples (e.g. `?,N,M,R,W,S,Y,K,V,H,D,B` or any other letters) are ignored and will not be extracted:

`python bin/snps_from_uce_alignments.py --input example_files/alignments/empirical_alignments --config example_files/config_file_empirical.txt --output example_files/snps/empirical_data --snps_per_locus one --phased --delimiter _`

In this case, the command leads to no positions being extracted, which becomes increasingly likely for alignments with many samples, with some of the samples missing in each alignment. To still extract SNP data, we can add the `--include_missing` flag, which allows for the extraction of sites even if these contain missing data for some of the samples.

`python bin/snps_from_uce_alignments.py --input example_files/alignments/empirical_alignments --config example_files/config_file_empirical.txt --output example_files/snps/empirical_data --snps_per_locus one --phased --delimiter _ --include_missing`

The config file makes it easy to select specific subsets of the taxa in your data, e.g. to only extract SNPs for a given clade or only for those taxa with good data. Here I selected the taxa with the most sequence information throughout the alignments and stored them in a new config file at `example_files/config_file_empirical_reduced.txt`

Focusing on only these taxa allows us to extract a fair number of high-quality SNPs, even when excluding sites with missing data:

`python bin/snps_from_uce_alignments.py --input example_files/alignments/empirical_alignments --config example_files/config_file_empirical_reduced.txt --output example_files/snps/empirical_data_reduced --snps_per_locus one --phased --delimiter _`
