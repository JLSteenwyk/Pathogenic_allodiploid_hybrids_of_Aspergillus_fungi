#!/bin/bash

###
#usage bash script.bash BBH_output_file
# This script will create mfastas for the best blast hit (BBH) between single genes from jakobsenii (genome A) to 
# singularis (genome B) for later use in determining Ks between the two genes after alignment
###

###
# Input file:
# BBH_output_file is the output of each gene singly blasted against the db genome
# and then the output was sorted by best hit according to bitscore and sorted into
# the following summary file titled best_hit_per_gene.txt in directory 
# /data/steenwj/HANSENIASPORA/orthomcl_v1.4/May_17/orthomcl_mfastas/FILES_jakobensii_singular_Ks
# The first few lines of the file look like the following where the 
# first two columns are most important:
# Hjakobsenii_35692-160519@1000_1	Hsingularis_35696-160519@2466_1	32.973	370	244	2	6	373	6	373	1.93e-54	184
# Hjakobsenii_35692-160519@1001_1	Hsingularis_35696-160519@3793_1	33.333	27	18	0	38	64	298	324	10.0	21.9
# Hjakobsenii_35692-160519@100_1	Hsingularis_35696-160519@1178_1	75.000	348	87	0	1	348	30	377	0.0	536
# ...					...			...	...	..	..	..	..	..	..	..	..	..	..
#
# The blast command used was:
# blastp -query  ./FILES_Hjakobsenii_single_gene_fasta/Hjakobsenii_35692-160519@1000_1.fa  -db yHMPu5000035696_hanseniaspora_singularis_160519.max_MOD.pep.fasta \
# -evalue 10 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 12 -out \
# ./FILES_Hjakobsenii_single_gene_fasta/Hjakobsenii_35692-160519@1000_1.fa.fmt6.blast
###

# save BBH argument as BBH_file; save two pep fasta files for genome A (jakobsenii) and genome B (singularis)
BBH_file=$1
A_pep_fasta=$2
B_pep_fasta=$3

# create subdirectory for multi-fasta files
mkdir FILES_BBH_mfastas_pep

# open BBH file and print out the first two columns and remove protein tag "_1"
cat $BBH_file | awk '{print $1, $2}' | \

# loop through lines BBH output
while read line col
do
	echo "Creating mfasta file for $line"
	# loop through columns
	for column in $col
	do
		A_gene=$(echo "$line")
		B_gene=$(echo "$column")
		samtools faidx $A_pep_fasta $A_gene >> ./FILES_BBH_mfastas_pep/$A_gene.mfasta
		samtools faidx $B_pep_fasta $B_gene >> ./FILES_BBH_mfastas_pep/$A_gene.mfasta
	done
	echo "..."
done < /dev/stdin
