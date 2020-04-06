# Scripts and pipelines from Pathogenic allodiploid hybrids of <i>Aspergillus</i> fungi

This repository provides command-line arguments and script for reproducibility of Pathogenic allodiploid hybrids of <i>Aspergillus</i> fungi in accordance with STAR methods of Current Biology.

If you use contents from this repository, please cite:

<br /><br />

### Determination of gene set completeness
Determines the number of complete, duplicated, fragmented, and missing BUSCO genes in a proteome.<br />
Software: [BUSCO](https://busco.ezlab.org/)<br />
```
bash run_busco.sh proteome.faa
```
<br />

### Determination of copy number variable regions
Exemplary Control-FREEC configuration file. Calculating p-values of each copy number variable region was assessed using the original software author's script [assess_significance.R](http://boevalab.inf.ethz.ch/FREEC/)<br />
Software: [Control-FREEC](http://boevalab.inf.ethz.ch/FREEC/)<br />
```
freec -conf config.freec |& tee freec.log
```
<br />
In addition to using Control-FREEC, we assessed if CNVnator was an accurate software for our case. To do so, we used the following wrapper script to run CNVnator. 1st agrument should be the root file; 2nd argument should be the tree file; 3rd argument should be the genome file; 4th argument should be the scaffolds directory path; 5th argument should be the window size.<br />
Software: [CNVnator](https://github.com/abyzovlab/CNVnator)<br />
```
bash CNVnator_wrapper.sh file.root ./path_to_sorted_bam_file ./path_to_file_with_scaffold_lengths ./path_to_directory_with_scaffolds window_size
```
<br />

### Predicting gene boundaries 
Augustus is a powerful and popular tool for predicting gene boundaries. We ran augustus with default parameters with species training on <i>Aspergillus nidulans</i>.<br />
Software: [Augustus](http://augustus.gobics.de/)<br />
<br />
```
augustus --species=aspergillus_nidulans
```
<br />

### Genome assembly 
To simply the genome assembly process, we used the wrapper utility iWGS. iWGS was run with the recommended default parameters using Kmergenie, Trimmomatic, SPAdes, MaSuRCA, and QUAST.<br />
Software: [iWGS](https://github.com/zhouxiaofan1983/iWGS)<br />
<br />

### Strain determination using taxonomically informative loci
To determine the evolutionary history of taxonomically informative loci -- i.e., determine that the parents of the hybrids are <i>Aspergillus spinulosporus</i> and a close relative of <i>Aspergillus quadrilineatus</i> -- we used the maximum likelihood software RAxML. To determine bipartition support, we used rapid bootstrap analysis. Below we provide an exemplary command used during tree search.<br /> 
Software: [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/)<br />

```
raxmlHPC -f a -m GTRGAMMAX -x 12345 -p 12345 -N 1000
```
<br />

### Predicting orthologous groups of genes
To predict groups of orthologous genes for downstream phylogenetic analyses, we used a sequence similarity-based cluster approach. The following is an exemplary command of how we did so.<br />
Sofware: [OrthoFinder](https://github.com/davidemms/OrthoFinder/releases)<br />
```
./orthofinder -os -M msa -I 1.5 -S blast -f directory_of_proteomes/
```
<br />

### Sequence alignment and trimming
To align and trim sequences for downstream analysis, we first created nucleotide multi-fasta files of single copy orthologous genes predicted by OrthoFinder. Then, we aligned and trimmed the sequences in the multi-fasta files using the commands shown here.<br />
Software: [Mafft](https://mafft.cbrc.jp/alignment/software/), [trimAl](http://trimal.cgenomics.org/)<br />
```
# Sequence alignment
mafft --maxiterate 1000 --genafpair input > output 

# Alignment trimming
trimal -in input -out output -gappyout
```
<br />

### Creating a concatenated genome-scale data matrix of multiple sequence alignments
A custom script ([link](https://raw.githubusercontent.com/JLSteenwyk/Phylogenetic_scripts/master/create_concat_matrix.py)) was used to concatenated the aligned and trimmed sequences described in section 'Sequence alignment and trimming'. The input files and parameters include a list of alignment files to concatenate, a list of taxa to include, whether the sequences are proteins or nucleotides, and a prefix for output files. Output files include a fasta file of concatenated sequence with '.fa' appended to the end, a RAxML style partition file, and a file that summarizes the occupancy of each gene from each alignment. 
Software: [biopython](https://biopython.org/)
```
python create_concat_matrix.py -a alignment.list -c sequence_character -t taxa.list -p output_prefix
```
Original author: [Jacob Steenwyk](https://jlsteenwyk.github.io/)
<br /><br />

### Genome-scale phylogenies of each parental genome and topology tests
Genome-scale phylogenies to predict the evolutionary history of each subgenome were examined using the exemplary command described here.<br />
Software: [IQ-TREE](http://www.iqtree.org/)<br />
```
iqtree -s input.fa -seed 86924356 -st DNA -pre output -nt 24 -nbest 10 -m TEST -bb 5000
```
In addition, we conducted topology tests. In brief, these topology tests were used to determine if the topology inferred from one data matrix (from one parental genome) was equivalent to the topology inferred from the other data matrix (the other parental genome).
```
iqtree -s data_matrix_from_one_parent.fa -z Phylogenies_inferred_from_both_data_matrices.tres -n 0 -zb 10000 -zw -au -m GTR+F+I+G4
```
<br />

### Reciprocal best blast hit (RBBH)
To conduct reciprocal best blast analysis, we used the following custom script. Part of the script relies on a resource from Harvard ([link](http://archive.sysbio.harvard.edu/csb/resources/computational/scriptome/UNIX/Protocols/Sequences.html)). The exemplary script was used for RBBH between nucleotide sequences. Changing lines 13,14,17,18 can allow for RBBH between protein sequences.<br />
Sofware: [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [perl](https://www.perl.org/)<br />
```
bash RBBH.bash fasta_file_a fasta_file_b
```
<br />

### Ks values for every gene
A custom pipeline to calculate Ks for every gene in one genome compared to its best blast hit in another genome. See ./Ks_pipeline/README for a detailed explanation of the concept and usage. Note, this pipeline is explicitly designed to function with the slurm job scheduler at Vanderbilt University's high performance computing cluster, ACCRE (Advanced Computing Center for Research and Education).<br />
Sofware: [Blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download), [paml](http://abacus.gene.ucl.ac.uk/software/paml.html), [pal2nal](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1538804/), [samtools](http://www.htslib.org/download/), [perl](https://www.perl.org/)<br />
```
bash ./Ks_pipeline/execute.sh
```




