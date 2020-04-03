# Scripts and pipelines from Pathogenic allodiploid hybrids of <i>Aspergillus</i> fungi

This repository provides command-line arguments and script for reproducibility of Pathogenic allodiploid hybrids of <i>Aspergillus</i> fungi in accordance with STAR methods of Current Biology.

If you use contents from this repository, please cite:


### BUSCO
Determines the number of complete, duplicated, fragmented, and missing BUSCO genes in a proteome.<br />
```
bash run_busco.sh proteome.faa
```

### Control-FREEC
Exemplary Control-FREEC configuration file. Calculating p-values of each copy number variable region was assessed using the original software author's script assess_significance.R (http://boevalab.inf.ethz.ch/FREEC/)<br />
```
freec -conf config.freec |& tee freec.log
```

### CNVnator
Wrapper script to run CNVnator. 1st agrument should be the root file; 2nd argument should be the tree file; 3rd argument should be the genome file; 4th argument should be the scaffolds directory path; 5th argument should be the window size<br />
```
bash CNVnator_wrapper.sh file.root ./path_to_sorted_bam_file ./path_to_file_with_scaffold_lengths ./path_to_directory_with_scaffolds window_size
```

### AUGUSTUS PLACE HOLDER 
<br />
```
```

### iWGS PLACE HOLDER 
<br />
```
```

### RAxML PLACE HOLDER 
<br />
```
```

### IQ-TREE 
Genome-scale phylogenies to predict the evolutionary history of each subgenome were examined using the following exemplary command:<br />
```
iqtree -s input.fa -seed 86924356 -st DNA -pre output -nt 24 -nbest 10 -m TEST -bb 5000
```
The following is an exemplary command for the topology tests conducted. These topology tests were used to determine if the topology inferred from one data matrix (from one parental genome) was equivalent to the topology inferred from the other data matrix (the other parental genome).
```
iqtree -s data_matrix_from_one_parent.fa -z Phylogenies_inferred_from_both_data_matrices.tres -n 0 -zb 10000 -zw -au -m GTR+F+I+G4
```

### Ks values for every gene 
A custom pipeline to calculate Ks for every gene in one genome compared to its best blast hit in another genome. See ./Ks_pipeline/README for a detailed explanation of the concept and usage.<br />
```
bash ./Ks_pipeline/execute.sh
```




