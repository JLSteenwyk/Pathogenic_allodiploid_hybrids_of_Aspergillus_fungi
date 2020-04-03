# usage: bash script.bash spp1_short_name spp1.cds spp1.pep spp2_short_name spp2.cds spp2.pep
# bash determine_Ks_for_every_gene_between_two_species.bash Tcaseinolytica tortispora_caseinolytica.max.cds tortispora_caseinolytica.max.pep Tganteri yHMPu5000035654_tortispora_ganteri_160519.max.cds yHMPu5000035654_tortispora_ganteri_160519.max.pep

# store arguments as variables
spp1_short_name=$1
spp1_cds=$2
spp1_pep=$3
spp2_short_name=$4
spp2_cds=$5
spp2_pep=$6

# parse blast for BBH
echo -e "Parsing blast for best blast hit\n..."
for i in $(ls *fa.fmt6.blast); do sort -k12,12 -n $i -r | sort -k1,1 -u; done > best_hit_per_gene.txt
echo -e "PASS"

# create mfasta for pep and cds with spp1 gene BBH
echo -e "Creating mfastas for pep with spp1 BBH to spp2\n..."
bash ./scripts/parse_BBHs_and_create_mfasta_pep.bash best_hit_per_gene.txt $spp1_pep.MOD $spp2_pep.MOD
echo "PASS"
echo -e "Creating mfastas for cds with spp1 BBH to spp2\n..."
bash ./scripts/parse_BBHs_and_create_mfasta_cds.bash best_hit_per_gene.txt $spp1_cds.MOD $spp2_cds.MOD
echo "PASS"

# align protein fasta files
echo -e "Aligning protein mfasta files using mafft\n..."
for i in $(ls -l ./FILES_BBH_mfastas_pep/*mfasta | awk '{print $9}'); do mafft --reorder --bl 62 --op 1.0 --maxiterate 1000 --retree 1 --genafpair $i > $i.mafft; done
echo -e "PASS"

# make protein fasta list
echo -e "Creating protein fasta list\n..."
ls -l ./FILES_BBH_mfastas_pep/*mafft | awk '{print $9}' > protein_fasta.list
echo "PASS"

# cds fasta list
echo -e "Creating cds fasta list\n..."
ls -l ./FILES_BBH_mfastas_cds/*fasta | awk '{print $9}' > cds_fasta.list
echo "PASS"

# reorder cds list to match protein list
echo -e "Reordering cds list to match protein list\n..."
cat protein_fasta.list | sed 's/.\/FILES_BBH_mfastas_pep\///g' | sed 's/_1././g' | sed 's/.mafft//g' | xargs -Ihello grep -w hello cds_fasta.list > cds_fasta.list.reordered
echo "PASS"

# execute pal2nal
echo -e "Executing pal2nal for protein list and cds list"
bash ./scripts/execute_pal2nal.bash protein_fasta.list cds_fasta.list.reordered 2>> pal2nal_errors.txt
echo "PASS"

# remove failed pal2nal conversions
echo -e "removing failed pal2nal conversions\n..."
for i in $(ls -l ./FILES_BBH_mfastas_pep/*aln_nuc.fasta | awk '{if ($5=="0") print $9}'); do rm $i; done
echo "PASS"

# create yn00 ctl files
echo -e "Creating yn00 ctl files\n..."
for i in $(ls ./FILES_BBH_mfastas_pep/*aln_nuc.fasta | sed 's/.\/FILES_BBH_mfastas_pep\///g'); do echo $i; cat ./scripts/yn00.ctl | sed "s|yn.test|$i.YN.OUT|g" | sed "s|test|./FILES_BBH_mfastas_pep/$i|g" > $i.yn00.ctl; echo "..."; done
echo "PASS"

# calculate Ks
echo -e "Calculating Ks\n..."
for i in $(ls *mafft.aln_nuc.fasta.yn00.ctl); do echo "determing Ks using ctl file $i"; yn00 $i; echo "..."; done
echo "PASS"

# parse yn00 output files
echo -e "Parsing yn00 output files\n..."
for i in $(ls *YN.OUT); do echo -n -e "$i\t" | sed 's/.mfasta.mafft.aln_nuc.fasta.YN.OUT//g'; cat $i | grep LWL85m: | awk '{print $4}' | sed 's/nan/0/g' | sed 's/-[0-9][0-9]*/0/' ; echo -e "\n"; done | sed '/^$/d' | awk -F'\t' '$2!=""' | awk -F'\t' '$2!="inf"' > Ks_val_per_gene.txt
echo "PASS"

# cleaning
echo -e "Moving single gene fasta files\n..."
mv $spp1_short_name@*.fa FILES_single_gene_fasta
echo "PASS"

echo -e "Moving single gene blast results\n..."
mv $spp1_short_name@*.fa.fmt6.blast FILES_single_gene_blast
echo "PASS"

echo -e "Moving yn00 ctl files\n..."
mv *.mfasta.mafft.aln_nuc.fasta.yn00.ctl FILES_yn00_ctl
echo "PASS"

echo -e "Moving yn00 output files\n..."
mv *mfasta.mafft.aln_nuc.fasta.YN.OUT FILES_yn00_out
echo -e "PASS"

echo -e "Moving gene name files\n..."
mv *cds.gene_names *pep.gene_names FILES_gene_names
echo "PASS"

echo -e "Removing split files\n..."
rm x0*
echo "PASS"

echo -e "Moving slurm output files\n..."
mkdir FILES_slurm_outfiles
mv slurm* FILES_slurm_outfiles
echo "PASS"

echo -e "Removing yn00 extra files\n..."
rm 2YN* rst rst1 rub
echo "PASS"

echo -e "Moving fas files\n..."
mv *.MOD* FILES_fas
echo "PASS"

echo -e "Moving extra files\n..."
mkdir FILES_extras
mv *.out cds_fasta.list* best_hit_per_gene.txt *jobID* blast_every_gene.bash pal2nal_errors.txt protein_fasta.list FILES_extras
echo "PASS"

echo -e "Moving all FILES directories to intermediates directory\n..."
mkdir intermediates
mv FILES* intermediates
echo "PASS"

echo -e "Tar zipping intermediates directory and removing intermediates directory"
tar -zcvf intermediates.tar.gz intermediates
rm -rf intermediates
echo "PASS"

