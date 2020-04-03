# usage: bash script.bash spp1_short_name spp1.cds spp1.pep spp2_short_name spp2.cds spp2.pep
# bash determine_Ks_for_every_gene_between_two_species.bash Tcaseinolytica tortispora_caseinolytica.max.cds tortispora_caseinolytica.max.pep Tganteri yHMPu5000035654_tortispora_ganteri_160519.max.cds yHMPu5000035654_tortispora_ganteri_160519.max.pep

# store arguments as variables
spp1_short_name=$1
spp1_cds=$2
spp1_pep=$3
spp2_short_name=$4
spp2_cds=$5
spp2_pep=$6

# create dirs to store results
mkdir FILES_fas
mkdir FILES_gene_names
mkdir FILES_aln_nuc
mkdir FILES_yn00_ctl
mkdir FILES_yn00_out
mkdir FILES_single_gene_blast
mkdir FILES_single_gene_fasta

### Rename fasta files
# check for conflicts between naming schemes
echo -e "checking if gene names follow the same order\ngenes are assumed to be in the same order and named using the same convention in each file\n..."

number_of_conflicts_spp1=$(paste <(cat $spp1_pep | grep '>' | awk '{print $2}') <(cat $spp1_cds | grep '>' | awk '{print $2}') | awk '{if ($2!=$2) print $0}' | wc -l)
number_of_conflicts_spp2=$(paste <(cat $spp2_pep | grep '>' | awk '{print $2}') <(cat $spp2_cds | grep '>' | awk '{print $2}') | awk '{if ($2!=$2) print $0}' | wc -l)
if [[ "$number_of_conflicts_spp1" != 0 ]]; then
echo "Number of gene conflicts is $number_of_conflicts_spp1"
echo -e "The gene order between the pep and cds file is not the same for $spp1_short_name\nFAIL"
exit
fi
if [[ "$number_of_conflicts_spp2" != 0 ]]; then
echo "Number of gene conflicts is $number_of_conflicts_spp1"
echo -e "The gene order between the pep and cds file is not the same for $spp2_short_name\nFAIL"
exit
fi
echo -e "PASS\n"

# rename using convention >shortname@ for spp1 cds
echo "Genes will be renamed using the following convention for $spp1_cds: $spp1_short_name@geneID"
awk -v shortname1="$spp1_short_name" '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp1_cds > $spp1_cds.MOD
echo "..."
awk '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp1_cds | grep '>' > temp.txt
echo "..."
grep '>' $spp1_cds | paste temp.txt - > $spp1_short_name.cds.gene_names
echo "..."
echo -e "PASS\n"

# rename using convention >shortname@ for spp1 pep
echo "Genes will be renamed using the following convention for $spp1_pep: $spp1_short_name@geneID"
awk -v shortname1="$spp1_short_name" '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp1_pep > $spp1_pep.MOD
echo "..."
awk '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp1_cds | grep '>' > temp.txt
echo "..."
grep '>' $spp1_cds | paste temp.txt - > $spp1_short_name.pep.gene_names
echo "..."
echo -e "PASS\n"

# rename using convention >shortname@ for spp2 cds
echo "Genes will be renamed using the following convention for $spp2_cds: $spp2_short_name@geneID"
awk -v shortname1="$spp2_short_name" '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp2_cds > $spp2_cds.MOD
echo "..."
awk '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp2_cds | grep '>' > temp.txt
echo "..."
grep '>' $spp2_cds | paste temp.txt - > $spp2_short_name.cds.gene_names
echo "..."
echo -e "PASS\n"

# rename using convention >shortname@ for spp2 pep
echo "Genes will be renamed using the following convention for $spp2_pep: $spp2_short_name@geneID"
awk -v shortname1="$spp2_short_name" '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp2_pep > $spp2_pep.MOD
echo "..."
awk '/^>/{print ">"shortname1"@" ++i; next}{print}' $spp2_cds | grep '>' > temp.txt
echo "..."
grep '>' $spp2_cds | paste temp.txt - > $spp2_short_name.pep.gene_names
echo "..."
echo -e "PASS\n"

# clean
rm temp.txt

### Finished renaming files to new convention


### blast every gene of spp1 against the genome of spp2
# split spp1 fasta into separate gene fastas
perl ./scripts/split_fasta_file_by_header.pl $spp1_pep.MOD

# make blastdb from spp2 pep
echo -e "Creating blastDB from $spp2_pep\n..."
makeblastdb -in $spp2_pep.MOD -dbtype prot
echo -e "PASS\n"

# create script to blast every gene
ls -l $spp1_short_name@*.fa | awk -v blastDB="$spp2_pep" '{print "blastp -query ", $9, " -db " blastDB".MOD -evalue 10 -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore\" -num_threads 12 -out", $9".fmt6.blast" }' > blast_every_gene.bash
echo -e "Creating script to blast every gene to send as a job array\n..."
total_lines=$(wc -l blast_every_gene.bash | awk '{print $1}')
echo "$total_lines in blast_every_gene.bash"
lines_per_part=$(echo "$total_lines" | awk '{print int($1/9)-1}')
echo "splitting blast_every_gene.bash into 10 equal parts with ~$lines_per_part genes per blast"
echo "..."
split -d -l $lines_per_part blast_every_gene.bash
echo "PASS"


