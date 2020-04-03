# usage: bash script.bash list_of_aligned_protein_fastas list_of_nucleotide_fasta_files

# save arguments
alnPfasta=$1
nuclfasta=$2

# loop through both list files at once
paste $alnPfasta $nuclfasta | \
while read Pfasta Nfasta
do
	echo "execute pal2nal on $Pfasta $Nfasta"
	echo "..."
	perl /scratch/steenwj/SOFTWARE/pal2nal.v14/pal2nal.pl $Pfasta $Nfasta -codontable 1 -nomismatch -output fasta > $Pfasta.aln_nuc.fasta
	echo "complete"
	echo ""
done

