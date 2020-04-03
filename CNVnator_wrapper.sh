###########
# This script will run CNVnator for you using arguments of 
# a root file, tree file, genome file, scaffolds directory path, and window size
#
# 1st agrument should be the root file
# 2nd argument should be the tree file
# 3rd argument should be the genome file
# 4th argument should be the scaffolds directory path
# 5th argument should be the window size
#
###########

# save arguments
ROOT=$1
TREE=$2
GENOME=$3
SCAFFOLDS_DIR=$4
WS=$5

WD=$(pwd)
CNVnator=$(echo "/home/steenwj/SOFTWARE/CNVnator_v0.3.2/src/cnvnator")
OUT_FILE=$(echo "$TREE" | awk '{print $1".CNVnator"}')
OUT_FINAL=$(echo "$OUT_FILE.WS$WS")
echo $OUT_FINAL

$CNVnator -root $ROOT -tree $TREE -genome $GENOME -d $SCAFFOLDS_DIR
$CNVnator -root $ROOT -genome $GENOME -his $WS -d $SCAFFOLDS_DIR
$CNVnator -root $ROOT -genome $GENOME -stat $WS -d $SCAFFOLDS_DIR
$CNVnator -root $ROOT -genome $GENOME -partition $WS -d $SCAFFOLDS_DIR
$CNVnator -root $ROOT -genome $GENOME -call $WS -d $SCAFFOLDS_DIR > $OUT_FINAL
