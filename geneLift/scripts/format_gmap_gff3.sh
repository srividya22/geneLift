#!/usr/bin/env bash 
############################################################
#
# Script to filter gff3 from gmap to remove overlapping genes
# Written by Srividya Ramakrishnan
# Affiliation : Johns Hopkins University 
#
###########################################################

set -e

# Minimum args to call the script
min_args=4

helpFunction()
{
   echo ""
   echo "Usage: $0 <gff_file> <red annotation txt> <path to output dir>"
   echo -e "\t - GFF file to format alignments"
   echo -e "\t - Reference annotation text file ( ; seperated )"
   echo -e "\t - Output dir path for output files"
   exit 1 # Exit script after printing help
}

trap 'catch $? $LINENO' EXIT
catch() {
  if [ "$1" != "0" ]; then
    # error handling goes here
    echo "Error $1 occurred on $2"
  fi
}


#check if user specified the requried parameters
if [ $# -eq $min_args ]; then
   s_path=$1
   gff_file=$2
   annot=$3
   out=$4

   mkdir -p ${out}
   ${s_path}/format_child_ids.sh ${gff_file} ${out}/gmap_format.gff
   gawk '$3=="mRNA" { split($9,a,";"); split(a[2],gid,"="); split(a[1],b,"="); split(b[2],c,"-") ; split(b[2],d,".") ; print c[1]"."d[length(d)]";"gid[2]";"b[2] }' ${out}/gmap_format.gff > ${out}/gmap_mRNAs.txt
   awk -F";" 'NR==FNR { a[$1] = $3 ; next } ($1 in a ) { print $3";"$2";"a[$1] } ' ${annot} ${out}/gmap_mRNAs.txt > ${out}/gmap_mapped_func.txt
   if [ ! -f "${out}/gmap_mapped_func.txt" ]; then
  	echo "ERROR: Formatting gene models !"
	exit
   fi
else
  helpFunction 
fi
