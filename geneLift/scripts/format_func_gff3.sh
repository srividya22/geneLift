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
min_args=3

helpFunction()
{
   echo ""
   echo "Usage: $0 <gff_file> <path to output dir>"
   echo -e "\t - GFF file to format alignments"
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
   out=$3

   mkdir -p ${out}
   sed -i '/^;$/d' ${gff_file}
   awk '!seen[$0]++' ${gff_file} > ${out}/gmap_func_dedup.gff
   ${s_path}/remove_note_from_genes.sh ${out}/gmap_func_dedup.gff ${out}/gmap_gene_models.gff
   if [ ! -f "${out}/gmap_gene_models.gff" ]; then
  	echo "ERROR: Formatting functional gene models !"
	exit
   fi
else
  helpFunction 
fi
