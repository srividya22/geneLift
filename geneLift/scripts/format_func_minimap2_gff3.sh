#!/usr/bin/env bash 
############################################################
#
# Script to filter gff3 from minimap2 to remove overlapping genes
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
   echo "Usage: $0 <s_path> <gff_file> <cov> <path to output dir>"
   echo -e "\t - GFF file to format alignments"
   echo -e "\t - Coverage"
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
   cov=$3
   out=$4

   mkdir -p ${out}
   sed -i '/^;$/d' ${gff_file}
   #awk '!seen[$0]++' ${gff_file} > ${out}/minimap2_func_dedup.gff
   ${s_path}/filter_overlaps_minimap2_gff3.sh ${gff_file} ${cov} ${out}
   ${s_path}/remove_note_from_genes.sh ${out}/final_nonoverlapping_minimap2.gff3 ${out}/minimap2_gene_models.gff
   if [ ! -f "${out}/minimap2_gene_models.gff" ]; then
  	echo "ERROR: Formatting functional gene models !"
	exit
   fi
else
  helpFunction 
fi
