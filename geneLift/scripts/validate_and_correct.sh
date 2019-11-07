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
   echo "Usage: $0 <gff_file> <fasta> <path to output dir>"
   echo -e "\t - GFF file to format alignments"
   echo -e "\t - Genome fasta file"
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
   ref=$3
   out=$4

   mkdir -p ${out}
   gff_valid.py -g ${gff_file} -f ${ref} -r ${out}/report.txt
   ${s_path}/remove_failed_genes.sh ${gff_file} ${ref} ${out}/report.txt ${out}
   if [ ! -f "${out}/gmap_func_final_filt.gff" ]; then
  	echo "ERROR: Validating gene models !"
	exit
   fi
else
  helpFunction 
fi
