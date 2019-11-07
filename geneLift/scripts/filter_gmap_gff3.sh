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
min_args=5

helpFunction()
{
   echo ""
   echo "Usage: $0 <gff_file> <coverage> <identity> <path to output dir>"
   echo -e "\t - GFF file to filter alignments"
   echo -e "\t - Coverage of cDNA alignments to be retained ; default :  90"
   echo -e "\t - Identity of cDNA alignments to be retained ; default : 95"
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

#Order Gene Models
GENE_PAT="TomatoPan" # pattern for less priority genes

#check if user specified the requried parameters
if [ $# -eq $min_args ]; then
   s_path=$1
   gff_file=$2
   cov=$3
   idty=$4
   out=$5

   COVERAGE=${cov}
   IDENTITY=${idty}
   mkdir -p ${out}

   # Code to priorotize gene ids
   #awk '$3 == "mRNA" { print $0 }' ${gff_file} | sed 's/[;,=]/\t/g' | awk '{ print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$16"\t"$18"\t"$20 }' |sort -k9,9nr -k7,7nr -k8,8nr | grep -v ${GENE_PAT} > ${out}/priority1_ids.txt

   #awk '$3 == "mRNA" { print $0 }' ${gff_file} | sed 's/[;,=]/\t/g' | awk '{ print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$16"\t"$18"\t"$20 }' |sort -k9,9nr -k7,7nr -k8,8nr | grep ${GENE_PAT} > ${out}/priority2_ids.txt

   #cat ${out}/priority1_ids.txt ${out}/priority2_ids.txt > ${out}/all_mRNAs.txt
    
   awk '$3 == "mRNA" { print $0 }' ${gff_file} | sed 's/[;,=]/\t/g' | awk '{ print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$16"\t"$18"\t"$20 }' |sort -k9,9nr -k7,7nr -k8,8nr > ${out}/gmap_all_mRNAs.txt
   
   ${s_path}/../remove_overlapping_genes.py ${out}/gmap_all_mRNAs.txt 90 ${out}/gmap_non_overlapping_mRNA.txt
   
   if [ ! -f "${out}/gmap_non_overlapping_mRNA.txt" ]; then
  	echo "ERROR: Removing overlapping gene models !"
	exit
   fi
   grep -Ff ${out}/gmap_non_overlapping_mRNA.txt ${gff_file} | sed '/###$/d' | sed 's/;Target=[A-Z,a-z,0-9, ,.,+,-]*//g' |  sed -E 's/.([0-9]).path([0-9])/-gm\2/g' | sed -E 's/.([0-9]).mrna([0-9])/-gm\2.\1/g' | sed '1i ##gff-version 3' > ${out}/gmap_filtered_nonoverlapping.gff
   
   if [ ! -f "${out}/gmap_filtered_nonoverlapping.gff" ]; then
  	echo "ERROR: Filtering overlapping gene models !"
	exit
   fi
else
  helpFunction 
fi
