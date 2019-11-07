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
   awk -F"\t" '{  if ( !/^#/ ) { $2="geneLiftOver" } ; OFS="\t" ; print }' ${gff_file} > ${out}/aln_minimap2_cds_UTRS_sorted_modified.gff

   # Extract and Rename Ids
   awk -F"\t" '!/^#/ { split($9,a,";") ;  if ($3 == "gene" ) { print a[1]"\t"a[2] } if ( $3 == "mRNA") { print a[1]"\t"a[3] } }' ${out}/aln_minimap2_cds_UTRS_sorted_modified.gff  | sed 's/[A-Z,a-z]*=//g'> ${out}/minimap2_fids_mapped_ids.txt

   # Rename ids
   ${s_path}/map_gff_ids.pl ${out}/minimap2_fids_mapped_ids.txt ${out}/aln_minimap2_cds_UTRS_sorted_modified.gff

   # Rename child feature ids
    ${s_path}/format_child_ids.sh ${out}/aln_minimap2_cds_UTRS_sorted_modified.gff ${out}/aln_minimap2_cds_UTRS_sorted_modified_ids.gff

    gawk '$3=="mRNA" { split($9,a,";"); split(a[2],gid,"="); split(a[1],b,"="); split(b[2],c,"-") ; split(b[2],d,".") ; print c[1]"."d[length(d)]";"gid[2]";"b[2] }' ${out}/aln_minimap2_cds_UTRS_sorted_modified_ids.gff > ${out}/minimap2_mRNAs.txt

    awk -F";" 'NR==FNR { a[$1] = $3 ; next } ($1 in a ) { print $3";"$2";"a[$1] } ' ${annot} ${out}/minimap2_mRNAs.txt > ${out}/minimap2_mapped_func.txt 


   if [ ! -f "${out}/minimap2_mapped_ids.txt" ]; then
  	echo "ERROR: Formatting gene models !"
	exit
   fi
   if [ ! -f "${out}/minimap2_mapped_func.txt" ]; then
  	echo "ERROR: Formatting gene models !"
	exit
   fi
else
  helpFunction 
fi
