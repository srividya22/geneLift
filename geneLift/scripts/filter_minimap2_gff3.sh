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
min_args=6

helpFunction()
{
   echo ""
   echo "Usage: $0 <s_path> <gff_file> <coverage> <identity> <path to output dir>"
   echo -e "\t - Path to scripts "
   echo -e "\t - paf file to filter alignments"
   echo -e "\t - Genome fasta file"
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

#check if user specified the requried parameters
if [ $# -eq $min_args ]; then
   s_path=$1
   paf_file=$2
   ref=$3
   cov=$4
   idty=$5
   out=$6

   COVERAGE=${cov}
   IDENTITY=${idty}
   mkdir -p ${out}

   # Code to priorotize gene ids
   sort -k6,6 -V -k8,8n ${paf_file} > ${out}/aln_minimap2_sorted.paf 

   # Filter alignments
   gawk '{s=$1; sub(/.[0-9]$/, "", s) ; split($1,n,"."); c[$1]+= 1 ; idty=($10/$2 * 100 ) ;OFS="\t" ; $1=s"-m"c[$1]"."n[length(n)] ; print $0"\t"idty }' ${out}/aln_minimap2_sorted.paf | awk -v x=${IDTY} '$NF >= x { OFS="\t" ;print $0}' > ${out}/aln_minimap2_sorted_filtered.paf

   #Convert paf to bed12 with multiple mappings
   paftools.js splice2bed -m ${out}/aln_minimap2_sorted_filtered.paf | awk '!seen[$0]++' > ${out}/aln_minimap2_sorted_filtered.bed12

   #Filter overlapping coords from bed12
   bedtools merge -i ${out}/aln_minimap2_sorted_filtered.bed12 -s -c 4 -o count  | awk 'NR==FNR { a[$1":"$2":"$3":"$4] ; next } ( ($1":"$2":"$3":"$6)  in  a) {OFS="\t"; print $0 }' - ${out}/aln_minimap2_sorted_filtered.bed12 > ${out}/aln_minimap2_non_overlapping.bed12

   #Convert bed12 to gff3
   ${s_path}/../bed_to_gff_converter.py ${out}/aln_minimap2_non_overlapping.bed12 ${out}/aln_minimap2.gff

   #Sort the gff3
   gt gff3 -sort -retainids yes ${out}/aln_minimap2.gff > ${out}/aln_minimap2_sorted.gff

   # Add CDS
   gt cds -seqfile ${ref} -matchdesc yes ${out}/aln_minimap2_sorted.gff > ${out}/aln_minimap2_cds_sorted.gff

   cd ${out}
   # Clean GFF to include IDS for everything
   GFFcleaner --add-missing-ids aln_minimap2_cds_sorted.gff > aln_minimap2_cds_clean.out
   cd - 
   if [ ! -f "${out}/aln_minimap2_cds_sorted.gff" ]; then
  	echo "ERROR: Filtering Minimap2 gene models !"
	exit
   fi
   if [ ! -f "${out}/aln_minimap2_cds_sorted_clean.gff" ]; then
  	echo "ERROR: GFFCleaner Cleaning up  Minimap2 gene models !"
	exit
   fi
   
else
  helpFunction 
fi
