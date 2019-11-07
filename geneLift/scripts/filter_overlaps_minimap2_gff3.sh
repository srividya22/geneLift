#!/usr/bin/env sh 

gff_file=$1
cov=$2
out=$3

COVERAGE=${cov}
IDENTITY=95

GENE_PAT="TomatoPan" # pattern for less priority genes

mkdir -p ${out}


#awk '$3 == "mRNA" { print $0 }' ${gff_file} | sed 's/[;,=]/\t/g' | awk '{ print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$16"\t"$16"\t"$16 }' | sort -k7,7nr | grep -v ${GENE_PAT} > ${out}/priority1_ids.txt

#awk '$3 == "mRNA" { print $0 }' ${gff_file} | sed 's/[;,=]/\t/g' | awk '{ print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$16"\t"$16"\t"$16 }' | sort -k7,7nr | grep ${GENE_PAT} > ${out}/priority2_ids.txt


awk '$3 == "mRNA" { print $0 }' ${gff_file} | sed 's/[;,=]/\t/g' | awk '{ print $1"\t"$4"\t"$5"\t"$7"\t"$10"\t"$14"\t"$16"\t"$16"\t"$16 }' | sort -k7,7nr > ${out}/all_mRNAs.bed

# Change it to the python code
/seq/schatz/sramakri/sources/geneLift/geneLift/remove_overlapping_genes.py ${out}/all_mRNAs.bed 90 ${out}/non_overlapping_mRNA.txt


#/seq/schatz/sramakri/sources/geneLift/geneLift/remove_overlapping_with_mRNA.pl ${out}/mRNA_gmap_filtered.gff3 ${out}/non_overlapping_mRNA.txt

grep -Ff ${out}/non_overlapping_mRNA.txt ${gff_file} | sed '/###$/d' | sed '1i ##gff-version 3' > ${out}/final_nonoverlapping_minimap2.gff3

#awk '$3=="gene"{ print $9 }' ${final_out} | awk -F";" '{ print $1}' |sed 's/.path[0-9];//g' | sed 's/[A-Z,a-z]*=[a-z]*://g' | awk -F"." '{ print $1}'  > final_all_gmap_annotated.txt

#awk 'NR==FNR { a[$1]=$2 ; next } ($1 in a ) { print $1"\t"a[$1] }' /isilon/seq/schatz/sramakri/assembly/vole_genomes/broad_assembly/protein_coding_genes.txt final_all_gmap_annotated.txt | awk '{ print $2}' | sort | uniq -c | wc -l

#sed '/###$/d' ${final_out} | sed 's/.mrna[0-9]*//g'  | sed 's/.[0-9].path[0-9,.]*//g' > ${final_out}
