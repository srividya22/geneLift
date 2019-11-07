#!/bin/bash
#submit_script.sh
                           # The following are options passed to qsub
#$ -V                      # Inherit the submission environment
#$ -cwd                    # Start job in submission directory
#$ -e $JOB_NAME.e$JOB_ID   # Combine stderr and stdout
#$ -o $JOB_NAME.o$JOB_ID   # Name of the output file (eg. myMPI.oJobID)
#$ -m bes                  # Email at Begin and End of job or if suspended
#$ -l m_mem_free=4g

## 
#Script to remover gene ids from gff based on validation report

# Set path to gff_valid.py script
EXE=${HOME}/bin

gff_file=${1} 
fasta=$2
report=$3
outdir=$4

base=$( basename ${gff_file} ${gff_file##*.} | cut -f1 -d"." )

# Count mRNAs and genes in the gff_file

total_mRNAs=$( awk '$3 =="mRNA" { sum += 1 }END  { print sum }' ${gff_file} )
total_genes=$( awk '$3 =="gene" { sum += 1 }END  { print sum }' ${gff_file} )

# Extract gene ids to remove genes
grep -B 1 "N_COUNT" ${report} | grep -o "ID=CDS:[A-Z,a-z,0-9,.,-:]*" | awk -F":" '{ print $2}' | sed -n 's/\([^"]]*\).[0-9]$/\1/p' > ${outdir}/to_be_removed_genes.txt

count=$( awk '{ sum += 1 }  END { print sum}' ${outdir}/to_be_removed_genes.txt )

# Remove selected genes from the files 
grep -vFf ${outdir}/to_be_removed_genes.txt ${gff_file} > ${outdir}/${base}_filt.gff

# Count final 

f_mRNAs=$( awk '$3 =="mRNA" { sum += 1 }END  { print sum }' ${outdir}/${base}_filt.gff )
f_genes=$( awk '$3 =="gene" { sum += 1 }END  { print sum }' ${outdir}/${base}_filt.gff )

# Rerun gff_valid.py on the new gff
${EXE}/gff_valid.py -g ${outdir}/${base}_filt.gff -f ${fasta} -t CDS -r ${outdir}/validation_report.txt
 
echo "Total genes : ${total_genes} Total mRNAS : ${total_mRNAs}" > ${outdir}/final_report.txt
echo "Filtered genes : ${count} " >> ${outdir}/final_report.txt
echo "Final genes : ${f_genes} Final mRNAS : ${f_mRNAs}" >> ${outdir}/final_report.txt
echo "DONE" >> ${outdir}/final_report.txt
