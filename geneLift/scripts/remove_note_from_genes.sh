#!/usr/bin/env bash
# Script to modify the exon id format in gff to exon:<mRNAID>:count , CDS:<mRNAID>:count , five_prime_UTR:<mRNAID>:count

gff=${1}
out=${2}

awk -F"\t" '{if ( $3 == "gene" )
        {
        gsub(/Alias=([A-Za-z0-9\._:,\-]+);/,"",$9) ; 
        gsub(/ref_id=([A-Za-z0-9\._:,\-]+);/,"",$9) ; 
        gsub(/;Note=([A-Za-z0-9\._:, \-]+)/,"",$9) ; 
        OFS="\t" ; print }
     else if ( $3 == "mRNA" )
        {
        gsub(/Alias=([A-Za-z0-9\._:,\-]+);/,"",$9) ; 
        gsub(/ref_id=([A-Za-z0-9\._:,\-]+);/,"",$9) ; 
        OFS="\t" ; print }
     else { OFS="\t" ; print } }' ${gff} > ${out}
