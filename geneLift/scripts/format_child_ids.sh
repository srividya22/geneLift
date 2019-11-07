#!/usr/bin/env bash
# Script to modify the exon id format in gff to exon:<mRNAID>:count , CDS:<mRNAID>:count , five_prime_UTR:<mRNAID>:count

gff=${1}
out=${2}

awk -F"\t" '/^#/ { print } 
  { if( $3 == "mRNA" ) { gsub(/Name=[A-Za-z0-9\.:_\-]+;/,"Name="id";", $9); OFS="\t"; print } 
  if ( $3 == "gene" )
     {
        match($9, /ID=([A-Za-z0-9\._:\-]+)/,m) ; id = m[1] ;
        gsub(/Name=([A-Za-z0-9\.:_\-]+)/,"Name="id";", $9); 
        OFS="\t" ; print }
  if ( $3 == "exon" ) 
     { 
        match($9, /Parent=([A-Za-z0-9\._:\-]+)/,m) ; mid = m[1] ; ec[mid] += 1 ; 
        gsub(/ID=[A-Za-z0-9\.:_\-]+;/,"ID=exon:"mid":"ec[mid]";", $9); gsub(/Name=[A-Za-z0-9\.:_\-]+;/,"", $9); OFS="\t" ; print } 
  if ( $3 == "CDS" ) 
      { 
        match($9, /Parent=([A-Za-z0-9\._:\-]+)/,m) ; cdd = m[1] ; cc[cdd] += 1 ; 
       gsub(/ID=[A-Za-z0-9\.:_\-]+;/,"ID=CDS:"cdd":"cc[cdd]";", $9); gsub(/Name=[A-Za-z0-9\.:_\-]+;/,"", $9); OFS="\t" ; print } 
  if ( $3 == "five_prime_UTR" ) 
      { 
        match($9, /Parent=([A-Za-z0-9\._:\-]+)/,m) ; fpd = m[1] ; fpc[fpd] += 1 ; 
        gsub(/ID=[A-Za-z0-9\.:_\-]+;/,"ID=five_prime_UTR:"fpd":"fpc[fpd]";", $9); OFS="\t" ; print } 
  if ( $3 == "three_prime_UTR" ) 
      { 
        match($9, /Parent=([A-Za-z0-9\._:\-]+)/,m) ; tpd = m[1] ; tpc[tpd] += 1 ; 
        gsub(/ID=[A-Za-z0-9\.:_\-]+;/,"ID=three_prime_UTR:"tpd":"tpc[tpd]";", $9); OFS="\t" ; print }}' ${gff} > ${out}
