#!/usr/bin/env python 
"""
####################
geneLiftover.py

Script to liftover reference gene models into a new assembly using cDNA alignments

Author :  Srividya Ramakrishnan
Affliation : Johns Hopkins University

####################

Usage: 
geneLift.py [options] -cDNA <FILE> -g <FILE> -aligner [gmap, minimap 2] -o <folder> 

-cDNA      # cDNA fasta with alignments
-g         #  Genome fasta
-aligner   # Aligner of choice either gmap or minimap2 ; gmap is default
-c         #  Coverage threshold for cDNA alignment
--report-duplicates  # Report gene duplications in the final gene models.
-i         #  Identity threshold for cDNA alignment
-o         # Output folder 
-t         # threads used for cDNA alignments
-h         # help          

"""

from __future__ import (division, print_function, absolute_import,unicode_literals)
import sys, re
import argparse
from geneLift import *
from geneLift.utils import *
from gmap_pipeline import Gmap
from minimap2_pipeline import Minimap2


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to liftover gene models from reference cDNA alignments")
    parser.add_argument("-cDNA",dest="cdna",type=str, required=True, help="fasta file containing transcript sequences")
    parser.add_argument("-g",dest="faa",type=str, required=True, help="Assembly fasta file to lift over gene models")
    parser.add_argument("-func",dest="annot",type=str, required=True, help="Semicolon (;) seperated text file with transcript_id,gene_id,functional annotation for the cDNA file")
    parser.add_argument("-o",dest="out_p", type=str, default="geneLift_output", help="output directory name")
    parser.add_argument("-x",dest="csp",action='store_false', help="enable cross species cDNA alignment ; only supported with aligner gmap")
    parser.add_argument("-aligner", dest="aligner",type=str, default="gmap", help='Aligner used for cDNA alignments options: [gmap,minimap2]')
    #parser.add_argument("--merge", dest="merge", action='store_false', help='Merge gmap and minimap2 gff files')
    parser.add_argument("--report-duplications", dest="dup",action='store_true', default=False, help="Report gene duplications")
    parser.add_argument("-mm", dest="m_path", type=str, default="minimap2", help='path to minimap2 executable')
    parser.add_argument("-gm", dest="gm_path", type=str, default="gmap", help='path to gmap executable')
    parser.add_argument("-c", metavar="90",dest="cov", type=int, default=90, help='Coverage threshold for transcript alignments')
    parser.add_argument("-i", metavar="95", dest="idty",type=int, default=95, help='Identity threshold for transcript alignments')
    parser.add_argument("-t", metavar="1", dest="threads",type=int, default=1, help='Number of threads')
    args = parser.parse_args()
    if args.aligner == "minimap2": 
          # Call minimap2 pipeline       
          mmap_out=os.path.join(args.out_p,"minimap2")
          setupWorkingDir(log,mmap_out)
          mins= Minimap2(args.m_path,args.cdna,args.faa,args.cov,args.idty,args.annot,args.threads,mmap_out)
          log(" Running minimap2")
          mins.runMinimap2()
          log(" Filter minimap2 gene models ")
          mins.filterMinimap2()
          log(" Adding UTRs to the gene models ")
          mins.addUTRs()
          log(" Cleaning and renaming gene models ")
          mins.formatMinimap2()
          log(" Transfer Functional annotation to the gene models")
          mins.addFunc()
          log(" Making final gene models")
          mins.formatFunc()
          log(" Validating final gene models")
          #gins.validateGFF()
          log(" Finished ")
          log(" Check " + mmap_out + "/minimap2_gene_models.gff" )
    else: 
          # Call gmap pipeline
          print(args.cdna)
          gmap_out=os.path.join(args.out_p,"gmap")
          setupWorkingDir(log,gmap_out)
          gins= Gmap(args.gm_path,args.cdna,args.faa,args.csp,args.cov,args.idty,args.dup,args.annot,args.threads,gmap_out)
          log(" Running gmap")          
          gins.runGmap()
          log(" Filter gmap gene models ")          
          gins.filterGmap()
          log(" Adding UTRs to the gene models ")
          gins.addUTRs()
          log(" Cleaning and renaming gene models ")
          gins.formatGmap()
          log(" Transfer Functional annotation to the gene models")
          gins.addFunc()
          log(" Making final gene models")
          gins.formatFunc()
          log(" Validating final gene models")
          #gins.validateGFF()
          log(" Finished ")
          log(" Check " + gmap_out + "/gmap_gene_models.gff" ) 
    log('goodbye')
