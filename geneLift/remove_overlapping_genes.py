#!/usr/bin/env python
# 
# Author : Srividya Ramakrishnan
# Affliation : Johns Hopkins University

from __future__ import print_function
import csv
import sys, os, time
import collections
import argparse
#from bx.intervals.intersection import IntervalTree

from intervaltree_bio import GenomeIntervalTree

def getmRNAlengths(gff3_file_path):
    with open(gff3_file_path, 'r') as annotations_file:
        lmRNA = dict()
        reader = csv.reader(annotations_file, delimiter='\t')
        for row in reader:
            if len(row) == 9 and row[2] == "exon" and not row[0].startswith('##') :
                 attrb = dict(item.split("=") for item in row[8].split(";") if item)
                 mid = attrb['Parent']
                 e_len = int(row[4]) - int(row[3]) + 1
                 if mid in lmRNA.keys():
                    lmRNA[mid] += e_len
                 else:
                    lmRNA[mid] = e_len
    return lmRNA
       
def create_gene_tree(bed_file_path):

    # dictionary mapping chromosome names to interval trees
    models = dict()
    #dmRNA = getmRNAlengths(gff3_file_path)
    tree = GenomeIntervalTree()
    # parse the annotations file (GFF3) and build the interval trees
    with open(bed_file_path, 'r') as annotations_file:
        reader = csv.reader(annotations_file, delimiter='\t')
        for row in reader:
            if len(row) == 9 and not row[0].startswith('##'):
                seqid = row[0]
                start = int(row[1])
                end = int(row[2])
                strand = row[3]
                m_id = row[4]
                g_id = row[5]
                cov = float(row[6])
                idty = float(row[7])
                matches = int(row[8])
                #tree = None
                if tree[seqid].overlaps(start, end):
                      continue
                else:
                      models[m_id] = 1
                      models[g_id] = 1
                      tree[seqid].addi(start,end, data = ({"ID" : m_id , "Parent" : g_id }))
    return models

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Script to remove overlapping genes.' ) 
    parser.add_argument("bed_file",metavar="<bed file>", type=str, help="bed file to filter overlapping gene models")
    parser.add_argument("cov",metavar="<int>",default="90", type=int, help="Percentage of overlap allowed in the final gff")
    parser.add_argument("output",metavar="<output gff file>",type=str, help="Output gff file")

    # Get the command line arguments
    args = parser.parse_args()
    selected_models=create_gene_tree(args.bed_file)
    with open(args.output, "w") as f:
      for key in selected_models:
        print(key, file=f)
    
