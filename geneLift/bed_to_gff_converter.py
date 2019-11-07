#!/usr/bin/env python
# This code exists in 2 places: ~/datatypes/converters and ~/tools/filters
import sys
import re

assert sys.version_info[:2] >= ( 2, 4 )

def __main__():
    input_name = sys.argv[1]
    output_name = sys.argv[2]
    skipped_lines = 0
    first_skipped_line = 0
    out = open( output_name, 'w' )
    out.write( "##gff-version 3\n\n" )
    #out.write( "##bed_to_gff_converter.py\n\n" )
    i = 0
    with open(input_name) as f:
     for i, line in enumerate( f ):
        complete_bed = False
        line = line.rstrip( '\r\n' )
        if line and not line.startswith( '#' ) and not line.startswith( 'track' ) and not line.startswith( 'browser' ):
            try:
                elems = line.split( '\t' )
                if len( elems ) == 12:
                    complete_bed = True
                chrom = elems[0]
                if complete_bed:
                    feature = "mRNA"
                    exonl = elems[10].split(",")[:-1]
                    mRNAl = sum(int(x) for x in exonl)
                else:
                    try:
                        feature = elems[3]
                    except:
                        feature = 'feature%d' % ( i + 1 )
                start = int( elems[1] ) + 1
                end = int( elems[2] )
                try:
                    score = elems[4]
                except:
                    score = '0'
                try:
                    strand = elems[5]
                except:
                    strand = '+'
                try:
                    group = elems[3]
                    regexp = re.compile(r'-c[0-9]*')
                    if regexp.search(group):
                        gene_id = ".".join(group.split(".")[:-1])
                        gene_id="%s-%s" % (gene_id,group.split("-")[-1])
                        #print(gene_id)
                    else:
                        gene_id = ".".join(group.split(".")[:-1])
                        #print(gene_id)
                except:
                    group = 'group%d' % ( i + 1 )
                if complete_bed:
                    out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\tID=%s;Name=%s;\n' % ( chrom,"gene", start, end, score, strand, gene_id, gene_id ) )
                    out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\tID=%s;Parent=%s;Name=%s;length=%s;\n' % ( chrom, feature, start, end, score, strand,group,gene_id,group,mRNAl ) )
                else:
                    print(feature)
                    out.write( '%s\tbed2gff\t%s\t%d\t%d\t%s\t%s\t.\t%s;\n' % ( chrom, feature, start, end, score, strand, group  ) )
                if complete_bed:
                    # We have all the info necessary to annotate exons for genes and mRNAs
                    block_count = int( elems[9] )
                    block_sizes = elems[10].split( ',' )
                    block_starts = elems[11].split( ',' )
                    for j in range( block_count ):
                        exon_start = int( start ) + int( block_starts[j] )
                        exon_end = exon_start + int( block_sizes[j] ) - 1
                        out.write( '%s\tbed2gff\texon\t%d\t%d\t%s\t%s\t.\tID=exon:%s;Parent=%s;\n' % ( chrom, exon_start, exon_end, score, strand, group,group ) )
            except:
                skipped_lines += 1
                if not first_skipped_line:
                    first_skipped_line = i + 1
        else:
            skipped_lines += 1
            if not first_skipped_line:
                first_skipped_line = i + 1
    out.close()
    info_msg = "%i lines converted to GFF version 2.  " % ( i + 1 - skipped_lines )
    if skipped_lines > 0:
        info_msg += "Skipped %d blank/comment/invalid lines starting with line #%d." %( skipped_lines, first_skipped_line )
    print(info_msg)

if __name__ == "__main__": __main__()
