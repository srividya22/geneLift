#!/usr/bin/env python

"""
App to add explicit UTR exon records (as inferred from the gene, mRNA, exon and CDS features) to GFF3 data.

Example command:
add_utrs_to_gff.py input.gff3 > output_with_utrs.gff3
"""

############################################################################

__license__ = """\
                       PUBLIC DOMAIN NOTICE
          National Center for Biotechnology Information

This software/database is a "United States Government Work" under the
terms of the United States Copyright Act.  It was written as part of
the author's official duties as a United States Government employee and
thus cannot be copyrighted.  This software/database is freely available
to the public for use. The National Library of Medicine and the U.S.
Government have not placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy
and reliability of the software and data, the NLM and the U.S.
Government do not and cannot warrant the performance or results that
may be obtained by using this software or data. The NLM and the U.S.
Government disclaim all warranties, express or implied, including
warranties of performance, merchantability or fitness for any particular
purpose.

Please cite the author in any work or product based on this material.
"""

__copyright__ = "Public Domain"
__author__ = "David Managadze"
__version__ = "1.0.0"
__status__ = "Production"

############################################################################

# region imports
from sys import stdout, argv, exit
from collections import namedtuple
import gzip
import re
from textwrap import dedent
# from argparse import ArgumentParser
# endregion

# region constants and variables
GFF_FIELDS = ["seq_id", "source", "seq_type", "start", "end", "score", "strand", "phase", "attributes"]
GFFRecord = namedtuple("GFFRecord", GFF_FIELDS)
cur_gene = None
cur_mrna = None
# endregion


# region classes

class GFFGene(object):
    """ Create data structure for storing genes. """

    def __init__(self, gff_record, transcripts=None, cdss=None, utrs=None, misc_features=None):
        self.gff_record = gff_record
        self.name = parse_gff_attributes(self.gff_record.attributes)["Name"]
        self.id = parse_gff_attributes(self.gff_record.attributes)["ID"]
        self.transcripts = transcripts if isinstance(transcripts, list) else []
        self.cdss = cdss if isinstance(cdss, list) else []
        self.utrs = utrs if isinstance(utrs, list) else []
        self.misc_features = misc_features if isinstance(misc_features, list) else []
        self.processed = False


class GFFTranscript(object):
    """ Create data structure for storing transcripts. """

    def __init__(self, gff_record, exons=None):
        self.gff_record = gff_record
        self.id = parse_gff_attributes(self.gff_record.attributes)["ID"]
        self.exons = exons or []
        # self.transcript = transcript


class GFFCds(object):
    """ Create data structure for storing CDS. """

    def __init__(self, gff_record, exons=None):
        self.gff_record = gff_record
        self.exons = exons or []
        # self.cds = cds

# endregion: classes


# region functions

def parse_gff_attributes(attr_str):
    """
    Parse GFF3 attribute column and return a dict with all attributes.
    :param attr_str:
    :type attr_str:
    :return: dictionary of attributes
    :rtype: dict
    """
    retval = dict()     # return value
    if attr_str == ".":
        return {}
    for attr in attr_str.split(";"):
        if not attr: continue
        key, value = attr.split("=")
        # if we want decoding of the URL-encoded string we can use
        # retval[urllib.unquote(key)] = urllib.unquote(value)
        # but we do not need it
        retval[key] = value
    return retval


def gene_utrs(gene):
    """
    Process data of the gene, return UTR records.
    :param gene: gene to analyze and find UTRs
    :type gene: GFFGene
    :return: UTRs of this gene
    :rtype: list
    """
    def add_utr(utr_start, utr_end, utr_type, exon, utrs):
        """
        Create GFFRecords of UTRs from: utr_start, utr_end, utr_type, exon
        Add it to the list: utrs (so it appends to the external list variable, does not return anything)
        """
        utr_attrs = list()
        exon_attrs = parse_gff_attributes(exon.attributes)
        #print exon_attrs
        #id_str = 'ID=utr' + re.sub(r'\D', "", exon_attrs['ID'])
        id_str = 'ID=' + utr_type + ":" + exon_attrs['Parent']
        
        utr_attrs.append(id_str)
        for attr in ('Parent', 'transcript_id', 'Dbxref', 'partial', 'start_range', 'end_range'):
            if attr in exon_attrs:
                attr_str = '{}={}'.format(attr, exon_attrs[attr])
                utr_attrs.append(attr_str)
        utr_attrs_str = ';'.join(utr_attrs)
        utr = create_gff_record(
            seq_id=exon.seq_id,
            source=exon.source,
            seq_type=utr_type,
            start=utr_start,
            end=utr_end,
            score=".",
            strand=exon.strand,
            phase=".",
            attributes=utr_attrs_str
        )
        utrs.append(utr)

    if not (gene and gene.transcripts and gene.cdss):
        return []
    utrs = list()
    cdss_start = None
    cdss_end = None
    # get leftmost and rightmost CDS boundaries
    for cds in gene.cdss:
        if cdss_start is None or cds.gff_record.start < cdss_start:
            cdss_start = cds.gff_record.start
        if cdss_end is None or cds.gff_record.end > cdss_end:
            cdss_end = cds.gff_record.end
    # based on strand, set GFF record type
    left_utr_type = "five_prime_UTR"
    right_utr_type = "three_prime_UTR"
    if gene.transcripts[0].gff_record.strand == '-':
        left_utr_type = "three_prime_UTR"
        right_utr_type = "five_prime_UTR"

    # analyze transcripts to find their UTRs
    for transcript in gene.transcripts:
        # analyze exons
        for exon in transcript.exons:
            # check left terminus
            if exon.start < cdss_start:
                utr_start = exon.start
                utr_end = exon.end
                utr_type = left_utr_type
                if exon.end >= cdss_start:
                    utr_end = cdss_start - 1
                add_utr(utr_start, utr_end, utr_type, exon, utrs)
            # check right terminus
            if exon.end > cdss_end:
                utr_start = exon.start
                utr_end = exon.end
                utr_type = right_utr_type
                if exon.start <= cdss_end:
                    utr_start = cdss_end + 1
                add_utr(utr_start, utr_end, utr_type, exon, utrs)
    return utrs


def create_gff_record(seq_id, source, seq_type, start, end, score, strand, phase, attributes):
    """ Create GFF record from the provided variables. """
    gff_dict = {
        # if we want decoding of the URL-encoded strings we can use
        # "seq_id": None if seq_id == "." else urllib.unquote(seq_id),
        # but we do not need it
        "seq_id": None if seq_id == "." else seq_id,
        "source": None if source == "." else source,
        "seq_type": None if seq_type == "." else seq_type,
        "start": None if start == "." else int(start),
        "end": None if end == "." else int(end),
        "score": None if score == "." else float(score),
        "strand": None if strand == "." else strand,
        "phase": None if phase == "." else phase,
        "attributes": attributes
    }
    gff_rec = GFFRecord(**gff_dict)
    return gff_rec


def gff_record_to_str(record):
    """
    Convert GFFRecord namedtuple to GFF string and return final string.
    :param record: GFF record to convert to string
    :type record: GFFRecord
    :return: string representation of the GFFRecord
    :rtype: str
    """
    seq_id = record.seq_id or "."
    source = record.source or "."
    seq_type = record.seq_type or "."
    start = record.start or "."
    end = record.end or "."
    score = record.score or "."
    strand = record.strand or "."
    phase = record.phase or "."
    attributes = "."
    if record.attributes != "." and isinstance(record.attributes, str):
        attributes = record.attributes
    elif isinstance(record.attributes, dict):
        # just in case we decide that input attributes should be dictionary this will convert it to string
        attr_list = ['='.join([key, value]) for key, value in list(record.attributes.items())]
        attributes = ';'.join(attr_list)
    else:
        raise Exception("ERROR: provided attribute is neither string nor dictionary:" + record.attributes)
    retval = '\t'.join([str(r) for r in [seq_id, source, seq_type, start, end, score, strand, phase, attributes]])
    return retval


def reset_gene():
    """ Reset global variables that are used to keep track of the current gene data. """
    global cur_gene, cur_mrna
    cur_gene = None
    cur_mrna = None


def utrs_to_str(utr_list):
    """ Print UTRs to stdout in GFF3 format. """
    # if different sorting becomes necessary, I can use
    # for utr in sorted(utr_list, key=attrgetter('attributes', 'seq_type', 'start')):
    utrs_str = '\n'.join(map(gff_record_to_str, utr_list))
    return utrs_str


def create_gene_gff_record_from_line(line):
    """ Create and return GFF record from the line (string) only if line contains records that we need. """
    fields = line.strip().split("\t")
    if len(fields) != len(GFF_FIELDS):
        return None
    if fields[2] not in ("gene", "mRNA", "CDS", "exon", "C_gene_segment", "V_gene_segment", "D_gene_segment",
                         "J_gene_segment"):
        return None
    gff_record = create_gff_record(*fields)
    return gff_record


def analyze_gff_record(gff_record):
    """ Analyze GFF record. If this is a new gene then analyze and process previous gene. """
    global cur_gene, cur_mrna
    # check type of the record and act accordingly
    if gff_record.seq_type == "gene":
        # this is gene: process old gene or create new empty (first?) gene
        if cur_gene is not None:
            utrs = gene_utrs(cur_gene)
            cur_gene.processed = True
            reset_gene()
            # create new gene based on gff record
            cur_gene = GFFGene(gff_record=gff_record)
            return utrs
        else:
            reset_gene()
            cur_gene = GFFGene(gff_record=gff_record)
            return None
    elif gff_record.seq_type in ("mRNA", "C_gene_segment", "V_gene_segment", "D_gene_segment", "J_gene_segment"):
        # this is mRNA or equal type of record: create transcript and add it to current gene's transcripts
        cur_mrna = GFFTranscript(gff_record=gff_record)
        cur_gene.transcripts.append(cur_mrna)
    elif gff_record.seq_type == "exon":
        # this is exon: add it to current transcript's exons
        if cur_mrna is not None:
            cur_mrna.exons.append(gff_record)
    elif gff_record.seq_type == "CDS":
        # this is CDS: create GFFCds object and add it to current gene's CDS list
        cds = GFFCds(gff_record=gff_record)
        cur_gene.cdss.append(cds)
    else:
        raise Exception("ERROR: Something went wrong while processing GFF record:" + str(gff_record))


def main(fh_in, fh_out=stdout):
    """
    Main method that reads from input file handle, starts UTR analysis and prints GFF3 to output file handle.
    :param fh_in: input file handle
    :type fh_in: file
    :param fh_out: output file handle (writable)
    :type fh_out: file
    """
    gff_record = None
    # parse lines
    for line in fh_in:
        if line.startswith("#") and not line.startswith("###"):
            # this is COMMENT/PRAGMA, just print it
            fh_out.write(line)
            continue
        elif line.startswith("###"):
            # this is SEGMENT END, process current gene and print the '###' to denote the end of segment
            if cur_gene is not None:
                utrs = gene_utrs(cur_gene)
                if utrs:
                    fh_out.write(utrs_to_str(utrs) + '\n')
            reset_gene()
            fh_out.write(line)
        elif line.strip() != '':
            # this is data line, process it
            gff_record = create_gene_gff_record_from_line(line)     # treat only if record type is what we want
            if isinstance(gff_record, GFFRecord):
                # if gff_record is gene, get and print UTRs of previous gene
                utrs = analyze_gff_record(gff_record)
                if utrs:
                    fh_out.write(utrs_to_str(utrs) + '\n')
            fh_out.write(line)
        else:
            fh_out.write(line)
    # finally,
    # at the end of the file, if unprocessed gene record is left, process it
    if cur_gene is not None:
        if not cur_gene.processed:
            utrs = gene_utrs(cur_gene)
            cur_gene.processed = True
            if utrs:
                fh_out.write(utrs_to_str(utrs) + '\n')

# endregion: functions


if __name__ == "__main__":
    message = """
    Application to add explicit five_prime_UTR and three_prime_UTR features (as
    inferred from the gene, mRNA, exon and CDS features) to GFF3 data.
    The input file should be a GFF3 or gziped GFF3.
    The application was developed for NCBI's GFF3 files, and may not work with all GFF3 files.

    Usage:
    add_utrs_to_gff.py input.gff3[.gz] > output_with_utrs.gff3"
    """
    if len(argv) == 1 or len(argv) > 2 or argv[1] == '-h' or argv[1] == '-help':
        print(dedent(message))
        exit(1)
    fin = argv[1]   # input file name
    file_open = gzip.open if fin.endswith(".gz") else open
    fh_input = file_open(fin)   # input file handle
    main(fh_input)
    fh_input.close()
