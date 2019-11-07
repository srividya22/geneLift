#!/usr/bin/env python

from __future__ import print_function
import os, sys, re
from geneLift import *
from geneLift.utils import *
#from utils import *

scripts_path="/seq/schatz/sramakri/sources/geneLift/geneLift/scripts/"
class Minimap2:
    """ 
    Generate gene models from minimap2
    """
    def __init__(self, mm_path,cDNA,faa,cov,idty,annot,threads,out):
        """
        :param cDNA: input cDNA fasta file
        :param faa: input assembly fasta file
        :param cov: coverage threshold
        :param idty: identity threshold
        :param out: Output path
        """
        self.cDNA = cDNA
        self.faa = faa
        self.cov = str(cov)
        self.idty = str(idty) 
        self.annot = annot
        self.out = out
        self.num_threads = str(threads)
        self.paf = os.path.join(out,"aln_minimap2.paf")
        self.final_gff3 = os.path.join(out,"minimap2_gene_models.gff3")
        
    def runMinimap2(self):
        """
        Run minimap2
        """
        dname = "".join(os.path.basename(self.faa).split(".")[:-1])
        ptool = "minimap2"
        minimap2_cmd = ( 'minimap2 -t '+ self.num_threads + ' -x splice -c -k 14 ' + self.faa + ' ' + self.cDNA + ' > ' + self.paf )
        run(minimap2_cmd)

    def filterMinimap2(self):
        """
        Filter minimap2 alignments based on identity and coverage 
        """
        script_file = os.path.join(scripts_path,"filter_minimap2_gff3.sh")
        filt_cmd = (script_file + ' ' + scripts_path  +' ' + self.paf + ' ' + self.faa  + ' ' + str(self.cov) + ' ' + str(self.idty) + ' '+ self.out ) 
        run(filt_cmd)

    def addUTRs(self): 
        """
        Add UTRs to the final annotation file
        """
        script_file = os.path.join(scripts_path,"../add_utrs_to_gff.py")
        utr_cmd = ( script_file  + ' ' + self.out + '/aln_minimap2_cds_sorted_clean.gff' + ' > ' + self.out + '/aln_minimap2_cds_UTRS_sorted.gff' )
        run(utr_cmd)
    
    def formatMinimap2(self):
        """
        Formatminimap2 alignments based on identity and coverage
        """
        script_file = os.path.join(scripts_path,"format_minimap2_gff3.sh")
        format_cmd = (script_file + ' ' + scripts_path  + ' ' + self.out + '/aln_minimap2_cds_UTRS_sorted.gff' + ' ' + self.annot  + ' ' + self.out )
        run(format_cmd) 
 
    def addFunc(self):
        """
        Add functional annotation to the gff file
        """ 
        script_file = os.path.join(scripts_path,"update_func_gff")
        func_cmd = ( script_file + ' ' + self.out + '/minimap2_mapped_func.txt' + ' ' +  self.out + '/aln_minimap2_cds_UTRS_sorted_modified_ids.gff' + ' ' + self.out + '/minimap2_mapped_ids.txt' + ' > ' + self.out + '/minimap2_func.gff' )
        run(func_cmd)

    def formatFunc(self):
        """
        Format functional annotation to the gff file
        """
        script_file = os.path.join(scripts_path,"format_func_minimap2_gff3.sh")
        format_func_cmd = (script_file + ' ' + scripts_path  + ' ' + self.out + '/minimap2_func.gff' + ' ' + self.cov + ' ' +  self.out)
        run(format_func_cmd)
      
    #def validateGFF(self):
    #    """
    #    Validate functional annotation to the gff file
    #    """
    #    script_file = os.path.join(scripts_path,"validate_and_correct.sh")
    #    validate_cmd = (script_file + ' ' + self.out + '/minimap2_func_final.gff' + ' ' + self.faa + ' ' + self.out)
    #    run(validate_cmd)
