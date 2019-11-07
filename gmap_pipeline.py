#!/usr/bin/env python

from __future__ import print_function
import os, sys, re
from geneLift import *
from geneLift.utils import *

scripts_path="/seq/schatz/sramakri/sources/geneLift/geneLift/scripts/"
class Gmap:
    """ 
    Generate gene models from gmap
    """
    def __init__(self, gm_path,cDNA,faa,csp,cov,idty,dup,annot,threads,out):
        """
        :param cDNA: input cDNA fasta file
        :param faa: input assembly fasta file
        :param cov: coverage threshold
        :param idty: identity threshold
        :param out: Output path
        """
        self.cDNA = cDNA
        self.faa = faa
        self.csp = csp
        self.cov = str(cov/100)
        self.idty = str(idty/100) 
        self.dup = dup
        self.annot = annot
        self.out = out
        self.num_threads = str(threads)
        self.gff3 = os.path.join(out,"aln_gmap.gff3")
        self.final_gff3 = os.path.join(out,"gmap_gene_models.gff3")
        
    def runGmap(self):
        """
        Run gmap
        """
        dname = "".join(os.path.basename(self.faa).split(".")[:-1])
        gmap_build_cmd = ('gmap_build -d ' + dname  + ' ' + self.faa )  
        run(gmap_build_cmd)
        ptool = "gmap"
        if self.csp:
           gmap_cmd = ( 'gmap -d '+ dname + ' --cross-species -B 5 -A -O -n 5 --no-chimeras --min-identity ' + self.idty + ' --min-trimmed-coverage ' + self.cov  + ' -f gff3_gene --gff3-swap-phase 1 -t ' + self.num_threads + ' ' + self.cDNA + ' > ' + self.gff3 )
        else:
           gmap_cmd = ( 'gmap -d '+ dname + ' -B 5 -A -O -n 5 --no-chimeras --min-identity ' + self.idty + ' --min-trimmed-coverage ' + self.cov  +' -f gff3_gene --gff3-swap-phase 1 -t ' + self.num_threads + ' ' + self.cDNA + ' > ' + self.gff3 )
        run(gmap_cmd)

    def filterGmap(self):
        """
        Filter gmap alignments based on identity and coverage 
        """
        script_file = os.path.join(scripts_path,"filter_gmap_gff3.sh")
        filt_cmd = (script_file + ' ' + scripts_path  +' ' + self.gff3 + ' ' + str(self.cov) + ' ' + str(self.idty) + ' '+ self.out ) 
        run(filt_cmd)

    def addUTRs(self): 
        """
        Add UTRs to the final annotation file
        """
        script_file = os.path.join(scripts_path,"../add_utrs_to_gff.py")
        utr_cmd = ( script_file  + ' ' + self.out + '/gmap_filtered_nonoverlapping.gff' + ' > ' + self.out + '/gmap_all.gff' )
        run(utr_cmd)
    
    def formatGmap(self):
        """
        Formatgmap alignments based on identity and coverage
        """
        script_file = os.path.join(scripts_path,"format_gmap_gff3.sh")
        format_cmd = (script_file + ' ' + scripts_path  + ' ' + self.out + '/gmap_all.gff' + ' ' + self.annot  + ' ' + self.out )
        run(format_cmd) 
 
    def addFunc(self):
        """
        Add functional annotation to the gff file
        """ 
        script_file = os.path.join(scripts_path,"update_func_gff")
        func_cmd = ( script_file + ' ' + self.out + '/gmap_mapped_func.txt' + ' ' +  self.out + '/gmap_format.gff' + ' ' + self.out + '/gmap_mapped_ids.txt' + ' > ' + self.out + '/gmap_func.gff' )
        run(func_cmd)

    def formatFunc(self):
        """
        Format functional annotation to the gff file
        """
        script_file = os.path.join(scripts_path,"format_func_gff3.sh")
        format_func_cmd = (script_file + ' ' + scripts_path  + ' ' + self.out + '/gmap_func.gff' + ' ' +  self.out)
        run(format_func_cmd)
      
    #def validateGFF(self):
    #    """
    #    Validate functional annotation to the gff file
    #    """
    #    script_file = os.path.join(scripts_path,"validate_and_correct.sh")
    #    validate_cmd = (script_file + ' ' + self.out + '/gmap_func_final.gff' + ' ' + self.faa + ' ' + self.out)
    #    run(validate_cmd)
