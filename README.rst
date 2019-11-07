===============================
geneLift
===============================


.. image:: https://img.shields.io/travis/srividya22/geneLift.svg
        :target: https://travis-ci.org/srividya22/geneLift
.. image:: https://circleci.com/gh/srividya22/geneLift.svg?style=svg
    :target: https://circleci.com/gh/srividya22/geneLift
.. image:: https://codecov.io/gh/srividya22/geneLift/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/srividya22/geneLift
   

Gene model transfer from closely related reference genomes using cDNA alignments

usage: geneLift.py [-h] -cDNA CDNA -g FAA -func ANNOT [-o OUT_P] [-x]
                   [-aligner ALIGNER] [--report-duplications] [-mm M_PATH]
                   [-gm GM_PATH] [-c 90] [-i 95] [-t 1]

Script to liftover gene models from reference cDNA alignments

optional arguments:
  -h, --help            show this help message and exit
  -cDNA CDNA            fasta file containing transcript sequences
  -g FAA                Assembly fasta file to lift over gene models
  -func ANNOT           Semicolon (;) seperated text file with
                        transcript_id,gene_id,functional annotation for the
                        cDNA file
  -o OUT_P              output directory name
  -x                    enable cross species cDNA alignment ; only supported
                        with aligner gmap
  -aligner ALIGNER      Aligner used for cDNA alignments options:
                        [gmap,minimap2]
  --report-duplications
                        Report gene duplications
  -mm M_PATH            path to minimap2 executable
  -gm GM_PATH           path to gmap executable
  -c 90                 Coverage threshold for transcript alignments
  -i 95                 Identity threshold for transcript alignments
  -t 1                  Number of threads
