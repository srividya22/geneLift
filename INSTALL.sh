# Running this script should be enough to install geneLift,
# and set up the environment appropriately.  After this is executed,
# you should be able to get back to this environment with the command
#
#   source activate geneLift
#
# which sets up the conda environment, but also sources the geneLift
# environment scripts.  Note, however, that after the geneLift scripts are
# run, they cannot be undone (unlike conda's), so the environment will
# be affected until you start a new shell.

set -e  # Exit immediately if a command exits with a non-zero status.
set -x  # Print commands and their arguments as they are executed.

# Set up the conda environment, and update everything
conda create -y --name geneLift python=3.6
source activate geneLift
conda install -c bioconda bedtools gmap minimap2 genometools
conda install -c bioconda gffutils
conda update -y --all
pip install --upgrade pip
pip install intervaltree_bio
pip install gff3
chmod +x *.py geneLift/*.py 
chmod +x geneLift/geneLift/scripts/* 
