# ICR-2023-Phylogenetics

## Background

This project provides a comparison of the effectivness of different summarizing statistics used to map the spread of disease, most notably the direction of infection(who infected whom).


##Running the code
You must first make sure that the 'simulated_trees' file is in the same directory as the sumstat_implement.py file
To run the code you must first download and activate ete3(example using conda):


To install Miniconda(skip if already installed):
"""
wget http://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/anaconda_ete/
export PATH=~/anaconda_ete/bin:$PATH;
"""

To download ete3 using conda:
(further instructions for alternate installations at http://etetoolkit.org/download/)
"""
conda create -n ete3 python=3
conda activate ete3
conda install -c etetoolkit ete3 ete_toolchain
ete3 build check
"""


Then you must remember to activate your conda enviorment to use ete3:
"""
conda activate ete3
"""


From there run the program with:
"""
python3 sumstat_implement.py
"""


 The graphs produced by the program are automatically saved as pdf files under the newly created directory: sumstat_graphs.