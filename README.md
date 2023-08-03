# phylogenetics-sumstat-analysis

## Background

This project provides a comparison of the effectivness of different summarizing statistics used to map the spread of disease, most notably the direction of infection(who infected whom).


## Running the code
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

You must also have matplotlib installed. If you do not have it installed, install using pip by:
"""
pip install matplotlib
"""


From there run the program with:
"""
python3 sumstat_implement.py
"""

Note: It may give warnings about the Y axis being depricated. THIS DOES NOT AFFECT THE RESULTS OF THE PROGRAM. The code should still run and the files should still be created.

 The graphs produced by the program are automatically saved as pdf files under the newly created directory: sumstat_graphs. Within this new directory are two subdirectories called sumstat12graphs and sumstat_histograms.After entering the directory with:
"""
cd sumstat_graphs/sumstat12graphs/
"""

You can then open the pdf of your choice with:
"""
evince FILENAME.pdf
"""
It also creates a file in the main directory called correctsumstats.pdf which is an assessment of the accuracy for the first 3 summary statistics. It then shows that accuracy of a majority rules system where it counts how many trees return 2 out of the 3 sumstats as correct identifications.

## Building Paper and Slides
To construct the research paper and slide show from their latex files. First, make sure you have all the nesessary latex files. If you are running a bare-bones version of latex you likely will not be able to construct the paper and slides.
To get all of the necessary files you can do a
"""
sudo apt install texlive-full
"""
It is important to note that this installation is over 6GB so it may be worth your time to look into smaller installations if size is an issue.

To construct the paper you must first enter the paper directory using:
"""
cd loganmoore_paper
"""
from there you must run three commands to ensure that the bib file works with the latex file properly:
"""
pdflatex loganmoore_paper.tex
bibtex loganmoore_paper
pdflatex loganmoore_paper.tex
pdflatex loganmoore_paper.tex
"""
from there you can open the pdf with:
"""
evince loganmoore_paper.pdf
"""

To construct the slides you must enter the slide directory using:
"""
cd loganmoore_slides
"""
then create the pdf with:
"""
pdflatex loganmoore_slides.tex
"""
then open the pdf with:
"""
evince loganmoore_slides.pdf
"""