# RPKM_generator
A python script for generating RPKM (Reads Per Kilobase per Million) from raw count files. 

Usage:
python Generate_RPKM.py --help

  usage: Generate_RPKM.py [-h] [--biomart= BIOMART=] [--counts= COUNTS=]
                        [--out_prefix= OUT_PREFIX=] [--input_sep= INPUT_SEP=]
                        [--ignore= IGNORE=] [--lengthtype= LENGTHTYPE=]
                        [--org= ORG=]


Example Usage:
python Generate_RPKM.py --counts /home/molly/Desktop/Cibersort/GSE125197_lunde_et_al_2019_mouse_counts.txt 
--out_prefix /home/molly/Desktop/Cibersort/GSE125197_rpkm --input_sep tab --ignore "0;2" 

Uses the following python modules: pandas, sys, getopt, argparse, wget, os

Note:
Developed on Python2, but it should also work on Python3.
This script and the resulting RPKM files refer to the normalized output as RPKM. However, if your input represents fragments, rather than reads, than your output is actually FPKM, and you may choose to modify the filenames of your results accordingly.

Use of the "--lengthtype==3", which uses the maximum length of all available isoforms, may generate non-numeric RPKMS (e.g. "inf" instead of a number). It is recommended to use either "--lengthtype==1" (default) or "--lengthtype==2" options.


User-provided biomart files are permitted with the "--biomart=" flag, however, user-defined files are not fully checked for appropriate structure and datatypes. Therefore, when internet access is available, it is recommended to allow the script to automatically download the biomart database from Ensembl BioMart (https://m.ensembl.org/info/data/biomart/index.html)


Please report any issues or features you think would be useful.
