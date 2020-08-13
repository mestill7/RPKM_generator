# RPKM_generator
A python script for generating RPKM files. Built on python2, should also work on python3.

Usage:
python Generate_RPKM.py --help

Example Usage:
python Generate_RPKM.py --biomart /home/molly/Desktop/Cibersort/mm10_biomart_database.txt  
--counts /home/molly/Desktop/Cibersort/GSE125197_lunde_et_al_2019_mouse_counts.txt 
--out_prefix /home/molly/Desktop/Cibersort/GSE125197_rpkm --input_sep tab --ignore "0;2"

Note:
This resource is still a work in progress. Please maintain the noted order of the input arguments, as the argument parser code is not yet fully optimized. 
It was designed for mouse. Future improvements will include an option to process read counts from human datasets. 
Please report any issues or features you think would be useful.
