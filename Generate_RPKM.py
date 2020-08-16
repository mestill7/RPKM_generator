import pandas as pd
import sys, getopt
import argparse
import wget
import os

### Example uses

# Generate_RPKM.py --biomart /home/molly/Desktop/Cibersort/mm10_biomart_database.txt 
# --counts /home/molly/Desktop/Cibersort/GSE125197_lunde_et_al_2019_mouse_counts.txt
# --out_prefix /home/molly/Desktop/Cibersort/GSE125197_rpkm --input_sep tab --ignore "0;2" --org="mouse"

# Generate_RPKM.py --counts /home/molly/Desktop/Cibersort/GSE135114_REH.genelevel.count.csv 
# --out_prefix /home/molly/Desktop/Cibersort/GSE135114_REH_rpkm --input_sep csv --org="human"

parser = argparse.ArgumentParser(description='Transform raw read counts into RPKM.',
    epilog="Report issues or feature requests to Github (https://github.com/mestill7/RPKM_generator)")
parser.add_argument('--biomart=',type=str, help='(Optional) Murine Biomart reference file')
parser.add_argument('--counts=',type=str, help='Raw counts file')
parser.add_argument('--out_prefix=',type=str, help='Prefix for RPKM file name (eg. GSE125197_rpkm')
parser.add_argument('--input_sep=',type=str, help='column seperator for input count file options are "tab" or "csv"')
parser.add_argument('--ignore=',type=str, help='(Optional) Columns to ignore (zero-start), seperated by semicolon (eg. "1;3" will omit the 2nd and 4th column, "0" will omit the first column')
parser.add_argument('--lengthtype=',type=int, help='(Optional) Input options: 1, 2 or 3. \n1=mean of lengths of isoforms, \n2=median of lengths of isoforms.\n3=max of lengths of isoforms.\n Default is 1.')
parser.add_argument('--org=',type=str, help='(Optional) Input options: mouse or human. Default is mouse.')
args = parser.parse_args()

# print 'ARGV      :', sys.argv[1:]
try:
    opts, args = getopt.getopt(sys.argv[1:],"",["biomart=","counts=","out_prefix=","input_sep=","ignore=","lengthtype=","org="])
except getopt.GetoptError:
    print("python Generate_RPKM.py --counts GSE125197_lunde_et_al_2019_mouse_counts.txt --out_prefix GSE125197_rpkm --input_sep tab --ignore \"0;2\"")
    sys.exit(2)

## Helper functions
def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 
## End helper functions

for option_key, option_value in opts:
    if option_key in ('--biomart'):
        biomart = option_value
    elif option_key in ('--counts'):
        counts_file = option_value
    elif option_key in ('--out_prefix'):
        out_file = option_value
    elif option_key in ('--input_sep'):
        ct_sep = option_value
    elif option_key in ('--ignore'):
        ignore_cols = [int(x) for x in option_value.split(';')]
    elif option_key in ('--lengthtype'):
        length_type = option_value
    elif option_key in ('--org'):
        organism = option_value

try:
    organism
except NameError:
    organism="mouse"
else:
    if organism not in ("human","mouse"):
        print(("You have picked the unknown organism \""+organism+"\".\nOnly mouse and human data are currently supported."))
        exit()

try:
    length_type
except NameError:
    length_type = 1
else:
    length_type = int(float(length_type))
    if length_type in (1,2,3):
        print(("Chosen length type: "+dict(list(zip((1,2,3),("mean","median","max"))))[length_type]))
    else:
        print(("You chose length type "+str(length_type)+". Please select 1 (mean), 2 (median) or 3 (maximum). Now exiting..."))
        exit()

try:
    biomart
except NameError:
    print(("No biomart database specified.\nNow attempting to downloading the "+organism+" database from biomart..."))
    if organism=="mouse":
        url = 'http://www.ensembl.org/biomart/martservice?query=<?xml%20version="1.0"%20encoding="UTF-8"?><!DOCTYPE%20Query><Query%20virtualSchemaName="default"%20formatter="TSV"%20header="1"%20uniqueRows="0"%20count=""%20datasetConfigVersion="0.6"><Dataset%20name="mmusculus_gene_ensembl"%20interface="default"><Attribute%20name="ensembl_gene_id"/><Attribute%20name="ensembl_gene_id_version"/><Attribute%20name="ensembl_transcript_id"/><Attribute%20name="ensembl_transcript_id_version"/><Attribute%20name="external_gene_name"/><Attribute%20name="refseq_mrna"/><Attribute%20name="transcript_length"/></Dataset></Query>'
        if os.path.isfile("./biomart_download.txt"):
            os.remove("./biomart_download.txt")
            biomart_file = wget.download(url,"./biomart_download.txt")
        else:
            biomart_file = wget.download(url,"./biomart_download.txt")
        print(biomart_file)
        biomart = pd.read_csv(open(biomart_file, "r"),header=0,sep='\t')
        biomart.columns = ['Gene_stable_ID','Gene_stable_ID_version','Txt_stable_ID','Txt_stable_ID_version','Gene_symbol','Refseq', 'Length']
        print(("\n\nYour database: "+str(biomart.shape[0])+" rows, "+str(biomart.shape[1])+" columns"))
        f=biomart.copy()
        f['Gene_symbol_v2'] = f['Gene_symbol'].str.upper()
        # f['Gene_symbol_v3'] = f['Gene_symbol_v2'].str.replace("-","")
        print("Here are the first rows of your database:")
        print((f.head()))
        print("Here are the column headers of your database:")
        print((f.columns.values))
    elif organism=="human": 
        url = 'http://www.ensembl.org/biomart/martservice?query=<?xml%20version="1.0"%20encoding="UTF-8"?><!DOCTYPE%20Query><Query%20virtualSchemaName="default"%20formatter="TSV"%20header="1"%20uniqueRows="0"%20count=""%20datasetConfigVersion="0.6"><Dataset%20name="hsapiens_gene_ensembl"%20interface="default"><Attribute%20name="ensembl_gene_id"/><Attribute%20name="ensembl_gene_id_version"/><Attribute%20name="ensembl_transcript_id"/><Attribute%20name="ensembl_transcript_id_version"/><Attribute%20name="external_gene_name"/><Attribute%20name="refseq_mrna"/><Attribute%20name="transcript_length"/></Dataset></Query>'
        if os.path.isfile("./biomart_download.txt"):
            os.remove("./biomart_download.txt")
            biomart_file = wget.download(url,"./biomart_download.txt")
        else:
            biomart_file = wget.download(url,"./biomart_download.txt")
        biomart = pd.read_csv(open(biomart_file, "r"),header=0,sep='\t')
        biomart.columns = ['Gene_stable_ID','Gene_stable_ID_version','Txt_stable_ID','Txt_stable_ID_version','Gene_symbol','Refseq', 'Length']
        print(("\n\nYour database: "+str(biomart.shape[0])+" rows, "+str(biomart.shape[1])+" columns"))
        f=biomart.copy()
        f['Gene_symbol_v2'] = f['Gene_symbol'].str.upper()
        # f['Gene_symbol_v3'] = f['Gene_symbol_v2'].str.replace("-","")
        print("Here are the first rows of your database:")
        print((f.head()))
        print("Here are the column headers of your database:")
        print((f.columns.values))
else:
    print("Using the user-provided biomart database...")
    # ## Read in the biomart resource
    f = pd.read_csv(open(biomart, "r"),header=0,sep='\t')
    #Format the database
    f.rename({'Gene stable ID':'Gene_stable_ID', 
    'Gene stable ID version': 'Gene_stable_ID_version', 
    'Transcript stable ID': 'Txt_stable_ID', 
    'Transcript stable ID version': 'Txt_stable_ID_version', 
    'Gene name': 'Gene_symbol', 
    'Transcript length (including UTRs and CDS)': 'Length',
    'RefSeq mRNA ID': 'Refseq'},
      axis=1,inplace=True)
    # Create possible modifications of the gene symbol
    f['Gene_symbol_v2'] = f['Gene_symbol'].str.upper()
    # f['Gene_symbol_v3'] = f['Gene_symbol_v2'].str.replace("-","")
    print(("Your database:"+str(biomart.shape[0])+"rows and"+str(biomart.shape[1])+" columns"))
    print("Here are the first rows of your database:")
    print((f.head()))
    print("Here are the column headers of your database:")
    print((f.columns.values))

print("\n\nReading counts file...")
if ct_sep == "csv":
    counts = pd.read_csv(open(counts_file, "r"),header=0)
    if ignore_cols is NameError:
        print("No columns being omitted.")
    else:
        print(("you are omitting columns "+str(ignore_cols)))
    to_r = counts.columns.values[ignore_cols] #establish which column names correspond to the column indexes to remove
    counts.drop(columns=to_r,axis=1,inplace=True)
    new_columns = counts.columns.values; new_columns[0] = 'ID'; counts.columns  = new_columns #replace the first column's name with "ID"
    counts['ID'] = counts['ID'].str.upper() # turn the ID into uppercase
    #remove NAN and duplicated IDs
    counts.dropna(axis=0,inplace=True)
    counts.drop_duplicates(subset='ID', keep="first",inplace=True)
    print(("\n\nYour counts file: "+str(counts.shape[0])+" rows, "+str(counts.shape[1])+" columns"))
    print("Here are the first rows of your counts:")
    print((counts.head()))
    print("Here are the column headers of your counts:")
    print((counts.columns.values))
elif ct_sep == "tab":
    counts = pd.read_csv(open(counts_file, "r"),header=0,sep='\t')
    if ignore_cols is NameError:
        print("No columns being omitted.")
    else:
        print(("you are omitting columns "+str(ignore_cols)))
    to_r = counts.columns.values[ignore_cols] #establish which column names correspond to the column indexes to remove
    counts.drop(columns=to_r,axis=1,inplace=True)
    new_columns = counts.columns.values; new_columns[0] = 'ID'; counts.columns  = new_columns #replace the first column's name with "ID"
    counts['ID'] = counts['ID'].str.upper() # turn the ID into uppercase
    #remove NAN and duplicated IDs
    counts.dropna(axis=0,inplace=True)
    counts.drop_duplicates(subset='ID', keep="first",inplace=True)
    print(("Your counts file: "+str(counts.shape[0])+" rows, "+str(counts.shape[1])+" columns"))
    print("Here are the first rows of your counts:")
    print((counts.head()))
    print("Here are the column headers of your counts:")
    print((counts.columns.values))

#Begin matching
print("Now identifying the type of gene identifier...")
if organism == "human":
    genes_match = counts['ID'].str.startswith("ENSG")
else:
    genes_match = counts['ID'].str.startswith("ENSMUSG")
if sum(genes_match) > (0.5*counts.shape[0]):
    print("You are using Ensembl gene identifiers\nNow testing for gene versions...")
    genes_match = counts['ID'].str.endswith("\\.[0-9]*")
    if sum(genes_match) > (0.5*counts.shape[0]):
        print("You are using Ensembl gene versions.")
        current_overlap = intersection(counts['ID'].values,f['Gene_stable_ID_version'].values)
        outmess = "You seem to be using Ensembl gene versions.\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
        print(outmess)
        gene_type = "Gene_stable_ID_version"
    else:
        current_overlap = intersection(counts['ID'].values,f['Gene_stable_ID'].values)
        outmess = "You seem to be using Ensembl gene identifiers.\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
        print(outmess)
        gene_type = "Gene_stable_ID"

elif sum(counts['ID'].str.startswith("NM_")) > (0.5*counts.shape[0]):
    genes_match = counts['ID'].str.startswith("NM_")
    current_overlap = intersection(counts['ID'].values,f['Refseq'].values)
    outmess = "You seem to be using Refseq mRNA identifiers\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
    print(outmess)
    gene_type = "Refseq"

else:
    if organism == "human":
        genes_match = counts['ID'].str.startswith("ENST")
    else:
        genes_match = counts['ID'].str.startswith("ENSMUST")
    if sum(genes_match) > (0.5*counts.shape[0]) :
        print("You are using Ensembl transcript identifiers\nNow testing for transcript versions...")
        genes_match = counts['ID'].str.endswith("\\.[0-9]*")
        if sum(genes_match) > (0.5*counts.shape[0]):
            print("You are using Ensembl transcript versions.")
            current_overlap = intersection(counts['ID'].values,f['Txt_stable_ID_version'].values)
            outmess = "You seem to be using Ensembl transcript versions.\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
            print(outmess)
            gene_type = "Txt_stable_ID_version"
        else:
            current_overlap = intersection(counts['ID'].values,f['Txt_stable_ID'].values)
            outmess = "You seem to be using Ensembl transcript identifiers.\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
            print(outmess)
            gene_type = "Txt_stable_ID"
    else:
        print("You may be using Gene symbols, now testing for gene symbol overlap...")
        current_overlap = intersection(counts['ID'].values,f['Gene_symbol_v2'].values)
        if len(current_overlap) > 0.5*counts.shape[0]:
            outmess = "You seem to be using gene symbols.\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
            print(outmess)
            gene_type="Gene_symbol_v2"
        else:
            print("Your gene identifiers do not appear to be Ensembl identifiers or gene symbols.\nDid you select the correct organism?\nExiting the program.")
            exit()


# Subset the counts to match the overlap with the gene type
current_overlap = intersection(counts['ID'].values,f[gene_type].values)
counts.set_index(counts['ID'],inplace=True) #Set index of counts to be the ID
counts.shape

counts = counts.loc[current_overlap]
print(("Number of genes used for generating RPKM: "+str(counts.shape[0])))
print(("First 5 genes: "+current_overlap[0]+" "+current_overlap[1]+" "+current_overlap[2]+" "+current_overlap[3]+" "+current_overlap[4]))

#Subset the database to match the overlapping genes
f_sub = f[f[gene_type].isin(current_overlap)]
f_key = f_sub[[gene_type,'Gene_symbol_v2']].copy()
f_key.drop_duplicates(subset=gene_type,keep="first",inplace=True)
f_key.set_index(gene_type,inplace=True)
f_key.shape

f_sub = f_sub[[gene_type,'Length']]  #,'Gene_symbol_v2'
if length_type == 1:
    f_len = f_sub.groupby(gene_type, as_index=False).mean()
    f_len.set_index(gene_type,inplace=True)
    f_len = f_len.loc[current_overlap]
if length_type == 2:
    f_len = f_sub.groupby(gene_type, as_index=False).median()
    f_len.set_index(gene_type,inplace=True)
    f_len = f_len.loc[current_overlap]
if length_type == 3:
    f_len = f_sub.groupby(gene_type, as_index=False).max()
    f_len.set_index(gene_type,inplace=True)
    f_len = f_len.loc[current_overlap]

## Generate RPKM
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
counts_rev = counts.drop("ID", axis=1).copy() # remove identifier column
colsum = counts_rev.copy().sum()/(10**6) # Sum count columns
genelength = f_len['Length'].values/1000

rpkm = counts_rev.divide(genelength,axis=0).divide(colsum,axis=1)

#put the gene ids, symbols, and rpkm together
if gene_type == 'Gene_symbol_v2':
    rpkm_out = f_len.merge(rpkm,left_index=True,right_index=True)
else:
    rpkm_out = f_key.merge(f_len,left_index=True,right_index=True).merge(rpkm,left_index=True,right_index=True)

rpkm_out.reset_index(inplace=True)
new_columns = rpkm_out.columns.values; new_columns[0] = gene_type; rpkm_out.columns  = new_columns

print("Now saving the RPKM files...")
file_name = out_file+"_full.txt"
rpkm_out.to_csv(file_name,sep="\t",index=False)

#Drop certain columns to make a user-friendly simple RPKM
if gene_type == 'Gene_symbol_v2':
    rpkm_out_min = rpkm_out.drop('Length', axis=1).copy()
else:
    rpkm_out_min = rpkm_out.drop([gene_type,'Length'], axis=1).copy()
file_name = out_file+"_minimal.txt"
rpkm_out_min.to_csv(file_name,sep="\t",index=False)
