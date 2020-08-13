import pandas as pd
import sys, getopt
import argparse

parser = argparse.ArgumentParser(description='Transform raw read counts into RPKM.')
parser.add_argument('--biomart', help='Murine Biomart reference file')
parser.add_argument('--counts', help='Raw counts file')
parser.add_argument('--out_prefix', help='Prefix for RPKM file name (eg. GSE125197_rpkm')
parser.add_argument('--input_sep', help='column seperator for input count file options are "tab" or "csv"')
parser.add_argument('--ignore', help='Columns to ignore (zero-start) (eg. "1;3" will omit the 2nd and 4th column, use \'None\' to keep all the columns)')
args = parser.parse_args()

# try:
# 	opts, args = getopt.getopt(sys.argv,"",["--biomart ","--counts ","--out_prefix ","--input_sep ","--ignore "])
# except getopt.GetoptError:
# 	print 'test.py -i <inputfile> -o <outputfile>'
# 	sys.exit(2)

# Generate_RPKM.py --biomart /home/molly/Desktop/Cibersort/mm10_biomart_database.txt 
# --counts /home/molly/Desktop/Cibersort/GSE125197_lunde_et_al_2019_mouse_counts.txt
# --out_prefix /home/molly/Desktop/Cibersort/GSE125197_rpkm --input_sep tab --ignore "0;2"

# Generate_RPKM.py --biomart /home/molly/Desktop/Cibersort/mm10_biomart_database.txt 
# --counts /home/molly/Desktop/Cibersort/GSE135114_REH.genelevel.count.csv 
# --out_prefix /home/molly/Desktop/Cibersort/GSE135114_REH_rpkm --input_sep csv --ignore None

def intersection(lst1, lst2): 
    lst3 = [value for value in lst1 if value in lst2] 
    return lst3 

biomart=sys.argv[2]
# biomart="/home/molly/Desktop/Cibersort/mart_export.txt"
counts_file=sys.argv[4]
# counts_file="/home/molly/Desktop/Cibersort/GSE125197_lunde_et_al_2019_mouse_counts.txt"
out_file=sys.argv[6]
# out_file="/home/molly/Desktop/Cibersort/GSE125197_rpkm"
ct_sep=sys.argv[8] # options are csv or tab
# ct_sep="tab"
if sys.argv[10] == "None":
	print("No columns are being omitted from your input file.")
else:
	print("you are omitting columns "+sys.argv[10])
	ignore_cols = [int(x) for x in sys.argv[10].split(';')]
# ignore_cols="0;2"

## Read in the biomart resource
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
f['Gene_symbol_v3'] = f['Gene_symbol_v2'].str.replace("-","")

if ct_sep == "csv":
	counts = pd.read_csv(open(counts_file, "r"),header=0)
	try:
		ignore_cols
	except NameError:
		print("No columns being omitted...")
	else:
		to_r = counts.columns.values[ignore_cols] #establish which column names correspond to the column indexes to remove
		counts.drop(columns=to_r,axis=1,inplace=True)
	new_columns = counts.columns.values; new_columns[0] = 'ID'; counts.columns  = new_columns #replace the first column's name with "ID"
	counts['ID'] = counts['ID'].str.upper() # turn the ID into uppercase
	#remove NAN and duplicated IDs
	counts.dropna(axis=0,inplace=True)
	counts.drop_duplicates(subset='ID', keep="first",inplace=True)
	counts.shape	
elif ct_sep == "tab":
	counts = pd.read_csv(open(counts_file, "r"),header=0,sep='\t')
	try:
		ignore_cols
	except NameError:
		print("No columns being omitted...")
	else:
		to_r = counts.columns.values[ignore_cols] #establish which column names correspond to the column indexes to remove
		counts.drop(columns=to_r,axis=1,inplace=True)
	new_columns = counts.columns.values; new_columns[0] = 'ID'; counts.columns  = new_columns #replace the first column's name with "ID"
	counts['ID'] = counts['ID'].str.upper() # turn the ID into uppercase
	#remove NAN and duplicated IDs
	counts.dropna(axis=0,inplace=True)
	counts.drop_duplicates(subset='ID', keep="first",inplace=True)
	counts.shape

#Begin matching
print("Now identifying the type of gene identifier...")
genes_match = counts['ID'].str.startswith("ENSMUSG")
if sum(genes_match) > 0.5*counts.shape[0] :
	print("You are using Ensembl gene identifiers\nNow testing for gene versions...")
	genes_match = counts['ID'].str.endswith("\\.[0-9]*")
	if sum(genes_match) > 0.5*counts.shape[0] :
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
elif sum(counts['ID'].str.startswith("NM_")) > 0.5*counts.shape[0] :
	genes_match = counts['ID'].str.startswith("NM_")
	current_overlap = intersection(counts['ID'].values,f['Refseq'].values)
	outmess = "You seem to be using Refseq mRNA identifiers\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
	print(outmess)
	gene_type = "Refseq"
else:
	genes_match = counts['ID'].str.startswith("ENSMUST")
	if sum(genes_match) > 0.5*counts.shape[0] :
		print("You are using Ensembl transcript identifiers\nNow testing for transcript versions...")
		genes_match = counts['ID'].str.endswith("\\.[0-9]*")
		if sum(genes_match) > 0.5*counts.shape[0] :
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
		if len(current_overlap) > 0.5*counts.shape[0] :
			outmess = "You seem to be using gene symbols.\nOverlap with Biomart database: "+str(round((float(len(current_overlap))/float(counts.shape[0]))*100,2))+"%"
			print(outmess)
			gene_type="Gene_symbol_v2"
		else:
			print("Your gene identifiers do not appear to be Ensembl identifiers or gene symbols.\nExiting the program.")
			exit()

# Subset the counts to match the overlap with the gene type
current_overlap = intersection(counts['ID'].values,f[gene_type].values)
counts.set_index(counts['ID'],inplace=True) #Set index of counts to be the ID
counts.shape

counts = counts.loc[current_overlap]
print("Number of genes used for generating RPKM: "+str(counts.shape[0]))
print("First 5 genes: "+current_overlap[0]+" "+current_overlap[1]+" "+current_overlap[2]+" "+current_overlap[3]+" "+current_overlap[4])

#Subset the database to match the overlapping genes
f_sub = f[f[gene_type].isin(current_overlap)]
f_key = f_sub[[gene_type,'Gene_symbol_v2']].copy()
f_key.drop_duplicates(subset=gene_type,keep="first",inplace=True)
f_key.set_index(gene_type,inplace=True)
f_key.shape


f_sub = f_sub[[gene_type,'Length']]  #,'Gene_symbol_v2'
f_mean = f_sub.groupby(gene_type, as_index=False).mean()
f_mean.set_index(gene_type,inplace=True)
f_mean = f_mean.loc[current_overlap]

## Generate RPKM
## RPKM =   numReads / ( geneLength/1000 * totalNumReads/1,000,000 )
counts_rev = counts.drop("ID", axis=1).copy() # remove identifier column
colsum = counts_rev.copy().sum()/(10**6) # Sum count columns
genelength = f_mean['Length'].values/1000

rpkm = counts_rev.divide(genelength,axis=0).divide(colsum,axis=1)

#put the gene ids, symbols, and rpkm together
if gene_type == 'Gene_symbol_v2':
	rpkm_out = f_mean.merge(rpkm,left_index=True,right_index=True)
else:	
	rpkm_out = f_key.merge(f_mean,left_index=True,right_index=True).merge(rpkm,left_index=True,right_index=True)

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
# col_specs=[1]
# col_specs.extend(range(3,rpkm_out.shape[1]))
file_name = out_file+"_minimal.txt"
rpkm_out_min.to_csv(file_name,sep="\t",index=False)

