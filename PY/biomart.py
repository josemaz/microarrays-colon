from apybiomart import find_marts, find_datasets, find_attributes, query
from termcolor import colored, cprint
import os
import pandas as pd
from pathlib import Path
from datetime import date

# pd.set_option('display.max_rows', None)

logprint = lambda x: cprint(x, 'red', attrs=["bold"])

def getBiomart(fname):
	# marts = find_marts()
	# # ENSEMBL_MART_ENSEMBL ==  Ensembl Genes 101
	# print(marts) 
	# ds = find_datasets(mart="ENSEMBL_MART_ENSEMBL")
	# # print(ds)
	# qry = ds["Dataset_name"].str.contains('[Hh][Uu][Mm][Aa][Nn]')
	# print(ds[qry])
	# attrs = find_attributes(dataset="hsapiens_gene_ensembl")
	# print(attrs)
	if not os.path.isfile(fname):
		print("Downloading Biomart ...")
		attrs = ["ensembl_gene_id","chromosome_name","start_position",
			"end_position","strand","band","percentage_gene_gc_content",
			"gene_biotype","external_gene_name"]
		chrs = [str(i) for i in range(1,22)]
		chrs.append('X')
		chrs.append('Y')
		bm = query(attributes=attrs,
			filters={"chromosome_name": chrs},
			dataset="hsapiens_gene_ensembl")
		# print(bm["Chromosome/scaffold name"].value_counts())
		bm = bm[bm["Gene type"] == "protein_coding"]
		bm.columns = ['stableID', 'chromName', 'gStart','gEnd', 'strand', 
			'band', 'gcCont', 'gBiotype', 'gName']
		bm.to_csv(fname,sep="\t", index=False)
	else:
		print("Reading Biomart ...")
		bm = pd.read_csv(fname,sep="\t")
	print("Biomart genes: ", bm.shape[0])
	return bm


###########################################################3
if __name__ == "__main__":
    today = date.today()
    Path("Data").mkdir(parents=True, exist_ok=True)
    os.chdir("Data")
    fn = "biomart-" + today.strftime("%Y%m%d") + ".tsv"
    getBiomart(fn)


