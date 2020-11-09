import pandas as pd
import numpy as np
import argparse, re, os


parse = argparse.ArgumentParser(description=".")

parse.add_argument('-i', dest='input_file', required=True,
				   help='input, vep tsv file')
parse.add_argument('-o', dest='output_file', required=True,
				   help="Path to output.tsv file")

try:
	args = parse.parse_args()

except:
	exit('parsing of arguments failed!')
na_value = np.nan
# parse args
infile_path = str(args.input_file)
outfile_path = str(args.output_file)
infile = pd.read_csv(infile_path, delimiter="\t")
gene_col = "SYMBOL"  # input
infile = pd.DataFrame(infile[gene_col])  # only consider symbol column, as DF for merging


gene_list_path = os.path.join(os.path.dirname(__file__))+"/data/drug_gene_interactions.tsv"
gene_list_file = pd.read_csv(gene_list_path, delimiter="\t")

res = infile.merge(gene_list_file, how="left", indicator=True)  # merge tables, use indicator to see wheter gene is present in both tables
druggable = res.apply(lambda x: 1 if x["_merge"] == "both" else na_value, axis=1)  # add info wheter gene is present in data
druggable.name = "gene_druggable"

druggable.to_csv(outfile_path, sep="\t", index=False, header=True)
