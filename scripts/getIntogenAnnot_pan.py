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

# parse args
infile_path = str(args.input_file)
outfile_path = str(args.output_file)
infile = pd.read_csv(infile_path, delimiter="\t")
gene_col = "SYMBOL"  # input
prot_col = "HGVSp"

# parse dbs
intogen_main = os.path.join(os.path.dirname(__file__))+"/data/"
gene_driver_path = intogen_main + "driver_genes.tsv"
var_af_path = intogen_main + "gene_var_pseudoAF_pan.tsv"
var_tier_path = intogen_main + "variant_driver_tier.tsv"
gene_driver_file = pd.read_csv(gene_driver_path, delimiter="\t")
var_af_file = pd.read_csv(var_af_path, delimiter="\t")
var_driver_file = pd.read_csv(var_tier_path, delimiter="\t")

tier_ranks = {"1": 2, "2": 1, }  # and NA = 0, ranks for Variant driver tiers
na_value = np.nan  # new in V4 annot

cType_map = {
	"COCA": ["COREAD"],
	"SK": ["CSCC", "SBCC", "CM"],
	"LUAD": ["LUAD"],
	"STAD": ["ST"],
	"H": ["HC"],
	"BLCA": ["BLCA"],
	"BRCA": ["BRCA"],
	"BLTC": ["BLCA"],
	"PRAD": ["PRAD"],
	"ED": ["UCEC"],
	"RCCC" : ["RCCC"],
	"BRCADC": ["BRCA"],
	"B": ["EPD","HGG"],
	"READ": ["COREAD"],
	"PA": ["PAIS", "PAAD"],
	"LUSC": ["LUSC"],
	"HC": ["HC"],
	"TH": ["THCA", "THYM"],
	"HNSC": ["HNSC"],
	"ESSC": ["ESCA"],
	#"BT": [],
	"SCLC": ["SCLC"],
	"SCC": ["SBCC", "CM"],
	"AML": ["AML"],
	"BCC": ["SBCC", "CM"],
	#"CANCER": [],
	"OVSE": ["OV"],
	"CESC": ["CESC"],
	"TCALL": ["ALL"],
	"L": ["LUAD", "LUSC", "LUNG_CANCER"],
	#"OCSC": [],
	"MDS": ["MDPS"],
	"BRCAL": ["BRCA"],
	"CLL": ["CLL"],
	"DLBCL": ["DLBCL"],
	"HEMATO": ["ALL", "AML", "CLL", "LY"],
	#"OC": [],
	"ESA": ["ESCA"],
	"HGP": ["ANGS"],
	"ESCA": ["ESCA"],
	"LY": ["LY"],
	"NSPH": ["NSPH"]
}

def transcript3to1Letters(table):
	prots = table[prot_col]
	amino_dict = {'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Asx': 'B', 'Cys': 'C', 'Glu': 'E', 'Gln': 'Q',
				  'Glx': 'Z', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F',
				  'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V', 'Ter': 'X', "fsX": "fs*"}

	missense_regex = "^[ARNDBCEQZGHILKMFPSTWYV]{1}[0-9]+[ARNDBCEQZGHILKMFPSTWYV]{1}"

	for three_letter, one_letter in amino_dict.items():
		# convert 3 to 1 letter code, if vars have ext entry -> cut off Lys234TerextTer41
		prots = [ x.replace(three_letter, one_letter).split("ext")[0] if type(x) is not float else x for x in prots]
		# TODO: improve classification -> N34Gfs*4 gets missense b/c of regex hit

	# classify converted vars into missoense and nonsense
	var_class = [x if type(x) is float else "Missense" if re.search(missense_regex,x) else "Nonsense" for x in prots]
	var_class = [x if type(x) is float else "Nonsense" if "fs" in x else "Nonsense" if x.endswith("X") else "Missense" for x in prots]

	table[prot_col] = prots
	table["prot_var_class"] = var_class

	return table

def getBothVarPseudoAF(table):
	return table.merge(var_af_file, how="left", on="SYMBOL")

def processSignleAFinto(gene,var_type):

	if type(var_type) is float:  # vars w/o protein annotation! FILTER BEFORE TRAINING!!!!!!
		return na_value

	hits = list(var_af_file[var_af_file[gene_col] == gene][var_type])  # search for hits, variant is missense or nonsense

	if hits:  # entry found for gene + variant combination
		return list(var_af_file[var_af_file[gene_col] == gene][var_type])[0]
	else:  # no hits - af of 0
		return na_value

def getSignleVarPseudoAF(table):
	collect_afs = []

	table.apply(lambda x: collect_afs.append(processSignleAFinto(x[gene_col], x["prot_var_class"])), axis=1)

	table["pseudo_som_af"] = collect_afs
	return table

def findBestVarTier(gene, var):

	# find protein and variant matches
	res = set(var_driver_file[(var_driver_file["gene"] == gene) & (var_driver_file["protein"] == var)]["driver_mut_prediction"])

	if res:  # if hits found - return best tier value
		res = [x[-1] for x in res]  # get only numbers of tiers
		return tier_ranks[min(res)]  # return translated rank
	else:  # if no hits found
		return na_value

def getBestVarDriverTierTable(table):

	collect_tiers = []
	table.apply(lambda x: collect_tiers.append(findBestVarTier(x[gene_col], x[prot_col])), axis=1)
	table["top_var_driver_tier"] = collect_tiers
	return table

def getGeneDriverTier(table):
	tese = [na_value if gene_driver_file[gene_driver_file[gene_col] == x].empty else  # if subset is emtpty -> return na_value
			max(gene_driver_file[gene_driver_file[gene_col] == x]["METHODS"])  # sub subset nonempty, return max
			for x in table[gene_col]]
	table["top_gene_driver_tier"] = tese

	return table

infile = transcript3to1Letters(infile)  # translate 3 to 1 letter code, also: add 'prot_var_class' column
infile = getSignleVarPseudoAF(infile)  # get max pseudoAF of gene & var type for same & diff
infile = getGeneDriverTier(infile)  # get highest gene driver tier from intogen 2019 data
infile = getBestVarDriverTierTable(infile)  # get highest scoring variant rank
infile.drop([gene_col, prot_col], axis=1).to_csv(outfile_path, sep="\t", index=False)  # output

