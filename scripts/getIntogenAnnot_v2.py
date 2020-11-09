import pandas as pd
import argparse, re, os
import numpy as np

parse = argparse.ArgumentParser(description=".")

parse.add_argument('-i', dest='input_file', required=True,
				   help='input, vep tsv file')
parse.add_argument('-c', dest='ctype', required=True,
				   help='cancertype of variants')
parse.add_argument('-o', dest='output_file', required=True,
				   help="Path to output.tsv file")

try:
	args = parse.parse_args()

except:
	exit('parsing of arguments failed!')

# parse args
infile_path = str(args.input_file)
outfile_path = str(args.output_file)
cType = str(args.ctype)
infile = pd.read_csv(infile_path, delimiter="\t")
gene_col = "SYMBOL"  # input
prot_col = "HGVSp"

# parse dbs
intogen_main = os.path.join(os.path.dirname(__file__))+"/data/"
gene_driver_path = intogen_main + "driver_genes.tsv"
var_af_path = intogen_main + "gene_var_pseudoAF_perCtype_v2.tsv"
var_tier_path = intogen_main + "variant_driver_tier.tsv"
gene_driver_file = pd.read_csv(gene_driver_path, delimiter="\t")
var_af_file = pd.read_csv(var_af_path, delimiter="\t")
var_driver_file = pd.read_csv(var_tier_path, delimiter="\t")

# translate tier to ranks
tier_ranks = {"1": 2, "2": 1, }  # and NA = 0, ranks for Variant driver tiers
na_value = np.nan  # new in V4 annot

cType_map_2019 = {
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

cType_map_2016 = {
	"COCA": ["COREAD"],
	"SK": ["CM"],
	"LUAD": ["LUAD","NSCLC", "LUSC", "SCLC"],
	"STAD": ["STAD"],
	"H": ["HC"],
	"BLCA": ["BLCA"],
	"BRCA": ["BRCA"],
	"BLTC": ["BLCA"],
	"PRAD": ["PRAD"],
	"ED": ["UCEC"],
	"RCCC" : ["RCCC"],
	"BRCADC": ["BRCA"],
	"B": ["GBM","LGG"],
	"READ": ["COREAD"],
	"PA": ["PA", "PAAD"],
	"LUSC": ["LUSC"],
	"HC": ["HC"],
	"TH": ["THCA"],
	"HNSC": ["HNSC"],
	"ESSC": ["ESCA"],
	#"BT": [],
	"SCLC": ["SCLC"],
	"SCC": ["CM"],
	"AML": ["AML"],
	"BCC": ["CM"],
	#"CANCER": [],
	"OVSE": ["OV"],
	#"CESC": ["CESC"],  # removed b/c no matching types in 2016 data
	"TCALL": ["ALL"],
	"L": ["LUAD", "NSCLC", "LUSC", "SCLC"],
	#"OCSC": [],
	#"MDS": ["MDPS"],  # no 2016 matching
	"BRCAL": ["BRCA"],
	"CLL": ["CLL"],
	"DLBCL": ["DLBC"],
	"HEMATO": ["MM", "AML", "CLL"],
	#"OC": [],
	"ESA": ["ESCA"],
	#"HGP": ["ANGS"],  # no 2016 match
	"ESCA": ["ESCA"],
	"LY": ["DLBC"],
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
		return 0

	hits = list(var_af_file[(var_af_file[gene_col] == gene) & (var_af_file["cType"] == cType)][var_type])  # search for hits, variant is missense or nonsense
	print(var_af_file)
	print(var_af_file[var_af_file["cType"] == cType])
	print(var_af_file[var_af_file[gene_col] == gene])
	print(hits)
	exit()

	if hits:  # entry found for gene + variant combination
		return list(var_af_file[var_af_file[gene_col] == gene][var_type])[0]
	else:  # no hits - af of 0
		return 0

def processSignleAFinto2(gene,var_type, cType_list):

	if type(var_type) is float:  # vars w/o protein annotation! FILTER BEFORE TRAINING!!!!!!
		return na_value

	hits = var_af_file[(var_af_file[gene_col] == gene)]  # search for hits, variant is missense or nonsense
	if hits.empty:  # no hits found - return NA
		return na_value
	else:  # entry found for gene, subset to relevant cTypes
		query = " or ".join(["cType == '%s'" % x for x in cType_list])  # construct query for rel. cTypes
		hits = hits.query(query)  # subset to relevant cTypes

		# if hits empty -> no cType match -> return NA, else return max found
		return na_value if hits[var_type].empty else max(hits[var_type])

def getSignleVarPseudoAF(table, cType):
	same_afs , diff_afs= [], []

	# get lists of same & diff cTypes to subset source table to
	same_cType_list = cType_map_2016[cType]  # list of same cType from intogen data
	cType_map_vals = [item for sublist in list(cType_map_2016.values()) for item in
					  sublist]  # s/o https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
	diff_cType_list = list(
		set(cType_map_vals) - set(same_cType_list))  # get all cTypes except ones described by input type

	# append afs to lists
	table.apply(lambda x: same_afs.append(processSignleAFinto2(x[gene_col], x["prot_var_class"], same_cType_list)), axis=1)
	#print("\n" + "-"*50 + "\n")
	table.apply(lambda x: diff_afs.append(processSignleAFinto2(x[gene_col], x["prot_var_class"], diff_cType_list)),
				axis=1)

	table["pseudo_som_af_same"] = same_afs
	table["pseudo_som_af_diff"] = diff_afs
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

def findBestGeneTierSamecType(gene, cType_list):

	gene_driver_subset = gene_driver_file[gene_driver_file[gene_col] == gene]  # subset to only relevant genes

	if len(gene_driver_subset) > 0:  # gene in subset

		if cType == "CANCER":  # if CANCER - take all results, pan cancer
			hits = list(gene_driver_subset["METHODS"])

		else:  # if not pan cancer, subset to relevant ones
			query = " or ".join(["cType == '%s'" % x for x in cType_list])  # construct query for rel. cTypes
			hits = list(gene_driver_subset.query(query)["METHODS"])  # query table to see relevant results

		return max(hits) if hits else na_value  # return max if hits found, else return 0

	else:  # gene not in set - return 0
		return na_value

def getGeneDriverTier(table, cType):

	same_cType_list = cType_map_2019[cType]  # list of same cType from intogen data
	same_cType_tiers = table.apply(lambda x: findBestGeneTierSamecType(x[gene_col], same_cType_list), axis=1)
	table["top_gene_driver_same"] = same_cType_tiers


	cType_map_vals = [item for sublist in list(cType_map_2019.values()) for item in sublist]  # s/o https://stackoverflow.com/questions/952914/how-to-make-a-flat-list-out-of-list-of-lists
	diff_cType_list = list(set(cType_map_vals) - set(same_cType_list))  # get all cTypes except ones described by input type

	diff_cType_tiers = table.apply(lambda x: findBestGeneTierSamecType(x[gene_col], diff_cType_list) ,axis=1)
	table["top_gene_driver_diff"] = diff_cType_tiers

	return table

infile = transcript3to1Letters(infile)  # translate 3 to 1 letter code, also: add 'prot_var_class' column
infile = getSignleVarPseudoAF(infile, cType)  # get max pseudoAF of gene & var type for same & diff
infile = getGeneDriverTier(infile, cType)  # get highest gene driver tier from intogen 2019 data
infile = getBestVarDriverTierTable(infile)  # get highest scoring variant rank

infile.drop([gene_col, prot_col], axis=1).to_csv(outfile_path, sep="\t", index=False)  # output

