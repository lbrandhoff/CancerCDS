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

#hotspot_path = "/mnt/users/ahbranl1/hotspots/prepped_hotspots.txt"
hotspot_path = os.path.join(os.path.dirname(__file__))+"/data/hotspots.txt"
hotspot_file = pd.read_csv(hotspot_path, delimiter="\t")

tolerance = 25
na_value = np.nan

def getTableHotspotPresence(table):
	hotspot_col = []
	table.apply(lambda x: getEntryHotSpotPresenceWithTolerance(x[gene_col], x[prot_col], hotspot_col, tolerance), axis=1)
	table["Hotspot"] = hotspot_col

	return table

def getEntryHotSpotPresence(gene, var, agg_list):
	position_regex = "[A-Z]{1}[0-9]+"

	if type(var) is not str:
		agg_list.append(0)
	else:
		var = re.search(position_regex, var)[0]
		res = hotspot_file[(hotspot_file[gene_col] == gene) & (hotspot_file["Mutation"] == var)]
		agg_list.append(1) if len(res) > 0 else agg_list.append(0)

def getEntryHotSpotPresenceWithTolerance(gene, var, agg_list, tolerance):

	# subset data table to relevant entries
	search_hotspots = hotspot_file[hotspot_file["SYMBOL"] == gene]  # subset table to only relevant gene
	search_hotspots["Mutation"] = [x[1:] for x in search_hotspots["Mutation"]]  # remove starting AS for range search

	if type(var) is not str:  # na protein annotation
		agg_list.append(na_value)
	elif len(search_hotspots) == 0:  # gene not in list
		agg_list.append(na_value)


	else:  # gene in list & viable variant
		pos = int(re.findall("[0-9]+", var)[0])  # get pos of query
		search_pos = range(pos-tolerance, pos+tolerance+1)  # get list of positions to search via tolerance
		query = " or ".join(["Mutation == '%s'" % x for x in search_pos])  # construct query for multiple positions
		hits = search_hotspots.query(query)["Mutation"].to_list()  # get position hits, vars further than distance not captured - get na anyway

		if not hits:  # if no hits found - end it
			agg_list.append(na_value)

		else:  # if hits found - evaluate best value
			hits = [abs(int(x)-pos) for x in hits]  # get absolute distance to hotspot from hits
			best_hit = min(hits)  # get best scoring hit
			best_hit = 1-(best_hit/(tolerance+1)) if best_hit > 0 else 1  # value transformation - the closer the better
			agg_list.append(best_hit)

infile = getTableHotspotPresence(infile)
infile.drop([gene_col, prot_col], axis=1).to_csv(outfile_path, sep="\t", index=False)

