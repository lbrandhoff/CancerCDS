import argparse, os, re
import numpy as np
import pandas as pd


parse = argparse.ArgumentParser(description=".")

parse.add_argument('-i', dest='input_file', required=True,
				   help='Path to input vcf')
parse.add_argument('-o', dest='output_name', required=True,
				   help="Path to output tsv ")
try:
	args = parse.parse_args()

except:
	exit('COSMIC_vep.py: parsing of arguments failed!')

# paths to files
input_path = str(args.input_file)
output_path = str(args.output_name)

# if output given -> check for .vcf format
if not output_path.endswith(".tsv"):
	exit("Please give output in .tsv format")

# open files
input_file = pd.read_csv(input_path, comment="#", delimiter="\t", header=None)

if len(input_file.columns) > 8:  # drop vcfsort cols if present
	input_file = input_file[[0,1,2,3,4,5,6,7]]

input_file.columns =["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]

# get header fields
header_file = open(input_path, "r")
header_regex = "([a-zA-Z_0-9+]*\|+[a-zA-Z_0-9]*)"
header_fields = []
exis_var = None  # bool 2 detect existing variation annot, only uebergangsweise - annot no longer in run vep script

for line in header_file:
	if line.startswith("##INFO=<ID=CSQ"):
		header = line
		header_fields = "".join(re.findall(header_regex,header)).split("|")
		header_fields = header_fields[1:]  # remove consequence, no required in output only for prioritisation

		if "Existing_variation" in header_fields:  # if present, remove - set bool for value removal
			exis_var = True
			header_fields = header_fields[:-1]
		break
header_file.close()

# none value
none_value = "NA"

consequence_rank = {
	"transcript_ablation" : 1,
	"splice_acceptor_variant" : 2,
	"splice_donor_variant" : 3,
	"stop_gained" : 4,
	"frameshift_variant" : 5,
	"stop_lost" : 6,
	"start_lost" : 7,
	"transcript_amplification" : 8,
	"inframe_insertion" : 9,
	"inframe_deletion" : 10,
	"missense_variant" : 11,
	"protein_altering_variant" : 12,
	"splice_region_variant" : 13,
	"incomplete_terminal_codon_variant" : 14,
	"start_retained_variant" : 15,
	"stop_retained_variant" : 16,
	"synonymous_variant" : 17,
	"coding_sequence_variant" : 18,
	"mature_miRNA_variant" : 19,
	"5_prime_UTR_variant" : 20,
	"3_prime_UTR_variant" : 21,
	"non_coding_transcript_exon_variant" : 22,
	"intron_variant" : 23,
	"NMD_transcript_variant" : 24,
	"non_coding_transcript_variant" : 25,
	"upstream_gene_variant" : 26,
	"downstream_gene_variant" : 27,
	"TFBS_ablation" : 28,
	"TFBS_amplification" : 29,
	"TF_binding_site_variant" : 30,
	"regulatory_region_ablation" : 31,
	"regulatory_region_amplification" : 32,
	"feature_elongation": 33,
	"regulatory_region_variant": 34,
	"feature_truncation": 35,
	"intergenic_variant": 36}

consequence_impact = {
	"transcript_ablation" : "1",
	"splice_acceptor_variant" : "1",
	"splice_donor_variant" : "1",
	"stop_gained" : "1",
	"frameshift_variant" : "1",
	"stop_lost" : "1",
	"start_lost" : "1",
	"transcript_amplification" :"1",  # last HIGH IMPACT
	"inframe_insertion" : "2",
	"inframe_deletion" : "2",
	"missense_variant" : "1",  # include, even though only medium, see ASCO guidelines
	"protein_altering_variant" : "2",
	"splice_region_variant" : "2",
	"incomplete_terminal_codon_variant" : "2",
	"start_retained_variant" : "2",
	"stop_retained_variant" : "2",
	"synonymous_variant" : "2",
	"coding_sequence_variant" : "2",
	"mature_miRNA_variant" : "2",
	"5_prime_UTR_variant" : "2",
	"3_prime_UTR_variant" : "2",
	"non_coding_transcript_exon_variant" : "2",
	"intron_variant" : "2",
	"NMD_transcript_variant" : "2",
	"non_coding_transcript_variant" : "2",
	"upstream_gene_variant" : "2",
	"downstream_gene_variant" : "2",
	"TFBS_ablation" : "2",
	"TFBS_amplification" : "2",
	"TF_binding_site_variant" : "2",
	"regulatory_region_ablation" : "2",
	"regulatory_region_amplification" : "2",
	"feature_elongation" : "2",
	"regulatory_region_variant" : "2",
	"feature_truncation" : "2",
	"intergenic_variant": "2"}

def processInfo(line, output_table):
	infos = line["INFO"][4:].split(",")  # get all csq annotations from line
	cons = [ x.split("|")[0] for x in infos]  # for each annot, get its consequence (first value)

	# get worst consequence for each of info fields in list of tuples [(id_infofield, conseq_rank, conseq_impact), ...]
	cons_numbers = [(id, consequence_rank[x], consequence_impact[x]) if "&" not in x else min([(id, consequence_rank[y], consequence_impact[y]) for  y in x.split("&")], key=lambda x:x[1]) for id, x in enumerate(cons)]

	# get consequence w/ highest rank, if multiple w/ same highest rank - choose first
	worst_cons = min(cons_numbers, key=lambda x:x[1])


	worst_info = infos[worst_cons[0]]  # get info from id via worst consequence tuple
	worst_info = [x if x else none_value for x in worst_info.split("|")]  # split & replace na's w/ variable
	worst_info = worst_info[1:]

	if exis_var:  # if exis var annot found remove it!
		worst_info = worst_info[:-1]

	output_table.loc[line.name] = worst_info

output_tsv = pd.DataFrame(columns=header_fields)  # create output table to fill

input_file.apply(lambda x: processInfo(x, output_tsv), axis=1)
output_tsv.to_csv(output_path, sep="\t", index=False)
