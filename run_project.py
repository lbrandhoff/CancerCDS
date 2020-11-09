import os, pickle, argparse
import pandas as pd

parse = argparse.ArgumentParser(description=".")

parse.add_argument('-i', dest='input_file', required=True,
				   help='Path to input vcf')
parse.add_argument('-o', dest='output_direcory', required=True,
				   help="Path to output directory")
parse.add_argument('-p', dest='proba_cut', required=True,
				   help="minimum decision probability (0-1) cut-off from output variants (min: 0.5)")
parse.add_argument('-c', dest='ctype', required=False,
				   help="Cancer type acronym associated with input file (if none chosen will only run pan classification)")

try:
	args = parse.parse_args()

except:
	exit('parsing of arguments failed!')

# paths to files
vcf_in = str(args.input_file)
out_dir = str(args.output_direcory)
proba_cutoff = float(args.proba_cut)
cType = str(args.ctype).upper() if args.ctype else None

annot_type_script = "./annotate_type.sh"
annot_pan_part_script = "./annotate_pan_partial.sh"
annot_pan_script = "./annotate_pan.sh"

classifier_dir = "./classifiers/"
classifs = [x[:-4] for x in os.listdir(classifier_dir)]  # parse available classifiers

# tier cols & af cols get filled w/ 0
tier_cols_type = ["top_var_driver_tier", "top_gene_driver_diff", "top_gene_driver_same", "gene_druggable", "Hotspot"]
tier_cols_pan = ["top_var_driver_tier", "top_gene_driver_tier", "gene_druggable", "Hotspot"]
af_cols = ["gnomAD_AF", "AF"]  # 4 NA filling
exclude_cols = ['SYMBOL', 'HGVSp', 'HGVSg', "prot_var_class"]  # cols from annot not required for training
output_cols = ['HGVSg', 'SYMBOL', 'HGVSp', 'proba_1']  # cols for output

# order of features in RFC different from order after annotation, reorder to match
type_feat_order = ['AF', 'CADD_phred', 'Eigenraw', 'GERP++_RS', 'Hotspot', 'LRT_score', 'MetaLR_score', 'MetaSVM_score',
				   'PolyPhen', 'SIFT', 'VEST3_score', 'gene_druggable', 'gnomAD_AF', 'pseudo_som_af_diff',
				   'pseudo_som_af_same', 'top_gene_driver_diff', 'top_gene_driver_same', 'top_var_driver_tier']

pan_feat_order = ['PolyPhen', 'SIFT', 'LRT_score', 'GERP++_RS', 'gnomAD_AF', 'AF', 'CADD_phred', 'Eigenraw',
				  'MetaSVM_score', 'VEST3_score', 'MetaLR_score', 'Hotspot', 'gene_druggable',
				  'pseudo_som_af', 'top_gene_driver_tier', 'top_var_driver_tier']

gene_treat_file = pd.read_csv("./scripts/data/gene_treatment.tsv", sep="\t")

# ANNOTATION FUNCTIONS
def runBothAnnots(taip):

	# run type-specific annotation
	command = "sh %s %s" % (annot_type_script, " ".join([vcf_in, taip, out_dir]))
	os.system(command)

	# run pan annotation, takes temp files from type specific, only reruns intogen annot, then deletes temp files!
	command = "sh %s %s" % (annot_pan_part_script, out_dir)
	os.system(command)

def runPanOnlyAnnot():
	command = "sh %s %s" % (annot_pan_script, " ".join([vcf_in, out_dir]))
	os.system(command)

#CLASSIFICATION FUNCTIONS
def fillNAs(data_table, zero_cols):

	for col in zero_cols:  # replace empty AF values w/ 0 - b/c they are not present
		data_table[col] = data_table[col].fillna(0)

	data_table.fillna(data_table.median(), inplace=True)  # fill other cols w/ mean values
	return data_table

def parseClassifier(cType):
	parse_path = classifier_dir + cType + ".rfc"

	with open(parse_path, "rb") as infile:
		forest = pickle.load(infile)
		infile.close()

	return forest

def predictData(forest, pred_data, out_data):

	# generate prediction data, min proba_1 >0.5 obviously
	probas = pd.Series(forest.predict_proba(pred_data)[:, 1], name="proba_1")  # get probas of predicted 1's
	preds = pd.Series(forest.predict(pred_data), name="label")  # predict labels
	prediction_data = pd.concat([preds, probas], axis=1)  # combine

	# get output info
	sufficient_hits = pd.concat([out_data, prediction_data], axis=1)  # combine infos
	sufficient_hits = sufficient_hits[(sufficient_hits["label"] == 1) & (sufficient_hits["proba_1"] >= proba_cutoff)][output_cols]  # subset
	sufficient_hits = sufficient_hits.sort_values(by=["proba_1"], ascending=False)  # sort by probas

	# return w/ new index
	return sufficient_hits.reset_index(drop=True)

def joinTreatments(outp_data):
	return pd.merge(left=outp_data, right = gene_treat_file, on ="SYMBOL", how="left")

def classifyData(forest, pred_data, out_data):

	# generate prediction data, min proba_1 >0.5 obviously
	probas = pd.Series(forest.predict_proba(pred_data)[:, 1], name="proba_1")  # get probas of predicted 1's
	preds = pd.Series(forest.predict(pred_data), name="label")  # predict labels
	prediction_data = pd.concat([preds, probas], axis=1)  # combine

	# get output info
	sufficient_hits = pd.concat([out_data, prediction_data], axis=1)  # combine infos
	sufficient_hits = sufficient_hits[(sufficient_hits["label"] == 1) & (sufficient_hits["proba_1"] >= proba_cutoff)][output_cols]  # subset
	sufficient_hits = sufficient_hits.sort_values(by=["proba_1"], ascending=False)  # sort by probas

	# rename proba cols
	sufficient_hits.rename(columns={'proba_1': 'Probability'}, inplace=True)

	# exclude synonymous & non-coding (if present) ('%3D' in HGVSp == '=' in HGVSp == synonymous var)
	sufficient_hits = sufficient_hits[~(sufficient_hits["HGVSp"].str.contains("%3D", na=True))]

	# return w/ new index
	return joinTreatments(sufficient_hits.reset_index(drop=True))


if cType in classifs:
	print("RUNNING TYPE-SPECIFIC & PAN ANNOTATION W/ %s" % cType)
	print()

	# start annotation
	runBothAnnots(cType)  # produces outdir/type_annotation.tsv & outdir/pan_annotation.tsv

	# parse annot data & prep for classification
	type_annot = pd.read_csv(out_dir + "type_annotation.tsv", sep="\t",header=0)
	type_pred_data = fillNAs(type_annot.drop(exclude_cols, axis=1), af_cols + tier_cols_type)

	pan_annot = pd.read_csv(out_dir + "pan_annotation.tsv", sep="\t", header=0)
	pan_pred_data = fillNAs(pan_annot.drop(exclude_cols, axis=1), af_cols + tier_cols_pan)

	# load classifiers
	type_forest = parseClassifier(cType)
	pan_forest = parseClassifier("pan_all")

	# classify variants, reorder cols pre classification correctly to be safe, join treatments
	type_classif = classifyData(type_forest, type_pred_data[type_feat_order], type_annot[exclude_cols])
	pan_classif = classifyData(pan_forest, pan_pred_data[pan_feat_order], pan_annot[exclude_cols])

	# output data
	type_classif.to_csv(out_dir + "%s_classification.tsv" % cType, sep="\t")
	print("wrote:", out_dir + "%s_classification.tsv" % cType)

	pan_classif.to_csv(out_dir + "pan_classification.tsv", sep="\t")
	print("wrote:", out_dir + "pan_classification.tsv")


else:
	print("RUNNING PAN-ONLY ANNOTATION")
	print()
	# annotate
	runPanOnlyAnnot()  # produces outdir/pan_annotation.tsv

	# parse annot data, prep for classification
	pan_annot = pd.read_csv(out_dir + "pan_annotation.tsv", sep="\t", header=0)
	pan_pred_data = fillNAs(pan_annot.drop(exclude_cols, axis=1), af_cols + tier_cols_pan)

	# parse RFC, classify data
	pan_forest = parseClassifier("pan_all")
	pan_classif = classifyData(pan_forest, pan_pred_data[pan_feat_order], pan_annot[exclude_cols])

	# output
	pan_classif.to_csv(out_dir + "pan_classification.tsv", sep="\t")
	print("wrote:", out_dir + "pan_classification.tsv")


