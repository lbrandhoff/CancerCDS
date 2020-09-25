import pandas as pd
import argparse, re
from pyfaidx import Fasta


parse = argparse.ArgumentParser(description=".")

parse.add_argument('-i', dest='input_file', required=True,
                   help='input, vep tsv file')
parse.add_argument('-o', dest='output_file', required=True,
                   help='input, vep tsv file')
parse.add_argument('-p', dest='protout_file', required=True,
                   help="Path to vep_prot.tsv file")
parse.add_argument('-g', dest='genomeout_file', required=True,
                   help="Path to vep_genome.tsv file")
try:
    args = parse.parse_args()

except:
    exit('parsing of arguments failed!')

# parse args
infile_path = str(args.input_file)
outfile_path = str(args.output_file)
protfile_path = str(args.protout_file)
genomefile_path = str(args.genomeout_file)

# other important stuff
ref_fasta_path = "/mnt/share/data/genomes/GRCh37.fa" # pyfaidx reference
ref_fasta_file = Fasta(ref_fasta_path)
none_value= "NA"
chr_indices = {
    "1":0,
    "2":1,
    "3":2,
    "4":3,
    "5":4,
    "6":5,
    "7":6,
    "8":7,
    "9":8,
    "10":9,
    "11":10,
    "12":11,
    "13":12,
    "14":13,
    "15":14,
    "16":15,
    "17":16,
    "18":17,
    "19":18,
    "20":19,
    "21":20,
    "22":21,
    "X":22,
    "Y":23,
    "23":22
}

def querySingle(pos, chr):
    # input: reference position (starts @1!!)
    chr_index = chr_indices[chr]
    return str(ref_fasta_file[chr_index][pos-1]).upper()

def queryRange(start, end, chr):
    # input: reference position (starts @1!!)
    chr_index = chr_indices[chr]
    return str(ref_fasta_file[chr_index][start-1:end]).upper()

def reformatDuplication(notation):

    chromosome = notation.split(":")[0][3:]
    start = int(re.search("[0-9]{3,}", notation)[0])
    end = int(re.search("_[0-9]{3,}", notation)[0][1:]) if re.search("_[0-9]{1,2}", notation) else None
    bases = queryRange(start, end, chromosome) if end else querySingle(start, chromosome)
    # generate new notation w/ duplicationnn
    new = "chr%s:g.%s_%sdup%s" % (chromosome, start, end, bases) if end else "chr%s:g.%sdup%s" % (chromosome, start, bases)
    return new

def reformatVEST3Score(scores):

    # if entry is not NA -> put max is list, else none value
    return [max([float(y) for y in x.split("&")]) if type(x) is not float else none_value for x in scores]

def dropEntriesWithoutHGVSg(table):
    no_hgvs_mask = table["HGVSg"].isna()
    if sum(no_hgvs_mask) > 0:
        table = table[~no_hgvs_mask]
        print("FOUND ENTRIES W/O HGVSg NOTATION, DROPPING %s ROWS" % sum(no_hgvs_mask))

    return table

# read input & update w/ cosmic presence
vep_file = pd.read_csv(infile_path, delimiter="\t")
vep_file = dropEntriesWithoutHGVSg(vep_file)

# ## genome table (CGI, keys for annotation) - update dup notations w/ bases to sllow for them to match
out = []
HGVSg_notations = pd.Series(["chr"+x if not x.startswith("chr") else x for x  in vep_file["HGVSg"]])
HGVSg_notations.apply(lambda x: out.append(reformatDuplication(x)) if "dup" in x else out.append(x))  # reformat dups b/c why not
HGVSg_notations = pd.Series(out)  # convert to series for easy output
HGVSg_notations.rename("HGVSg")
HGVSg_notations.to_csv(genomefile_path, sep="\t", index=False, header=["HGVSg"])

# ## protein table (Civic, OncoKB)
hgvsp = []
prot_table = vep_file[["SYMBOL", "HGVSp"]]  # this throws warning, b/c its just a reference (run script w/ -W ignore)
prot_table.apply(lambda x: hgvsp.append(x["HGVSp"].split("p.")[1] if type(x["HGVSp"]) is not float else "NA"), axis=1)
prot_table["HGVSp"] = hgvsp # update in df
prot_table.to_csv(protfile_path, sep="\t", index=False)

# ## reformat VEST3_score if present
if "VEST3_score" in vep_file.columns:
    vep_file["VEST3_score"] = reformatVEST3Score(vep_file["VEST3_score"])

# drop unnecessary cols & write updated final_vep.tsv
vep_file.drop(columns=["SYMBOL", "HGVSp", "HGVSg"], axis=1, inplace=True)
vep_file.to_csv(outfile_path, sep="\t", index=False, na_rep="NA")
