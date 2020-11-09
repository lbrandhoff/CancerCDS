"""
vep path -  /mnt/share/opt/ensembl-vep-release-98.3/vep
vep cache directory - /tmp/local_ngs_data/ensembl-vep-98/
"""
import argparse, os

parse = argparse.ArgumentParser(description=".")

parse.add_argument('-i', dest='input_file', required=True,
                   help='Path to input vcf')
parse.add_argument('-o', dest='output_name', required=False,
                   help="Path to output vcf, default ./annotated_'input'.vcf ")
parse.add_argument('-t', dest='threads', required=False,
                   help="# of threads to use (default: 1)")
parse.add_argument('-b', dest='buffer', required=False,
                   help="buffer size of vep (default: 5000)")
try:
    args = parse.parse_args()

except:
    exit('run_vep.py: parsing of arguments failed!')

# paths to files
input_vcf = str(args.input_file)
output_vcf = str(args.output_name)
num_threads = args.threads
buffer_size = args.buffer

# if no output given -> give standard suffix
if output_vcf == "None":
    output_vcf = "annotated_%s" % input_vcf

# num of threads, buffer size
num_threads = int(num_threads) if num_threads else 1
buffer_size = int(buffer_size) if buffer_size else 5000

# path to vep relevant stuff
vep_path = "/mnt/users/ahbranl1/data_vep/ensembl-vep-release-98.3/vep"
vep_cache = "/mnt/users/ahbranl1/data_vep/ensembl-vep-cache/cache/"
dbNSFP_file = "/mnt/users/ahbranl1/data_vep/dbNSFP/dbNSFP_hg19_3.5.gz"

reference_fasta_path = "/mnt/share/data/genomes/GRCh37.fa"  # reference genome fasta
db_path = "/mnt/share/data/dbs/"


# list for accumulating arguments
args = []

# general args
args.append("--format vcf -i %s" % input_vcf)  # input file + format
args.append("-o %s --vcf" % output_vcf)  # output file + format
args.append("--cache --dir_cache %s" % vep_cache)
args.append("--species homo_sapiens --assembly GRCh37")  # species + assembly type
args.append("--fork %s --offline --fasta %s" % (num_threads, reference_fasta_path))  # num of threads + path to reference fasta
args.append("--force_overwrite --no_stats")   # forces overwrite + no stats html

# annotation args
args.append("--af --af_gnomad")  # 1000g & gnomad frequence
args.append("--sift s --polyphen s")  # damage predictions
args.append("--hgvs --hgvsg")  # HGVS notation
args.append("--buffer_size %s" % buffer_size)  # size of vars for each thread

# custom annos (added to args)
args.append("--custom  %sgnomAD/gnomAD_genome_r2.1.1.vcf.gz,gnomADg,vcf,exact,0,AF" % db_path)

# plugin annos
args.append("--plugin dbNSFP,%s,LRT_score,GERP++_RS,Polyphen2_HVAR_score,SIFT_score,CADD_phred,Eigen-raw,MetaSVM_score,VEST3_score,MetaLR_score" % dbNSFP_file)

# customise fields (added to args)
fields = ["Consequence","PolyPhen", "SIFT", "LRT_score", "GERP++_RS",
          "gnomAD_AF", "AF", "SYMBOL", "HGVSp", "HGVSg", "CADD_phred", "Eigen-raw","MetaSVM_score","VEST3_score",
          "MetaLR_score"]


args.append("--fields " + ",".join(fields))

# collapse to string for bash
command_str = vep_path + " " + " ".join(args)

# exec
os.system(command_str)
