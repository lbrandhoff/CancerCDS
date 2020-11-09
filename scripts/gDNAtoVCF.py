import re
from pyfaidx import Fasta
import argparse

parse = argparse.ArgumentParser(prog="gDNAtoVCF Converter",
                                description="Requires newline separated gDNA notations w/o RefSeq sequence "
                                            "e.g. chr2:g.29443600_29443601delinsCG \\n chr11:g.108164053dupT \\n "
                                            "chr7:g.140453134T>C")

parse.add_argument('-i', dest='input_file', required= True,
                   help='Path to file containing gDNA notations')
parse.add_argument('-r', dest='reference_fasta', required= False,
                   help='Path to reference genome fasta file')
parse.add_argument('-o', dest='output_name', required= True,
                   help='Path w/ name of output vcf file')
try:
    args = parse.parse_args()

except:
    exit('gDNAtoVCF Converter: parsing of arguments failed!')

# paths to files
input_path = str(args.input_file)
ref_fasta_path = str(args.reference_fasta)
output_path = str(args.output_name)

# check if output has correct ending attached
output_path = output_path if output_path.endswith(".vcf") else output_path + ".vcf"

# if path given - use, else use GRCh37 reference
ref_fasta_path = "/mnt/share/data/genomes/GRCh37.fa" if ref_fasta_path == "None" else ref_fasta_path

# open files
input_file = open(input_path) # read
output_file = open(output_path, "w")
ref_fasta_file = Fasta(ref_fasta_path)  # check if chr order is same is in chr_indices - else change!

# dict for correct index accession in fasta file
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

current_var = []  # elements in specific order: chr, pos, id, ref, alt

# write header to output file
header_line = "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO\n"])
output_file.write(header_line)

# count of lines
total = 0
conv = 0
unconv = 0

# vars for variant info
chromosome = "0"
pos = 0
ref = "X"
alt = "X"

def querySingle(pos):
    # input: reference position (starts @1!!)
    if chromosome == '':
        print(line)
        exit()
    chr_index = chr_indices[chromosome]
    return str(ref_fasta_file[chr_index][pos-1]).upper()

def queryRange(start, end):
    # input: reference position (starts @1!!)
    chr_index = chr_indices[chromosome]
    return str(ref_fasta_file[chr_index][start-1:end]).upper()


# convert lines
for line in input_file:
    total += 1
    line = line.strip()

    # snps
    if ">" in line:

        conv += 1
        split = re.split("\.|>", line)
        chromosome = split[0][3:-2]
        pos = split[1][:-1]
        ref = split[1][-1]
        alt = split[2]

    # duplications
    elif "dup" in line:

        conv += 1
        split = re.split("\.|_|dup", line)

        chromosome = split[0][3:-2]
        pos = int(split[-2])  # pos: last position of duplication is first reference pos
        ref = split[-1]  # for now ref: bases after dup

        # if bases present after dup no query to reference
        if ref:
            ref = ref[-1]  # ref only last base of dupl in vcf!
            alt = ref + split[-1]  # ATT dupl -> TATT alt in vcf

        # else, query reference
        else:
            if len(split) == 3:  # only one base dup
                # -1: pos from ref, ref starts @ 1, faidx starts @ 0
                ref = querySingle(pos)
                alt = 2*ref

            elif len(split) == 4:  # more base dup
                start_refpos = int(split[1])
                ref = queryRange(start_refpos, pos)  # aGT (GT dup) -> ref: GT
                alt = ref[-1] + ref  # alt: TGT
                ref = ref[-1]  # ref: T


            else:
                print("UNCONVERTABLE 'dup' LINE, skipping;", line)
                unconv += 1
                continue  # skip writing to file

    # all ins lines
    elif re.search("[0-9]ins", line):

        conv += 1

        split = re.split("\.|_|ins", line)

        chromosome = split[0][3:-2]
        pos = int(split[1])
        ref = querySingle(pos)
        alt = ref + split[-1]

    # all delins e.g. delTTinsAGT, delinsATTG
    elif re.search("del[A,C,G,T]*ins", line):

        conv += 1

        split = re.split("\.|_|del|ins", line)

        chromosome = split[0][3:-2]
        pos = int(split[1])-1
        end_refpos = int(split[-3])
        ref = queryRange(pos, end_refpos)
        alt = ref[0] + split[-1]


    # del only lines - ins & delins already catched
    elif re.search("del[^ins]*", line):

        conv += 1
        split = re.split("\.|_|del", line)
        chromosome = split[0][3:-2]
        pos = int(split[1]) - 1

        # if bases present after del
        if split[-1]:
            alt = querySingle(pos)
            ref = alt + split[-1]  # cAG w/ AG deleted -> ref:CAG, alt: C

        # if no bases after del
        else:
            # if only one pos del
            if len(split) == 3:
                end_refpos = int(split[1])

                ref = queryRange(pos, end_refpos)
                alt = ref[0]


            # multiple bases deleted
            elif len(split) == 4:

                end_refpos = int(split[2])

                ref = queryRange(pos, end_refpos)
                alt = ref[0]

            else:
                print("UNCONVERTABLE 'del' LINE, skipping;", line)
                unconv += 1
                continue  # skip writing to file

    elif "inv" in line:
        print("Inversion not implemented, skipping line:", line)
        unconv += 1

    else:  # count unconvertable line
        print("UNCONVERTABLE LINE, skipping:", line)
        unconv += 1
        continue  # skip writing to file

    # create line write
    current_var = [chromosome, str(pos), ".", ref, alt, ".", ".", ".\n"]

    # write line to file
    output_file.write("\t".join(current_var))

output_file.close()

# print counters
print("# total lines in file:", total)
print("# converted lines:", conv)
print("# unconvertable lines:", unconv)
