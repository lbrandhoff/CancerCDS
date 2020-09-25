#!/bin/bash
GRCh37_fa=/mnt/share/data/genomes/GRCh37.fa

dir=./scripts/
gDNA2vcf_py=${dir}gDNAtoVCF.py
run_vep_py=${dir}run_vep.py
vcf2tsv_py=${dir}vcf_to_tsv.py
vep2prot_py=${dir}vep_tsv_reformat.py
hotspot_annot_py=${dir}getHotspotAnnot.py
drug_annot_py=${dir}getGeneDruggabilityAnnot.py
intogen_annot_py=${dir}getIntogenAnnot_pan.py


# read from input
vcf_in=$1
out_dir=$2

echo "------------------------------------------------------------------"
echo "----Annotating ${vcf_in} with Pan-Cancer"
echo "------------------------------------------------------------------"

temp_dir=${out_dir}  # dir to store temp files in


mkdir -p ${temp_dir}  # makes tempdir if not present

#echo "Sorting .vcf"
#/mnt/share/opt/ngs-bits-current/VcfSort -in ${vcf_in} -out ${temp_dir}trash_me.vcf  # sort
#echo "Running VEP annotation"
#python3 $run_vep_py -i ${temp_dir}trash_me.vcf -o ${temp_dir}annotated_trash_me.vcf -t 5 -b 1000
echo "Converting .vcf to .tsv file"
python3 $vcf2tsv_py -i ${temp_dir}annotated_trash_me.vcf -o ${temp_dir}vep.tsv
echo "Extracting protein and genome notations"
python3 -W ignore $vep2prot_py -i ${temp_dir}vep.tsv -o ${temp_dir}reformat_vep.tsv -p ${temp_dir}vep_prot.tsv -g ${temp_dir}vep_gdna.tsv

echo "Adding hotspot annotation"
python3 -W ignore $hotspot_annot_py -i ${temp_dir}vep_prot.tsv -o ${temp_dir}hotspot.tsv

echo "Adding gene druggability annotation"
python3 $drug_annot_py -i ${temp_dir}vep_prot.tsv -o ${temp_dir}gene_drug.tsv

echo "Adding intogen annotation"
python3 $intogen_annot_py -i ${temp_dir}vep_prot.tsv -o ${temp_dir}intogen.tsv

# combine data
echo "Combining data"
paste ${temp_dir}vep_gdna.tsv ${temp_dir}reformat_vep.tsv ${temp_dir}intogen.tsv ${temp_dir}hotspot.tsv ${temp_dir}gene_drug.tsv ${temp_dir}vep_prot.tsv > ${temp_dir}pan_annotation.tsv


echo "Removing temporary data"

#rm ${temp_dir}trash_me.vcf
#rm ${temp_dir}/annotated_trash_me.vcf
rm ${temp_dir}vep_gdna.tsv
rm ${temp_dir}vep_prot.tsv
rm ${temp_dir}vep.tsv
rm ${temp_dir}reformat_vep.tsv
rm ${temp_dir}hotspot.tsv
rm ${temp_dir}gene_drug.tsv
rm ${temp_dir}intogen.tsv
echo "Finished"


