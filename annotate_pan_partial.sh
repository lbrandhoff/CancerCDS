#!/bin/bash
GRCh37_fa=/mnt/share/data/genomes/GRCh37.fa

# scripts
dir=./scripts/
intogen_annot_py=${dir}getIntogenAnnot_pan.py


temp_dir=$1  # working directory

echo "Adding intogen annotation"
python3 $intogen_annot_py -i ${temp_dir}vep_prot.tsv -o ${temp_dir}intogen.tsv

# combine data
echo "Combining data"
paste ${temp_dir}vep_gdna.tsv ${temp_dir}reformat_vep.tsv ${temp_dir}intogen.tsv ${temp_dir}hotspot.tsv ${temp_dir}gene_drug.tsv ${temp_dir}vep_prot.tsv > ${temp_dir}pan_annotation.tsv


echo "Removing temporary data"
rm ${temp_dir}/trash_me.vcf
#rm ${temp_dir}/annotated_trash_me.vcf
rm ${temp_dir}/vep_gdna.tsv
rm ${temp_dir}/vep_prot.tsv
rm ${temp_dir}/vep.tsv
rm ${temp_dir}/reformat_vep.tsv
rm ${temp_dir}/hotspot.tsv
rm ${temp_dir}/gene_drug.tsv
rm ${temp_dir}/intogen.tsv

echo "Finished"


