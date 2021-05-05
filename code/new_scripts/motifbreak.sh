#!/usr/bin/env bash

## This script loops through all files in the directory
## and outputs the commands necessary to launch the slurm array jobs to be executed in parallel

wd="/data/projects/punim0586/dvespasiani/Archaic_introgression_in_PNG"
input_dir="$wd/Motifbreak/input_files"
output_dir="$wd/Motifbreak/output_files"
script_dir="$wd/scripts"

if [ ! -d "$output_dir/hocomoco" ]; then
  mkdir -p "$output_dir/hocomoco" &&  mkdir -p "$output_dir/jaspar" ;
fi

for file in "$input_dir"/*.bed; do

basename_file="$(echo $(basename $file)| cut -f 1 -d '.')"

echo "R --vanilla -f "$script_dir/motifbreak_analysis.R" \
 --args input=$file output_jaspar=$output_dir/jaspar/${basename_file}_jaspar.txt \
 output_hocomoco=$output_dir/hocomoco/${basename_file}_hocomoco.txt" 
done



