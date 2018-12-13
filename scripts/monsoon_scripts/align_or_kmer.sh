#!/bin/bash
input_dir="$1"
output_dir="$2"

start_path=$(pwd)
cd "$output_dir"
outdir_path=$(pwd)
cd "$start_path"/"$input_dir"

files=$(ls | grep -v "out" )
for item in ${files[@]}; do
    num_align=$(grep -c ">" "$item".aligned_out)
    num_kmer=$(grep -c ">" "$item"_out_R_1 )
    
    if [ "${num_align:-0}" -lt "$num_kmer" ]; then
        cp "$item".aligned_out "$outdir_path"
    else
        cp "$item"_out_R_1 "$outdir_path"
    fi
done
