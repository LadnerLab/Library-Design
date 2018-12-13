#!/bin/bash

input_dir="$1"
output_dir="$2"


for item in $(ls "$input_dir" | grep -v "out" ); do
    if [ ! -f "$input_dir"/"$item"_out_R_1 ]; then
        echo "$item"
        cp "$input_dir"/"$item" "$output_dir"
    fi
done
