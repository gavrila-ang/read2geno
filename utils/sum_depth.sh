#!/bin/bash

# Usage : ./sum_columns.sh your_file.txt output.txt
# def calculate_raw_snp_read_depth_summed():
#     # Sum depths of the variants contained within the STICH positions file across all samples
#     subprocess.run(
#         f"{pipeline_config['code']}/utils/sum_depth.sh "
#         f"{raw_prefix}.depth",
#         f"{raw_prefix}_depth_summed.txt",
#         shell=True, check=True
#     )



# Check for input and output file arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

input_file="$1"
output_file="$2"

# AWK script to sum columns starting from the 3rd, preserving the first two columns
awk -F'\t' '{
    sum = 0
    for (i = 3; i <= NF; i++) {
        if ($i ~ /^-?[0-9]+$/)
            sum += $i
        else if ($i == "")
            sum += 0
    }
    print $1 "\t" $2 "\t" sum
}' "$input_file" > "$output_file"
