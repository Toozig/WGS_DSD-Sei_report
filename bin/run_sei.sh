#!/bin/bash

#SBATCH --time=1-00:00:00     # Maximum runtime for the job (1 day)
#SBATCH -c16
#SBATCH -n 1                   # Number of tasks
#SBATCH --mem 100G             # Request 100 GB of memory

# module load hurcs bcftools

 
# # Check if the correct number of arguments is provided
# if [ $# -ne 3 ]; then
#     echo "Usage: $0 variant_file vcf_path_file output_dir"
#     exit 1
# fi

# # Input arguments
variant_file=$1
vcf_path_file=$2
output_dir=$3

# # activatin venv
source ~/new_venv/bin/activate

# Extract the base filename for BED file
bed_file=$(basename "$variant_file" ".csv").bed

# Prepare BED file from CSV file
awk -F',' 'BEGIN {OFS = FS} {print $2, $3 - 1, $3}' "$variant_file" |
    tr ',' '\t' |
    tail -n +2 > "$bed_file"

# Temporary VCF file
tmp_vcf="tmp.vcf"

# Merge VCF files based on the BED regions and normalize
bcftools merge -l "$vcf_path_file" -R "$bed_file" |
    bcftools norm -m +any -o "$tmp_vcf"


# Run variant effect prediction
sh 1_variant_effect_prediction.sh "$tmp_vcf" "hg38" "$output_dir" 

Clean up temporary files
rm -f "$bed_file" "$tmp_vcf"

# creates the effect score tsv
sh 2_varianteffect_sc_score.sh "$output_dir/chromatin-profiles-hdf5/tmp.ref_predictions.h5"  "$output_dir/chromatin-profiles-hdf5/tmp.alt_predictions.h5"  "$output_dir"
