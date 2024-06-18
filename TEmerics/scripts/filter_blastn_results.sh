#!/bin/bash

# Setting strict mode
set -euo pipefail

display_usage(){
    echo "Filter BLASTn results by coverage and seq-ref relationship."
    echo ""
    echo "Usage: filter_blastn_results.sh [[-l in-lst] [-t in-tmap] [-s stats] [-o out-dir] | [-h]]"
    echo ""
    echo "Args:"
    echo "     in-lst: input e-value list (.list)"
    echo "     in-tmap: input tab delimited file list (.tmap)"
    echo "     stats: statistics output directory path."
    echo "     out-dir: output path and file name(.tsv)"
    echo ""
    echo "The allowed relationships are:"
    echo ""
    echo "    - Complete, exact match of intron chain. (=)"
    echo "    - Multi-exon with at least one junction match. (j)"
    echo "    - Contained in reference (intron compatible). (k)"
    echo "    - Retained intron(s), all introns matched pr retained. (m)"
    echo "    - Retained intron(s), not all introns matched/covered. (n)"
    echo "    - Other same strand overlap with reference exons. (o)"
    echo "    - Exonic overlap on the opposite strand. (x)"
    echo "    - Contains a reference within its intron(s). (y)"
    echo ""
}


### Make sure all parameters are present or if help

if [[($# -ne 3) || ($1 == "-h")]]
then
    display_usage
    exit
fi


### Get CLI options

# Stats default val
stats=$PWD

while [$# -gt 0];do
    case $1 in
        -l | --list)     shift in_lst=$1 ;;
        -t | --tmap)     shift tmap=$1 ;;
        -s | --stats)    shift stats=$1   ;;
        -o | --out-tbl)    shift out_table=$1   ;;
        -h | --help)    display_usage   exit 0   ;;
        *)  display_usage   exit 1
    esac
    shift
done

cov_names="${stats}/names_filtered_bin90_100.tsv"
cov_table="${stats}/filtered_results_bin90_100.tsv"

### Main
grep -E "Bin_90|Bin_100" ${in_lst} | cut -f2,3  > "$cov_names"

while read transcript
do
    trinity_ID=$(echo "$transcript" | cut -f1)
    stringtie_ID=$(echo "$transcript" | cut -f2)

    ref_gene_ID=$(awk -v ID="$stringtie_ID" '$5 == ID' "$tmap" | cut -f1 | sort -u )
    ref_trans_ID=$(awk -v ID="$stringtie_ID" '$5 == ID' "$tmap" | cut -f2 | sort -u )
    relation_type=$(awk -v ID="$stringtie_ID" '$5 == ID' "$tmap" | cut -f3 | sort -u )

    echo -e "$trinityi_ID\t$stringtie_ID\t$relation_type\t$ref_gene_ID\t$ref_trans_ID" >> "$cov_table"
done < "$cov_names"


# awk '$3 == "=" || $3 == "j" || $3 == "k" || $3 == "m" || $3 == "n" || $3 == "o" || $3 == "x" || $3 == "y"' {input.tab}  | cut -f1 |  sort -u > ${out_file} 
awk '$3 ~ /[=|j|k|m|n|o|x|y]/' "$cov_table" | cut -f1 | sort -u > "$out_table"

if ["$stats" == $PWD];
then
    rm "$cov_names"
    rm "$cov_table"
fi
