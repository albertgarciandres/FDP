#!/bin/bash

# Using BEDOPS to turn gff to bed and using awk to modify bed to bed6 format.
# sed used to have a pretty ID field
#
# BED& format:
# chr start stop . ID . strand



gff2bed < $in_file | awk 'BEGIN {FS="\t";OFS="\t"}{print $1, $2, $3, $10, $4, $6}' | sed 's/ID=//' > $out_file

