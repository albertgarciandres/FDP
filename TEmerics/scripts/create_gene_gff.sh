#!/bin/bash/

# -i in_file
# -o gene_gff
# -l gene_trans_list

cut -f1,4 $in_file | sort -u > $gene_trans_list
cut -f2 $gene_trans_list | sed 's/^/ID=/' > $gene_list
grep -f $gene_list $ref_gff | sed -i 's/;.*//' > $gene_gff
