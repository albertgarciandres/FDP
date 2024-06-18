#!/bin/bash

# Setting strict mode
set -euo pipefail

insertions=$(grep $condition results/Chimerics.tab | \
awk  -F $'\t'  ' $22 >= 0.39 && $23 == "mRNA" '  |  \
cut -f  1,3,4,6,8,15,21,23,22,27 | sort -u )

i=1

while IFS= read -r insertion
do
trinity=$(echo "$insertion" | cut -f2)
TEfam=$(echo "$insertion" | cut -f3)
stringtieID=$(echo "$insertion" | cut -f4)
geneID=$(echo "$insertion" | cut -f5)
TEconsensusID=$(echo "$insertion" | cut -f6)
CP=$(echo "$insertion" | cut -f7)
expr=$(echo "$insertion" | cut -f9)
id=$(echo "$insertion" | cut -f10)

grep -w "$trinity" results/RepeatMasker/${condition}Trinity.hits_bin90.fasta.out.gff | \
grep -w "$TEconsensusID" | mergeBed -d 25 > $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff

n=$(wc -l $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff | cut -f1 -d' ')

if [[ $n -eq 0 ]]; then

consensus=$(grep -w "$trinity" results/RepeatMasker/${condition}Trinity.hits_bin90.fasta.out.gff | \
cut -f 9 | cut -f2 -d' ' | tr -d '"' | cut -f2 -d':' | sort -u)

# DATA PROVIDED ON THE CONFIG FILE (UP TO CHANGE)
fam=$(grep "$consensus" $/TEs/TE_library_family.csv | cut -f 8)

if [[ $TEfam == $fam ]]; then
grep -w "$trinity" results/RepeatMasker/${condition}Trinity.hits_bin90.fasta.out.gff | \
grep -w "$consensus" | mergeBed -d 25 > $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff
n=$(wc -l $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff | cut -f1 -d' ')
echo $trinity $condition ${n}R $lengthTE >> n
lengthTE=$(awk '{sum += $3 - $2 + 1} END {print sum}' $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff )
fi

else
lengthTE=$(awk '{sum += $3 - $2 + 1} END {print sum}' $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff )

echo $condition $n $lengthTE >> n

fi

for seq in `seq 1 $n`
do
sed "${seq}q;d"  $condition/tmp/gff/${trinity}_${TEfam}.${i}.gff  > $condition/tmp/gff/${trinity}_${TEfam}.${seq}.${i}.gff 

bedtools getfasta -fi results/RepeatMasker/${condition}Trinity.hits_bin90.fasta.masked \
    -bed $condition/tmp/gff/${trinity}_${TEfam}.${seq}.${i}.gff | \
sed "s/$trinity:.*/${id}_${condition}_${TEfam}_${lengthTE}_${trinity}/g" > $condition/tmp/fasta/${trinity}_${TEfam}.${seq}.${i}.fasta 
getorf -sequence $condition/tmp/fasta/${trinity}_${TEfam}.${seq}.${i}.fasta  \
    -outseq $condition/tmp/fasta/${trinity}_${TEfam}.ORF.${seq}.${i}.fasta

# NEED TO SPECIFY WHERE THE PFAM FILES ARE ON THE WORKFLOW
pfam_scan.pl -fasta $condition/tmp/fasta/${trinity}_${TEfam}.ORF.${seq}.${i}.fasta  \
    -dir pfamFiles -cpu 2 \
    -outfile $condition_${trinity}_${TEconsensusID}_${group}.${seq}.${i}.pfam.output
done
i=$(( i + 1 ))

done <<< "$insertions"
