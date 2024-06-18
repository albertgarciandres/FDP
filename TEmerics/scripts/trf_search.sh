#!/bin/bash

# Setting strict mode
set -euo pipefail

display_usage(){
    echo "DESCRIPTION."
    echo ""
    echo "Usage trf_search.sh [[-g in-gff] [-f in-masked-fa] [-o out-dir] [-m match-p] [-i indel-p] [-ms min-score] [-mp max-period] [-mw match-weigth] [-mmw mismatch-weigth] [-iw indel-weigth] | [-h]]"
    echo ""
    echo "Args:"
    echo "      in-gff: "
    echo "      in-fa: "
    echo "      out-dir: "
    echo "      match-p: "
    echo "      indel-p: "
    echo "      min-score: "
    echo "      max-period: "
    echo "      indel-weigth: "
    echo "      match-weigth: "
    echo "      mismatch-weigth: "
    echo ""
#    echo "Usage trf_search.sh [[OPTIONS] [-g in-gff] [-f in-masked-fa] [-o out-dir]| [-h]]"
#    echo ""
#    echo "OPTIONS"
#    echo "      -m   | --match-p         (Default 80)"
#    echo "      -i   | --indel-p         (Default 10)"
#    echo "      -ms  | --min-score       (Default 20)"
#    echo "      -mp  | --max-period      (Default 15)"
#    echo "      -iw  | --indel-weigth    (Default 5)"
#    echo "      -mw  | --match-weigth    (Default 2)"
#    echo "      -mmw | --mismatch-weigth (Default 3)"
#    echo ""
}

### Make sure input files present or if help

if [[($# -ne 10) || ($1 == "-h")]]
then
    display_usage
    exit
fi

### Get CLI options

while [$# -gt  0];do
  case $1 in
        -g | --in-gff)  shift gff=$1 ;; # ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.out.gff
        -f | --in-fa)   shift   fa=$1 ;; # ${DIR}/RepeatMasker/${strain}_${tissue}/Trinity.hits_bin90.fasta.masked
        -o | --out-dir)  shift   out_dir=$1  ;; # ${DIR}/trf/${strain}_${tissue}/tmp/
        -m | --match-p)    shift pm=$1 ;;   # 80
        -i | --indel-p)    shift pi=$1 ;;   # 10
        -ms | --min-score)  shift min_score=$1  ;;  # 20
        -mp | --max-period) shift   max_period=$1   ;;  # 15
        -iw | --indel-weigth)   shift   indel_w=$1  ;;  # 5
        -mw | --match-weigth)   shift   match_w=$1  ;;  # 2
        -mmw | --mismatch-weigth)   shift   mm_w=$1 ;;  # 3
        -h | --help)    display_usage   exit 0 ;;
        *)  display_usage   exit 1
    esac
    shift
done


### Main 
while read transcript_ID
	do
        inter_dir=$out_dir/$transcript_ID
        # Get transcript GFF
		awk -v ID="$transcript_ID" ' $1 == ID ' $gff > $inter_dir.intermediates.gff
        cat $inter_dir.intermediates.gff | cut -f 1,4,5,9 | awk '{print $1"\t"$2"\t"$3"\t"$5}' FS="[\"]" OFS="\t" | cut -f 1,2,3,5 > $inter_dir.gff
		
        # Sort and merge unique entries
		> $inter_dir.collapse.gff 
        cat $inter_dir.gff | cut -f 4 | sort | uniq | while read C
			do
				awk -v C=${C} '($4==C)' $inter_dir.gff | sort -t $'\t' -k1,1 -k2,2n | bedtools merge  -c 4 -o distinct >> $inter_dir.collapse.gff 
			done
        # Remove "Motif" keywords (like in vim)
		sed -i "s/Motif://g" $inter_dir.collapse.in_gff
        
        # Generate FASTA from GFF
		bedtools getfasta -fi $fa -bed $inter_dir.collapse.gff -name -fo $inter_dir.collapse.fa

        # Run Tandem Repeat Finder (trf)
		trf $inter_dir.collapse.fa $match_w $mm_w $indel_w $pm $pi $min_score $max_period -h -d

        # Classify/filter trf results
		i=0
		while read line ; do
			if echo "$line" | grep -q "Sequence"; then
				i=$(( $i + 1 ))
				echo $line > $inter_dir.trf.$i
				continue
			fi
			echo "$line" >> $inter_dir.trf.$i
		done < <(cat $inter_dir.collapse.fa.$match_w.$mm_w.$indel_w.$pm.$pi.$min_score.$max_period.dat)

        # Merge results
		> $inter_dir.collapse.fa.trf
		for seq in `seq 1 $i`
			do
				merge=$(cat $inter_dir.trf.$seq | grep "^[0-9]" | cut -f 1,2 -d' ' | sed  -e 's/^/TE /'  | tr " " "\t" | mergeBed | cut -f2,3 | tr "\t" " ")
				info=$(cat $inter_dir.trf.$seq | grep -E "^Sequence|^Parameter")
				echo -e "$info\n$merge"
			done >> $inter_dir.collapse.fa.trf

	done < <(cut -f1 $gff | grep -v "^#" | sort -u )

