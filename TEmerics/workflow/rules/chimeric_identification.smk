###############################################################################
### WHAT IT DOES?
###############################################################################
#
#   - Run Repeatmasker to detect TEs
#   - Run Tandem Repeat Finder to filter TEs
#

###############################################################################
### Imports
###############################################################################

import pandas as pd
from snakemake.utils import validate

from pathlib import Path

###############################################################################
### Configuration validation
###############################################################################

validate(config, Path("../../config/config_schema.json"))

###############################################################################
### Paths configuration
###############################################################################

OUT_DIR = Path(config["out_dir"])
SCRIPTS_DIR = Path(config["scripts_dir"])
ENV_DIR = Path(f"{workflow.basedir}/envs")

LOCAL_LOG = Path(config["local_log"])

###############################################################################
### Including functions
###############################################################################

include: "common.smk"


###############################################################################
### Global configuration
###############################################################################


localrules:
    finish_chimeric_id,


###############################################################################
### Local rules
###############################################################################


rule finish_chimeric_id:
    input:
        out = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.out",
        masked = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.masked",
        gff_out = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.out.gff",
        trans_trf = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "intermediates" / "{TRANS_ID}.collapse.fa.trf",
        table = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "percentageSRR.tab",
        fa = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.fasta",
        table= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics_filtered_unique.blastn",

###############################################################################
### Annotate TEs with Repeatmasker
# -no_is skips bacterial insertion element check
# -nolow does not mask low complexity DNA or simple repeats
# -norna does not mask small RNA genes while still masking SINEs
# -xsmall repetitive regions are not masked but in lowercase
# -s Slow search (0-5% more sensitive, 2-3 time slower
###############################################################################

rule annotate_te_repeatmasker:
    input:
        te_lib = config["te_lib"],
        fasta=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_filtered.fasta",
    output:
        fa = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_masked.fasta",
        gff = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_masked.gff",
    params:
        cutoff = config["cutoff"],
        out_dir= OUT_DIR / "chimeric_id" / "{SAMPLE}",
        out_fa = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.masked",
        out_gff = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.out.gff",
    resources:
        processors = 4,
    log:
        LOG_DIR / "annotate_{SAMPLE}_te_repeatmasker.log"
    conda:
        ENV_DIR / "repeatmasker.yalm"
    shell:
        """
        (RepeatMasker {input.fasta} \
        -engine rmblast \
        -lib {input.te_lib} \
        -dir {params.out_dir} \
        -cutoff {params.cutoff} \
        -pa {resources.processors} \
        -norna \
        -nolow \
        -no_is \
        -xsmall \
        -gff \
        -s && mv {params.out_gff} {output.gff} && mv {params.out_fa} {output.fa} \
        ) &> {log}
        """

###############################################################################
### Process RepeatMasker GFF output file
###############################################################################

rule find_tandem_repeats_per_transcript:
    input:
        fa = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_masked.fasta",
        gff = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_masked.gff",
        script = SCRIPTS_DIR / "trf_search.sh",
    output:
        final = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "intermediates" / "{TRANS_ID}.collapse.fa.trf"
    params:
        match_p=config["PM"],
        indel_p=config["PI"],
        min_score=config["min_score"],
        indel_weight=config["indel_w"],
        max_period=config["max_period"],
        mism_weight=config["mism_weight"],
        match_weight=config["match_weight"],
        out_dir= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "intermediates",
    log:
        LOG_DIR / "find_{SAMPLE}_tandem_repeats_{TRANS_ID}.log"
    conda:
        ENV_DIR / "trf.yalm"
    shell:
        """
        (bash {input.script} \
        --in-gff {input.gff} \
        --in-fa {input.masked} \
        --out-dir {params.out_dir} \
        --match-weight {params.match_weight} \
        --mismatch-weight {params.mism_weight} \
        --indel-weight {params.indel_weight} \
        --match-p {params.match_p} \
        --indel-p {params.indel_p} \
        --min-score {params.min_score} \
        --max-period {params.max_period} \
        ) &> {log}
        """

###############################################################################
### Get TRF table per transcript per condition
###############################################################################

rule get_trf_table_per_trans:
    input:
        trans_trf = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "intermediates" / "{TRANS_ID}.collapse.fa.trf"
        script = SCRIPTS_DIR / "trf_search.py",
    output:
        table = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "{TRANS_ID}.percentageSRR.tab",
    params:
        # Not sue this will work. Maybe you need to refactor the python script
        trans_id = [wildcards.TRANS_ID],
    log:
        LOG_DIR / "get_{SAMPLE}_trf_{TRANS_ID}_table.log",
    conda:
        ENV_DIR / "python.yaml"
    shell:
        """
        (python {input.script} \
        {input.trans_trf} \
        {params.trans_id} \
        > {output.table} \
        ) &> {log}
        """

###############################################################################
### Merge TRF tables
###############################################################################

rule merge_trf_results:
    input:
         # Pretty sure this is not the actual way to do it. Maybe you have to
         # use a lambda function.
        table = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "{TRANS_ID}.percentageSRR.tab",
    output:
        table = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "percentageSRR.tab"
    log:
        LOG_DIR / "merge_{SAMPLE}_trf_tables.log"
    shell:
        """
        (touch {output.table} && cat {input.file} >> {output.file}) &> {log}
        """

###############################################################################
### WHAT IT DOES?
###############################################################################
#
#   - Chimeric identification
#
###############################################################################
### Get a list of chimeric entries per condition
###############################################################################

rule get_chimeric_entries:
    input:
        gff = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.out.gff",
    output:
        lst = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.lst",
    log:
        LOG_DIR / "get_chimeric_in_{SAMPLE}.log",
    shell:
        """
        (grep -v "^#" {input.gff} | cut -f1 | sort -u > {output.lst}) &> {log}
        """


###############################################################################
### Convert list to a fasta file
###############################################################################

rule get_chimeric_fasta:
    input:
        masked = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "reads_filtered.fasta.masked",
        lst = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.lst"
    output:
        fa = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.fasta",
    log:
        LOG_DIR / "retrieve_{SAMPLE}_chimeric_seqs.log",
    conda:
        ENV_DIR / "seqtk.yalm"
    shell:
        """
        (seqtk subseq {input.file} {input.lst} > {output.file} )&> {log}
        """


###############################################################################
### Run BLAST: Retrieved chimeric sequences
# Output format specifications in https://www.metagenomics.wiki/tools/blast/blastn-output-format-6
###############################################################################

rule run_blast_chimeric:
    input:
        query = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.fasta",
        db_name = OUT_DIR / "guided_assembly" / "all_seqs_database.fasta",
    output:
        blast = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.blastn",
    params:
        e_val = config["e_val_chimerics"],
    resources:
        threads = 4,
    log:
        LOG_DIR / "run_blastn_{SAMPLE}_chimeric.log",
    conda:
        ENV_DIR / "blastn.yalm"
    shell:
        """
        (blastn \
        -query {input.query} \
        -db {input.db_name} \
        -out {output.blast} \
        -outfmt 6 \
        -evalue {params.e_val} \
        -num_threads {resources.threads} \
        ) &> {log}
        """


###############################################################################
### Generate percentages table (filter by length, coverage and bitscore)
# MODIFY SCRIPT TO ACCOMODATE THE ACTUAL COLNAMES 
# MAKE SURE SCOV IS NOT REQUIRED AND CAN BE COMPUTED
###############################################################################

rule generate_percentage_table:
    input:
        blast = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics.blastn",
        ref_trans = "NO_FILE??",
        script = SCRIPTS_DIR / "calculate_percentage_v2.R"
    output:
        table= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics_filtered.blastn",
    log:
        LOG_DIR / "generate_{SAMPLE}_percent_table.log"
    conda:
        ENV_DIR / "r.yalm"
    shell:
        """
        (Rscript {input.script} \
        --in-file {input.blast}
        --ref-trans {input.ref_trans} \
        --out-file {output.table}
        --verbose \
        ) &> {log}
        """


###############################################################################
### Filter BLAST information table
###############################################################################

rule filter_percent_table:
    input:
        table= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics_filtered.blastn",
        script = SCRIPTS_DIR / "filter_percent_table.sh"
    output:
        table= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics_filtered_unique.blastn",
    log:
        LOG_DIR / "filter_percent_table_{SAMPLE}.log"
    shell:
        """
        bash {input.script} {input.table} {output.table}) &> {log}
        """

###############################################################################
### Create gene list for chimeric transcripts
###############################################################################

rule create_gene_gff:
    input:
        table= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "chimerics_filtered_unique.blastn"
        ref = "dmel-all-no-analysis-r6.31.gene.gff"
        script = SCRIPTS_DIR / "create_gene_gff.sh"
    output:
        gff = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "gene.gff"
        list = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "gene_trans.lst"
    log:
        LOG_DIR / "create_{SAMPLE}_gene_list.log"
    shell:
        """
        (bash {input.script} \
        -i {input.table} \
        -gff {input.ref} \
        -o {output.gff} \
        -l {output.list} \
        ) &> {log}
        """

###############################################################################
### Convert GFF to BED6
###############################################################################

rule convert_gene_gff_to_bed6:
    input:
        gff= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "gene.gff"
        script = SCRIPTS_DIR / "gff_to_bed6.sh"
    output:
        bed = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "gene.bed"
    log:
        LOG_DIR / "convert_{SAMPLE}_gene_gff_to_bed.log"
    conda:
        ENV_DIR / "bedops.yalm"
    shell:
        """
        (bash {input.script} \
        -i {input.gff} \
        -o {output.bed} \
        ) &> {log}
        """
###############################################################################
### Get distance between closest genes
# -d add distance in additional field (if overlapping, val = 0)
# -N make sure different names
# -s account for strandness
# If no closest feature found: · -1 -1 · -1 · · 1 In half-right side
###############################################################################

rule get_distance_to_closest_gene:
    input:
        bed= OUT_DIR / "chimeric_id" / "{SAMPLE}" / "gene.bed"
    output:
        bed = OUT_DIR / "chimeric_id" / "{SAMPLE}" / "gene_distances.bed"
    log:
        LOG_DIR / "get_{SAMPLE}_dist_closest_gene.log"
    conda:
        ENV_DIR / "bedtools.yalm"
    shell:
        """
        (bedtools closest \
        -a {input.bed} \
        -b {input.bed} \
        -d \
        -N \
        > {output.bed} \
        ) &> {log}
        """


###############################################################################
### Extract chimerics that are not fussioned
###############################################################################

rule extract_not_fussioned_chim:
    input:
        tab = "check_fussion/output/{SAMPLE}/resultDistance.tab"

    output:
        lst = "check_fussion/output/{SAMPLE}/genes_no_fussioned.lst"

    log:
        LOG_DIR / "subworkflow4/extract_not_fussioned_chim/{SAMPLE}.log"

    shell:
        """
        (awk ' $3 == 1 ' {input.tab} | cut -f1 | sort -u > {output.lst}) &> {log} 
        """


###############################################################################
### Process the list with the fasta
###############################################################################

rule process_fasta:
    input:
        fasta = "trinity/trinity_{SAMPLE}/Trinity.hits_bin90.fasta",
        lst = "check_fussion/output/{SAMPLE}/genes_no_fussioned.lst"
    
    output:
        fa = "minimap2/{SAMPLE}/Trinity.hits_bin90.clean.fasta"

    log:
        LOG_DIR / "subworkflow4/process_fasta/{SAMPLE}.log"
    
    conda:
        ENV_DIR / "seqtk.yalm"

    shell:
        """
        (seqtk subseq {input.fasta} \
        {input.lst} > \
        {output.fa} \
        ) &> {log}
        """


###############################################################################
### Clean top hits
###############################################################################

rule clean_top_hits:
    input:
        lst = "check_fussion/output/{SAMPLE}/genes_no_fussioned.lst",
        tab = "trinity/trinity_{SAMPLE}/Trinity_refTrans.blastn.hist_bin90.list.names.keep.tab"

    output:
        tab = "minimap2/{SAMPLE}/Trinity_refTrans.blastn.hist_bin90.chimerics.keep.tab"

    log:
        LOG_DIR / "subworkflow4/clean_top_hits/{SAMPLE}.log"

    shell:
        """
        (grep -w -f {input.lst} {input.tab} > {output.tab}) &> {log}
        """


###############################################################################
### Execute minimap script
###############################################################################

rule execute_minimap:
    input:
        "minimap2/{STRAIN}_{TISSUE}/Trinity_refTrans.blastn.hist_bin90.chimerics.keep.tab"
        
    
    output:
        "minimap2/{STRAIN}_{TISSUE}/Trinity_refTrans.blastn.hist_bin90.chimerics.keep.tab"

    log:
        LOG_DIR /  "subworkflow4/execute_minimap/{SAMPLE}.log"

    shell:
        """
        (./scripts/minimap.sh) &> {log}
        """


###############################################################################
### Obtain gene status 
###############################################################################

rule obtain_gene_status:
    input:
        file = "input_files/fbgn_fbtr_fbpp_expanded_fb_2019_06.tsv"

    output:
        tab = "minimap2/genes_status.tab"

    log:
        LOG_DIR / "subworkflow4/obtain_gene_status.log"

    shell:
        """
        (grep "^Dmel" {input.file} | cut -f 3,7 | sort -u > {output.tab}) &> {log}
        """


###############################################################################
### Execute clean Minimap2
###############################################################################

rule execute_clean_minimap:
    input:
        "minimap2/genes_status.tab"

    output:
        "/minimap2/{STRAIN}_{TISSUE}/transcripts_status.lst"

    log:
        LOG_DIR / "subworkflow4/execute_clean_minimap/{SAMPLE}.log"

    shell:
        """
        (./clean_minimap.sh) &> {log}
        """


###############################################################################
### Extract chimerics that are not fussioned
###############################################################################

rule obtain_transcript_list:
    input:
        gtf = config["gene_ano"]
    
    output:
        lst = "list_of_transcripts.txt"

    log:
        LOG_DIR / "subworkflow4/obtain_transcript_list.log"

    shell:
        """
        (cut -f9 {input.gtf} \
        | cut -f6 -d' ' \
        | tr -d '"' \
        | tr -d ";" \
        | sort -u > {output.lst} \
        ) &> {log}
        """

###############################################################################
### Calculate number of exons
###############################################################################

rule calculate_num_exons:
    input:
        "list_of_transcripts.txt"
    
    output:
        "all.exons.gtf"
    
    log:
        LOG_DIR / "subworkflow4/calculate_num_exons.log"

    shell:
        """
        (for i in `less {input}`
        do
	        nExon=$(grep -e "transcript_id \"$i\";" input_files/dmel.gtf | grep CDS |wc -l)
	        if [ "$nExon" -gt 3 ]
		        then grep -e "transcript_id \"$i\";" input_files/dmel.gtf | grep CDS |  sed '1d;$d'
        fi
        done  > {output} \
        ) &> {log}
        """

###############################################################################
### Sort the exons
###############################################################################

rule store_unique_exons:
    input:
        "all.exons.gtf"
    output:
        "all_exons_unique.gtf"

    log:
        LOG_DIR / "subworkflow4/store_unique_exons.log"

    shell:
        """
        (cat {input} | cut -f1,2 -d' ' | sort -u > {output}) &> {log}
        """



###############################################################################
### Store random unique values
###############################################################################

rule store_random_unique:
    input:
        "all_exons_unique.gtf"

    output:
        "dmel_500_random_exons3.gtf"

    log:
        LOG_DIR / "subworkflow4/store_random_unique.log"

    shell:
        """
        (cat {input} | shuf -n 500 > {output}) &> {log}
        """

###############################################################################
### Iterate for information
###############################################################################

rule iterate_for_info:
    input:
        "dmel_500_random_exons3.gtf"
    
    output:
        "AGsites.bed",
        "GTsites.bed",
        "dmel_500_random_exons3.gtf"
    
    log:
        LOG_DIR / "subworkflow4/iterate_for_info.log"

    shell:
        """
        (while IFS= read -r exonLine
        do

            chr=$(echo $exonLine | cut -f1 -d' ' )
            start=$(echo $exonLine | cut -f4 -d' ' )
            end=$(echo $exonLine | cut -f5 -d' ' )
            strand=$(echo $exonLine | cut -f7 -d' ' )
            info=$(echo $exonLine | cut -f10 -d' ' | tr -d '"' |tr -d ';')

            if [ "$strand" == "+" ]
            then
                AGstartPlus=$(($start-11))
                AGendPlus=$(($start+1))
                GTstartPlus=$(($end-4))
                GTendPlus=$(($end+7))
                echo -e "$chr\t$AGstartPlus\t$AGendPlus\t$info\t.\t$strand" >> AGsites.bed
                echo -e "$chr\t$GTstartPlus\t$GTendPlus\t$info\t.\t$strand" >> GTsites.bed
            fi

            if [ "$strand" == "-" ]
            then
                GTstartPlus=$(($start-8))
                GTendPlus=$(($start+3))
                AGstartPlus=$(($end-2))
                AGendPlus=$(($end+10))
                echo -e "$chr\t$AGstartPlus\t$AGendPlus\t$info\t.\t$strand" >> AGsites.bed
                echo -e "$chr\t$GTstartPlus\t$GTendPlus\t$info\t.\t$strand" >> GTsites.bed
            fi

        done < dmel_500_random_exons3.gtf) &> {log}
        """ 


###############################################################################
### Extract sequences for MEME to run
###############################################################################

rule extract_sequences_4_meme:
    input:
        AG = "AGsites.bed",
        GT = "GTsites.bed",
        data = "input_files/dmel-chr.fasta"

    output:
        AG = "AGsites.fasta",
        GT = "GTsites.fasta"

    log:
        LOG_DIR / "subworkflow4/extract_sequences_4_meme"

    shell:
        """
        (bedtools getfasta -fi {input.data}-bed {input.AG} -s -name > {output.AG}
        bedtools getfasta -fi {input.data} -bed {input.GT} -s -name > {output.GT}) &> {log}
        """


###############################################################################
### Execute MEME program
###############################################################################

rule run_meme:
    input:
        AG = "AGsites.fasta",
        GT = "GTsites.fasta"
    output:
        "AGmotif.txt",
        "GTmotif.txt"
    log:
        LOG_DIR / "subworkflow4/run_meme.log"

    shell:
        """
        (meme {input.AG} -dna -oc \
        AGmotif -nostatus \
        -time 18000 -mod zoops \
        -nmotifs 1 -minw 6 \
        -maxw 12 -objfun classic \
        -revcomp -markov_order 0 \
        
        meme {input.GT} -dna -oc \
        GTmotif -nostatus \
        -time 18000 -mod zoops \
        -nmotifs 1 -minw 6 \
        -maxw 12 -objfun classic \
        -revcomp -markov_order 0 \
        ) &> {log}
        """
