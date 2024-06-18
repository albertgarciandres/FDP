###############################################################################
### WHAT IT DOES?
###############################################################################
#
#   This subworkflow performs the de novo assembly as well as the first
#   filtering steps for the identification of chimeric transcritps.
#
#   Steps:
#   - Perform de novo assembly with Trinity
#   - Search against the local database created with BLAST+
#   - Filter results by coverage and 80% similarity 
#   - Generate FASTA file with filtered sequences
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
ENV_DIR = Path(f"{workflow.basedir}/envs")
SCRIPTS_DIR = Path(config["scripts_dir"])
STATS_DIR = Path(Config["stats_dir"])

LOCAL_LOG = Path(config["local_log"])

###############################################################################
### Including functions
###############################################################################


include: "common.smk"


###############################################################################
### Global configuration
###############################################################################


localrules:
    finish_deNovo,


###############################################################################
### Local rules
###############################################################################


rule finish_deNovo:
    input:
        trim_fa=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_trimmed.fasta",
        filt_fa=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_filtered.fasta",
        list=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "blastn_eval.list",
        names=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "filtered_blastn_results.tsv",


###############################################################################
### Convert reads to FASTA
# If not using the "condition file" CLI options --left and --right to specify
# input files.
###############################################################################


rule convert_reads_to_fasta:
    input:
        sample_1=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_1.fq.gz",
        sample_2=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_2.fq.gz",
    output:
        fasta=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_trimmed.fasta",
    params:
        seq="fq",
        out_dir=OUT_DIR / "denovo_assembly" / "{SAMPLE}",
        out_file=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "Trinity.fasta",
    resources:
        max_mem=78,
        cpu=12,
    log:
        LOG_DIR / "trinity_reads_to_fasta_{SAMPLE}.log",
    conda:
        ENV_DIR / "trinity.yalm"
    shell:
        """
        (Trinity \
        --left {input.sample_1} \
        --right {input.sample_2} \
        --seqType {params.seq} \
        --output {params.out_dir} \
        --jaccard_clip \
        --max_memory {resources.max_mm} \
        --CPU {resources.cpu} && mv {params.out_file} {output.fasta} \
        ) &> {log}
        """


###############################################################################
### Obtain Trinity FASTA files stats
# Script:https://github.com/trinityrnaseq/trinityrnaseq/blob/master/util/TrinityStats.pl
###############################################################################

if Config["get_stats"] == 1:

    rule obtain_trinity_stats:
        input:
            file=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_trimmed.fasta",
        output:
            stats=STATS_DIR / "Trinity_{SAMPLE}.stats",
        log:
            LOG_DIR / "obtain_trinity_stats_{SAMPLE}.log",
        conda:
            ENV_DIR / "trinity.yalm"
        shell:
            """
           (TrinityStats.pl {input.file} > {output.stats}) &> {log}
           """


###############################################################################
### Run BLAST: DeNovo assembly against the guided assembly local data base
###############################################################################


rule run_blastn_denovo_local_db:
    input:
        db_name=OUT_DIR / "guided_assembly" / "all_seqs_database.fasta",
        query=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_trimmed.fasta",
    output:
        blast=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "Trinity_refTrans.blastn",
    params:
        e_val=config["e_val_assembly"],
        max_target_seqs=config["max_target_seqs"],
    resources:
        threads=12,
    log:
        LOG_DIR / "run_blastn_{SAMPLE}_denovo_local_db.log",
    conda:
        ENV_DIR / "blast.yalm"
    shell:
        """
        (blastn \
        -query {input.query} \
        -db {input.db} \
        -evalue {params.e_val} \
        -out {output.blast} \
        -outfmt  6 \
        -dust no \
        -task megablast \
        -num_threads {resources.threads} \
        -max_target_seqs {params.max_target_seqs} \
        ) &> {log}
        """


###############################################################################
### Evaluate the quality of the blast results
# Script in: https://github.com/macmanes/trinityrnaseq-1/blob/master/util/analyze_blastPlus_topHit_coverage.pl
###############################################################################


rule evaluate_blastn_results:
    input:
        db_name=OUT_DIR / "guided_assembly" / "all_seqs_database.fasta",
        query=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_trimmed.fasta",
        blast=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "Trinity_refTrans.blastn",
    output:
        list=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "blastn_eval.list",
    params:
        list=OUT_DIR
        / "denovo_assembly"
        / "{SAMPLE}"
        / "Trinity_refTrans.blastn.hist.list",
    log:
        LOG_DIR / "evaluate_blastn_{SAMPLE}_results.log",
    conda:
        ENV_DIR / "blast.yalm"
    shell:
        """
        (analyze_blastPlus_topHit_coverage.pl {input.blast} \
        {input.query} \
        {input.db_name} && mv {params.list} {output.list} \
        ) &> {log}
        """


###############################################################################
### Create file and keep the > 80% results
###############################################################################

if Config["get_stats"] == 1:

    rule filter_by_coverage_and_relation:
        input:
            lst=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "blastn_eval.list",
            tmap_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.tmap",
            script=SCRIPTS_DIR / "filter_blastn_results.sh",
        output:
            names=OUT_DIR
            / "denovo_assembly"
            / "{SAMPLE}"
            / "filtered_blastn_results.tsv",
        params:
            stats_dir=STATS_DIR / "{SAMPLE}",
        log:
            LOG_DIR / "filter_blastn_{SAMPLE}_by_coverage_relation.log",
        shell:
            """
            (bash {input.script} \
            -l {input.lst} \
            -t {input.tmap_gtf} \
            -s {params.stats_dir} \
            -o {output.names} \
            ) &> {log}
            """

else:

    rule filter_by_coverage_and_relation:
        input:
            lst=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "blastn_eval.list",
            tmap_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.tmap",
            script=SCRIPTS_DIR / "filter_blastn_results.sh",
        output:
            names=OUT_DIR
            / "denovo_assembly"
            / "{SAMPLE}"
            / "filtered_blastn_results.tsv",
        log:
            LOG_DIR / "filter_blastn_{SAMPLE}_by_coverage_relation.log",
        shell:
            """
            (bash {input.script} \
            -l {input.lst} \
            -t {input.tmap_gtf} \
            -o {output.names} \
            ) &> {log}
            """


###############################################################################
### Filter Trinity FASTA files with the blastn filtered results
###############################################################################


rule filter_fasta_with_blastn_results:
    input:
        query=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_trimmed.fasta",
        names=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "filtered_blastn_results.tsv",
    output:
        fasta=OUT_DIR / "denovo_assembly" / "{SAMPLE}" / "reads_filtered.fasta",
    log:
        LOG_DIR / "filter_{SAMPLE}_fasta_with_blastn_results.log",
    conda:
        ENV_DIR / "seqtk.yalm"
    shell:
        """
        (seqtk subseq {input.query} {input.names} > {output.fasta}) &> {log}
        """
