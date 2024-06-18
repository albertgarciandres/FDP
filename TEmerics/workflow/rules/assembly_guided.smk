###############################################################################
### WHAT IT DOES?
###############################################################################
#
#   This subworkflow aims to create a local database to be used in the de novo
#   assembly filtering with blastn.
#   
#   Steps:
#       - Preprocess data with fastp
#       - Map input reads to the reference genome with HISAT2
#       - Build transcriptome assembly with StringTie
#       - Create BLAST+ local database
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
STATS_DIR = Path(config["stats_dir"])

LOCAL_LOG = Path(config["local_log"])

###############################################################################
### Including functions
###############################################################################


include: "common.smk"


###############################################################################
### Reading samples table
###############################################################################

samples_table = pd.read_csv(
    config["samples"],
    header=0,
    index_col=0,
    comment="#",
    engine="python",
    sep="\t",
)

###############################################################################
### Obtain wildcards
###############################################################################

SAMPLE = (pd.unique(samples_table.index.values),)

###############################################################################
### Global configuration
###############################################################################


localrules:
    finish_guided,


###############################################################################
### Local rules
###############################################################################


rule finish_guided:
    input:
        sample_1=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_1.fq.gz",
        sample_2=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_2.fq.gz",
        hisat2_index=expand(
            OUT_DIR / "guided_assembly" / "hisat2_ref_index.{num}.ht2",
            num=[1, 2, 3, 4, 5, 6, 7, 8],
        ),
        maps=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_hisat2.bam",
        gtf=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_struct_defs.gtf",
        tmap_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.tmap",
        db_name=OUT_DIR / "guided_assembly" / "all_seqs_database.fasta",


###############################################################################
### Filter libraries with fastp
###############################################################################

if config["get_stats"] == 1:

    rule filter_fastq:
        input:
            sample_1=lambda wildcards: expand(
                pd.Series(samples_table.loc[wildcards.SAMPLE, "sample_file_1"]).values,
            ),
            sample_2=lambda wildcards: expand(
                pd.Series(samples_table.loc[wildcards.SAMPLE, "sample_file_2"]).values,
            ),
        output:
            sample_1=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_1.fq.gz",
            sample_2=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_2.fq.gz",
            html_report=STATS_DIR / "{SAMPLE}" / "fastp_report.html",
        params:
            min_read=config["min_read_len"],
            min_qual=config["min_qual_val_per_base"],
        resources:
            threads=2,
        log:
            LOG_DIR / "filter_fastq_{SAMPLE}.log",
        conda:
            ENV_DIR / "fastp.yalm"
        shell:
            """
            (fastp \
            --in1 {input.sample_1} \
            --in2 {input.sample_2} \
            -l {params.min_read} \
            -q {params.min_qual} \
            -w {resources.threads} \
            --out1 {output.sample_1} \
            --out2 {output.sample_2} \
            -h {output.html_report} \
            ) &> {log}
            """

else:

    rule filter_fastq:
        input:
            sample_1=lambda wildcards: expand(
                pd.Series(samples_table.loc[wildcards.SAMPLE, "sample_file_1"]).values,
            ),
            sample_2=lambda wildcards: expand(
                pd.Series(samples_table.loc[wildcards.SAMPLE, "sample_file_2"]).values,
            ),
        output:
            sample_1=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_1.fq.gz",
            sample_2=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_2.fq.gz",
        params:
            min_read=config["min_read_len"],
            min_qual=config["min_qual_val_per_base"],
        resources:
            threads=2,
        log:
            LOG_DIR / "filter_fastq_{SAMPLE}.log",
        conda:
            ENV_DIR / "fastp.yalm"
        shell:
            """
            (fastp \
            --in1 {input.sample_1} \
            --in2 {input.sample_2} \
            -l {params.min_read} \
            -q {params.min_qual} \
            -w {resources.threads} \
            --out1 {output.sample_1} \
            --out2 {output.sample_2} \
            ) &> {log}
            """


############################# GUIDED ASSEMBLY #################################
#
# HISAT2 docu: http://daehwankimlab.github.io/hisat2/manual/
#
###############################################################################
### Get exons with HISAT2
###############################################################################


rule get_exons_hisat2:
    input:
        ref_anno=config["ref_anno"],
    output:
        exons=OUT_DIR / "guided_assembly" / "hisat2_exons.exon",
    log:
        LOG_DIR / "get_exons_hisat2.log",
    conda:
        ENV_DIR / "hisat2.yalm"
    shell:
        """
        (hisat2_extract_exons.py \
        {input.ref_ano} > {output.exons} \
        ) &> {log}
        """


###############################################################################
### Get splice sites with HISAT2
###############################################################################


rule get_splice_sites_hisat2:
    input:
        ref_anno=config["ref_anno"],
    output:
        splice_sites=OUT_DIR / "guided_assembly" / "hisat2_splice_sites.ss",
    log:
        LOG_DIR / "get_splice_sites_hisat2.log",
    conda:
        ENV_DIR / "hisat2.yalm"
    shell:
        """
        (hisat2_extract_splice_sites.py \
        {input.ref_ano} > {output.splice_sites} \
        ) &> {log}
        """


################################################################################
### Generate HISAT2 reference index
###############################################################################


rule generate_hisat2_ref_index:
    input:
        ref_genome=config["ref_genome"],
        splice_sites=OUT_DIR / "guided_assembly" / "hisat2_splice_sites.ss",
        exons=OUT_DIR / "guided_assembly" / "hisat2_exons.exon",
    output:
        hisat2_index=expand(
            OUT_DIR / "guided_assembly" / "hisat2_ref_index.{num}.ht2",
            num=[1, 2, 3, 4, 5, 6, 7, 8],
        ),
    params:
        index_basename=OUT_DIR / "guided_assembly" / "hisat2_ref_index",
    resources:
        threads=4,
    log:
        LOG_DIR / "generate_hisat2_ref_index.log",
    conda:
        ENV_DIR / "hisat2.yalm"
    shell:
        """
        (hisat2-build {input.ref_genome} \
        {params.index_basename} \
        --ss {input.splice_sites} \
        --exon {input.exons} \
        -p {resources.threads} \
        -q \
        ) &> {log}
        """


###############################################################################
###  Map reads with HISAT2 index (output for StringTie --dta)
###############################################################################


rule map_reads_hisat2_index:
    input:
        sample_1=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_1.fq.gz",
        sample_2=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "reads_trimmed_2.fq.gz",
    output:
        maps=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_hisat2.sam",
    params:
        index_basename=OUT_DIR / "guided_assembly" / "hisat2_ref_index",
    resources:
        threads=4,
    log:
        LOG_DIR / "map_{SAMPLE}_hisat2.log",
    conda:
        ENV_DIR / "hisat2.yalm"
    shell:
        """
        (hisat2 \
        -1 {input.sample_1} \
        -2 {input.sample_2} \
        -x {params.index_basename} \
        -S {output.maps} \
        --dta \
        -p {resources.threads} \
        -q \
        ) &> {log}
        """


###############################################################################
### Convert SAM to BAM
###############################################################################


rule convert_sam_to_bam:
    input:
        maps=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_hisat2.sam",
    output:
        maps=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_hisat2.bam",
    resources:
        threads=4,
    log:
        LOG_DIR / "convert_{SAMPLE}_sam_to_bam.log",
    conda:
        ENV_DIR / "samtools.yalm"
    shell:
        """
        (samtools \
        -sort \
        -@ {resources.threads} \
        -o {output.maps} \
        {input.maps} \
        ) &> {log}
        """


###############################################################################
### Get structural definitions of the assembled transcripts
###############################################################################


rule get_aln_transcripts_gtf:
    input:
        ref_anno=config["ref_anno"],
        maps=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_hisat2.bam",
    output:
        gtf=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_struct_defs.gtf",
    params:
        min_gap=config["min_gap"],
        min_bp=config["min_bp_around"],
        max_mm=config["max_multimappers"],
        min_read_cov=config["min_read_cov"],
        min_iso=config["min_iso_abundance"],
        min_spliced=config["min_spliced_reads"],
    resources:
        threads=4,
    log:
        LOG_DIR / "get_{SAMPLE}_aln_transcripts_gtf.log",
    conda:
        ENV_DIR / "stringtie.yalm"
    shell:
        """
        (stringtie {input.bam} \
        -G {input.ref_anno} \
        -a {params.min_bp} \
        -M {params.max_mm} \
        -g {params.min_gap} \
        -f {params.min_iso} \
        -j {params.min_spliced} \
        -c {params.min_read_cov} \
        -p {resources.threads} \
        -o {output.gtf} \
        ) &> {log}
        """


###############################################################################
### Merge all structural definition files
###############################################################################


rule merge_stringtie_gtf:
    input:
        ref_anno=config["ref_anno"],
        gtf=OUT_DIR / "guided_assembly" / "{SAMPLE}" / "alignments_struct_defs.gtf",
        files=expand(
            OUT_DIR / "guided_assembly" / "{smp}" / "alignments_struct_defs.gtf",
            smp=SAMPLE,
        ),
    output:
        merged_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.gtf",
    params:
        min_gap=config["min_gap"],
        min_tpm=config["min_tpm"],
        min_fpkm=config["min_fpkm"],
        min_iso=config["min_iso_abundance"],
        min_read_cov=config["min_read_cov"],
    resources:
        threads=4,
    log:
        LOG_DIR / "merge_stringtie_gtf.log",
    conda:
        ENV_DIR / "stringtie.yalm"
    shell:
        """
        (stringtie \
        --merge {input.files} \
        -T {params.min_tpm} \
        -g {params.min_gap} \
        -f {params.min_iso} \
        -G {input.ref_anno} \
        -F {params.min_fpkm} \
        -p {resources.threads} \
        -o {output.merged_gtf} \
        -c {params.min_read_cov} \
        ) &> {log}
        """


###############################################################################
### Compare reference genome with gtf
# More information in: http://ccb.jhu.edu/software/stringtie/gffcompare.shtml
###############################################################################
### Check -G parameter
# Should not matter (gotta test)
###############################################################################

if config["get_stats"] == 1:

    rule compare_merged_gtf_accuracy:
        input:
            ref_anno=config["ref_anno"],
            merged_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.gtf",
        output:
            tmap_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.tmap",
        params:
            prefix=STATS_DIR / "guided_assembly" / "gffcmp_merged",
            tmap_def=OUT_DIR / "guided_assembly" / "stringtie_merged.gtf.tmap",
        log:
            LOG_DIR / "compare_merged_gtf_accuracy.log",
        conda:
            ENV_DIR / "gffcompare.yalm"
        shell:
            """
            (gffcompare {input.merged_gtf} \
            -r {input.ref_anno} \
            -G \
            -o {params.prefix} && mv {params.tmap_def] {output.tmap_gtf} \
            ) &> {log}
            """
else:
    rule compare_merged_gtf_accuracy:
        input:
            ref_anno=config["ref_anno"],
            merged_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.gtf",
        output:
            tmap_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.tmap",
        params:
            prefix=OUT_DIR / "tmp" / "gffcmp_merged",
            tmap_def=OUT_DIR / "guided_assembly" / "stringtie_merged.gtf.tmap",
        log:
            LOG_DIR / "compare_merged_gtf_accuracy.log",
        conda:
            ENV_DIR / "gffcompare.yalm"
        shell:
            """
            (gffcompare {input.merged_gtf} \
            -r {input.ref_anno} \
            -G \
            -o {params.prefix} && mv {params.tmap_def] {output.tmap_gtf} && \
            rm -rf {params.prefix} \
            ) &> {log}
            """

###############################################################################
### Retrieve FASTA sequence from merged GTF
###############################################################################


rule get_transcripts_seqs:
    input:
        merged_gtf=OUT_DIR / "guided_assembly" / "stringtie_merged.gtf",
        ref_genome=config["ref_genome"],
    output:
        fasta=OUT_DIR / "guided_assembly" / "stringtie_merged.fasta",
    log:
        LOG_DIR / "get_transcript_seqs.log",
    conda:
        ENV_DIR / "gffread.yalm"
    shell:
        """
        (gffread {input.merged_gtf} \
        -g {input.ref_genome} \
        -w {output.fasta} \
        -F \
        ) &> {log}
        """


###############################################################################
### Create BLAST database
# NOTE!!!
# I added the extra param "db_name" so you can choose the db NAME, remove it
# if you want and change the `-db` CLI option in the blast rule(s).
###############################################################################


rule create_blast_local_database:
    input:
        fasta=OUT_DIR / "guided_assembly" / "stringtie_merged.fasta",
    output:
        db_name=OUT_DIR / "guided_assembly" / "all_seqs_database.fasta",
    log:
        LOG_DIR / "create_blast_local_db.log",
    conda:
        ENV_DIR / "blast.yalm"
    shell:
        """
        (makeblastdb \
        -in {input.fasta} \
        -out {output.db_name} \
        -dbtype nucl \
        ) &> {log}
        """
