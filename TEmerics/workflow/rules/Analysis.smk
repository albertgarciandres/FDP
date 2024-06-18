###############################################################################
### WHAT IT DOES?
###############################################################################
#   - Perform coding profile assesment
#   - Perform expression analysis
#   - Perform pfam domains analysis


from pathlib import Path
OUT_DIR_ANNO = Path(config["out_dir_anno"])
OUT_DIR_IND = Path(config["out_dir_index"])
OUT_DIR_MAP = Path(config["out_dir_map"])
OUT_DIR_IMP = Path(config["out_dir_improved"])
LOG_DIR = Path(config["local_log"])
include: "common.smk"

###############################################################################
### Global configuration
###############################################################################


localrules:
    finish_analysis,


###############################################################################
### Local rules
###############################################################################


rule finish_analysis:
    input:
        pfam = OUT_DIR / "analysis" / "{SAMPLE}" / "pfam.output",
        table = OUT_DIR / "analysis" / "Chimerics.tab",



###############################################################################
### Execute cpat analysis
###############################################################################

rule coding_profiles:
    input:
        fasta = "minimap2/{SAMPLE}/.fa",
        hexa = DATA / "fly_Hexamer.tsv",
        logit = DATA / "Fly_logitModel.RData",
    
    output:
        cpat = "analysis/{SAMPLE}",

    log:
        LOG_DIR / "subworkflow5/coding_profiles/{SAMPLE}.log"
    
    conda:
        ENV_DIR / "cpat.yalm"

    shell:
        """
        (cpat.py -x {input.hexa} -d {input.logit}\
			--top-orf=100 \
			--antisense \
			-g {input.fasta} \
			-o {output.cpat}\
        ) &> {log}
        """

###############################################################################
### Calculate expression
###############################################################################

rule expression_analysis:
    input:
        file = "trinity" / "{SAMPLE}" / "salmon_outdir " / "quant.sf"

    output:
        expr = "results" / "exp.tab"
    log:
        LOG_DIR / "subworkflow4/chimeric_set.log"

    shell:
        """
        (grep -w "$transcript_trinity" {input.file} | cut -f 4 >> {output.expr}) &> {log}
        """


###############################################################################
### Generate chimeric table
###############################################################################

rule chimeric_set:
    input:
        file = "resluts" / "{SAMPLE}" / "transcripts.lst"

    output:
        tab = "analysis" / "Chimerics.tab"

    log:
        LOG_DIR / "subworkflow4/chimeric_set.log"

    shell:
        """
        (echo -e "$condition\t$transcript_trinity\t$transcript_stringtie\t$class\t$transcript\t$gene\t$lengthTranscript\t$totalExons\t$posExon\t$exon\t$coordExon\t$coordTE\tTE_overlap_3\t$TE_consensus\t$TE_family\t$TE_superfamily\t$TE_order\t$TE_class\t$TE_length_incorporated\t$TE_length_total\t$CP\t$geneStatus\t$expr\t$SS\t$result\t$group\t$roo_type" >> chimerics.tab) &> {log}
        """



###############################################################################
### Perform pfam database analysis
###############################################################################

rule pfam_scan_analysis:
    input:
        script = SCRIPTS_DIR / "pfam_scan.sh",
        
    output:
        "analysis" / "{SAMPLE}" / "pfam.output"

    log:
        LOG_DIR /  "subworkflow4/execute_minimap/{SAMPLE}.log"

    shell:
        """
        ({input.script}) &> {log}
        """




