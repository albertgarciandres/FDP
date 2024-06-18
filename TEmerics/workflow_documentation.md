# TEmerics: workflow Documentation

This document describes all the individual steps of the workflow (rules) along with all the thir-party software used.

## Table of Contents

- [Third-party software used](#third-party-software-used)
- [Description of workflow steps](#description-of-workflow-steps)
  - [Preparatory](#preparatory)
    - [Read sample table](#read-sample-table)
  - [Snakefile](#snakefile)
    - [`finish`](#finish)
  - [Guided assembly workflow](#guided-assembly-workflow)
    - [`finish_guided`](#finish_guided)
    - [`filter_fastq`](#filter_fastq)
    - [`get_exons_hisat2`](#get_exons_hisat2)
    - [`get_splice_sites_hisat2`](#get_splice_sites_hisat2)
    - [`generate_hisat2_ref_index`](#generate_hisat2_ref_index)
    - [`map_reads_hisat2_index`](#map_reads_hisat2_index)
    - [`convert_sam_to_bam`](#convert_sam_to_bam)
    - [`get_aln_transcripts_gtf`](#get_aln_transcripts_gtf)
    - [`merge_stringtie_gtf`](#merge_stringtie_gtf)
    - [`compare_merged_gtf_accuracy`](#compare_merged_gtf_accuracy)
    - [`get_transcripts_seqs`](#get_transcripts_seqs)
    - [`create_blast_local_database`](#create_blast_local_database)
  - [_De novo_ assembly workflow](#de-novo-assembly-workflow)
    - [`finish_deNovo`](#finish_denovo)
    - [`convert_reads_to_fasta`](#convert_reads_to_fasta)
    - [`obtain_trinity_stats`](#obtain_trinity_stats)
    - [`run_blastn_denovo_local_db`](#run_blastn_denovo_local_db)
    - [`evaluate_blastn_results`](#evaluate_blastn_results)
    - [`filter_by_coverage_and_relation`](#filter_by_coverage_and_relation)
    - [`filter_fasta_with_blastn_results`](#filter_fasta_with_blastn_results)
  - [Chimeric identification workflow](#chimeric-identification-workflow)
    - [`finish_chimeric_id`](#finish-chimeric-id)
    - [`annotate_te_repeatmasker`](#annotate-te-repeatmasker)
    - [`find_tandem_repeats_per_transcript`](#find-tandem-repeats-per-transcipt)
    - [`get_trf_table_per_transcript`](#get-trf-table-per-transcript)
    - [`merge_trf_results`](#merge-trf-results)
    - [`get_chimeric_entries`](#get-chimeric-entries)
    - [`get_chimeric_fasta`](#get-chimeric-fasta)
    - [`run_blast_chimeric`](#run-blast-chimeric)
    - [`generate_percentage_table`](#generate-percentage-table)
    - [`filter_percent_table`](#filter-percent-table)
    - [`create_gene_gff`](#create-gene-gff)
    - [`convert_gene_gff_to_bed6`](#convert-gene-gff-to-bed6)
    - [`get_distance_to_closest_gene`](#get-distance-to-closest-gene)
  - [Analysis workflow](#analysis-workflow)
    - [`finish_analysis`](#finish-analysis)
    - [`expression_analysis`](#expression-analysis)
    - [`coding_profiles`](#coding-profiles)
    - [`pfam_scan_analyis`](#pfam-scan-analyis)
    - [`chimeric_set`](#chimeric-set)
- [Output file](#output-file)
    - [Main output](#main-output)
    - [Secondary outputs](#secondary-outputs)

## Third-party software used

> Tag lines were taken from the developers' websites

| Name | License | Tag line | More information |
|:---:|:---:|:---:|:---:|
| fastp | [MIT][license-mit] | _"provide fast all-in-one preprocessing for FastQ files"_ | [code][code-fastp] / [manual][manual-fastp] / [publication][publication-fastp] |
| HISAT2 | [GPLv3][license-gpl3] | _"HISAT2 is a fast and sensitive alignment program for mapping next-generation sequencing reads to a population"_ | [code][code-hisat2] / [manual][manual-hisat2] / [publication][publication-hisat2] |
| StringTie | [MIT][license-mit] | _"efficient transcript assembly and quantitation of RNA-Seq data"_ | [code][code-stringtie] / [manual][manual-stringtie] / [publication][publication-stringtie] |
| SAMtools | [MIT][license-mit] | _"mpileup and other tools for handling SAM, BAM, CRAM"_ | [code][code-samtools] / [manual][manual-samtools] / [publication][publication-samtools] |
| gffread | [MIT][license-mit] | _"GFF/GTF utility providing format conversions, filtering, FASTA sequence extraction and more"_ | [code][code-gffread] / [manual][manual-gffread] / [publication][publication-gffread] |
| gffcompare | [MIT][license-mit] | _"compare and evaluate the accuracy of RNA-Seq transcript assemblers"_ | [code][code-gffcompare] / [manual][manual-gffcompare] / [publication][publication-gffcompare] |
| BLAST+ | [U.S. Public Domain][license-us] | _"builds a BLAST database"_ | [code][code-blast] / [manual][manual-blast] / [publication][publication-blast] |
| RepeatMasker | [OSLv2.1][license-osl2] | _"screens DNA sequences for interspersed repeats and low complexity DNA sequences"_ | [code][code-repeatmasker] / [manual][manual-repeatmasker] / [publication][publication-repeatmasker] |
| BEDtools | [MIT][license-mit] | _"intersect, merge, count, complement, and shuffle genomic intervals from multiple files in widely-used genomic file formats such as BAM, BED, GFF/GTF, VCF"_ | [code][code-bedtools] / [manual][manual-bedtools] / [publications][publication-bedtools] |
| Minimap2 | [MIT][license-mit] | _"sequence alignment program that aligns DNA or mRNA sequences against a large reference database"_ | [code][code-minimap] / [manual][manual-minimap] / [publication][publication-minimap] | 
| Seqtk | [MIT][license-mit] | _"fast and lightweight tool for processing sequences in the FASTA or FASTQ format"_ | [code][code-seqtk] / [manual][manual-seqtk] / [publication][publication-seqtk] |
| CPAT | [GPLv2][license-gpl2] | _"predict RNA’s coding probability based on the RNA sequence"_ | [code][code-cpat] / [manual][manual-cpat] / [publication][publication-cpat] |
| PfamScan | [GPLv3][license-gpl3] | _"identify protein domains in one or more protein sequences"_ | [code][code-pfam] / [manual][manual-pfam] / [publication][publication-pfam] |
| FIMO | [MIT][license-mit] | _"scans a set of sequences for individual matches to each of the motifs you provide"_ | [code][code-fimo] / [manual][manual-fimo] / [publication][publication-fimo] |
| Trinity | [BSDv3][license-bsd3] | _"assembles transcript sequences from Illumina RNA-Seq data"_ | [code][code-trinity] / [manual][manual-trinity] / [publication][publication-trinity] |
| BEDmap | [GPLv2][license-gpl2] | _"retrieve and process signal or other features over regions of interest in BED files"_ | [code][code-bedmap] / [manual][manual-bedmap] / [publication][publication-bedmap] |
| Tandem Repeat FInder | [AGPLv3][license-agpl3] | _"locate and display tandem repeats in DNA sequences"_ | [code][code-trf] / [manual][manual-trf] / [publication][publication-trf] |


## Description of workflow steps

> The workflow consists of five Snakefile files. A main `Snakefile` and an
individual Snakemake file for each step in the workflow (guided assembly, 
_de novo_ assembly, chimeric gene-TE transcript identification and analyses).
The main `Snakefile` contains the configuration file validation along with the
inclusion of the sub-workflows. Each individual step of the workflow is 
described briefly along with some examples, and links to the respective 
software manuals are given. Parameters that can be modified by the user (via 
the samples table and the configuration file) are also described.


### Preparatory

> In order to correctly analyze the data we first have to previously prepare 
some files.

#### Read sample table

For a propper functionality of the workflow, the creation of a `.tsv` table 
providing the path to the input data files we want to analyze as well as a 
classification with its corresponding sample_id must be created.

|  PAIRED-END READ 1  |  PAIRED-END READ 2  | SAMPLE_ID |
|:-------------------:|:-------------------:|-----------|
| `CON1_read1.fastq.gz` | `CON1_read2.fastq.gz` |    CON1   |
| `CON2_read1.fastq.gz` | `CON2_read2.fastq.gz` |    CON2   |

#### Configuration file

This file is for the user to fill with its corresponding file paths as well as 
the corresponding parameters for each tool used depending on which type of 
species they are doing the analysis for.

### Snakefile

> This is the main `Snakefile` which contains the configuration file validation
along with the inclusion of all the subworkflows necessary to perform the 
corresponding analysis _TEmerics_ does.

#### `finish`

- **Input**

### Guided assembly workflow

> This module objective is to preprocess the input data with 
[fastp](#third-party-software-used), create a local database with 
[BLAST+](#third-party-software-used) combining the aligner tool 
[HISAT2](#third-party-software-used) - to map the input reads to a reference 
genome - and [StringTie](#third-party-software-used) to build the transcriptome
 assembly from the resulting alignments.

#### `finish_guided`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Paired-end trimmed read 1 (`.fq.gz`); from 
  [map_reads_hisat2_index](#map_reads_hisat2_index)
  - Paired-end trimmed read 2 (`.fq.gz`); from 
  [map_reads_hisat2_index](#map_reads_hisat2_index)
  - Reference index (`.ht2`); from 
  [map_reads_hisat2_index](#map_reads_hisat2_index)
  - Mapped alignments (`.bam`); from 
  [get_aln_transcripts_gtf](#get_aln_transcripts_gtf)
  - Structural definitions (`.gtf`); from 
  [merge_stringtie_gtf](#merge_stringtie_gtf)
  - Local database (`.fasta`); from 
  [create_blast_local_database](#create_blast_local_database)

#### `filter_fastq`

Preprocess (filter and trim) paired-end FASTQ input data files with [fastp](#third-party-software-used).

- **Input**
  - Paired-end read 1 (`fastq.gz`)
  - Paired-end read 2 (`fastq.gz`)
- **Output**
  - Paired-end trimmed read 1 (`.fq.gz`) used in 
  [map_reads_hisat2_index](#map_reads_hisat2_index) && [convert_reads_to_fasta](#convert_reads_to_fasta)
  - Paired-end trimmed read 2 (`.fq.gz`) used in 
  [map_reads_hisat2_index](#map_reads_hisat2_index) && [convert_reads_to_fasta](#convert_reads_to_fasta)
  - Fastp report (`.html`)
    > Optional

#### `get_exons_hisat2`

Extract exon information from the reference annotation file using the script 
_`hisat2_extract_exons.py`_ of [HISAT2](#third-party-software-used).

- **Input**
  - Gene annotations (`.gtf`)
- **Output**
  - Exon information file (`.exon`) used in 
  [genereate_hisat2_ref_index](#generate_hisat2_ref_index)

#### `get_splice_sites_hisat2`

Extract splice sites information from the reference annotation file using the 
script _`hisat2_extract_splice_sites.py`_ of 
[HISAT2](#third-party-software-used).

- **Input**
  - Gene annotations (`.gtf`)
- **Output**
  - Splice site information (`.ss`) used in 
  [genereate_hisat2_ref_index](#generate_hisat2_ref_index)

#### `generate_hisat2_ref_index`

Generation of the reference genome index given the splice sites and exon 
information with [HISAT2](#third-party-software-used).

- **Input**
  - Gene annotations (`.gtf`)
  - Splice site information (`.ss`)
  - Exon information file (`.exon`)
- **Output**
  - Reference index (`.ht2`) used in 
  [map_reads_hisat2_index](#map_reads_hisat2_index)
  
#### `map_reads_hisat2_index`

Mapping of the paired-end reads to the reference genome index with [HISAT2](#third-party-software-used).

- **Input**
  - Paired-end trimmed read 1 (`.fq.gz`)
  - Paired-end trimmed read 2 (`.fq.gz`)
  - Reference index (`.ht2`)
- **Output**
  - Mapped alignments (`.sam`) used in 
  [convert_sam_to_bam](#convert_sam_to_bam)

#### `convert_sam_to_bam`

Convert and sort mapped alingments from SAM to BAM format with 
[SAMtools](#third-party-software-used).

- **Input**
  - Mapped alignments (`.sam`)
- **Output**
  - Mapped alignments (`.bam`) used in 
  [get_aln_transcripts_gtf](#get_aln_transcripts_gtf)

#### `get_aln_transcripts_gtf`

Obtain structural definitions of the assembled transcripts with 
[StringTie](#third-party-software-used).

- **Input**
  - Gene annotations (`.gtf`)
  - Mapped alignments (`.bam`)
- **Output**
  - Structural definitions (`.gtf`) used in 
  [merge_stringtie_gtf](#merge_stringtie_gtf)

#### `merge_stringtie_gtf`

Merge all structural definition files with 
[StringTie](#third-party-software-used).

- **Input**
  - Gene annotations (`.gtf`)
  - Structural definitions (`.gtf`)
- **Output**
  - Merged structural definitions (`.gtf`) used in 
  [compare_merged_gtf_accuracy](#compare_merged_gtf_accuracy) && 
  [get_transcripts_seqs](#get_transcripts_seqs)

#### `compare_merged_gtf_accuracy`

Compare the merged definitions against a reference genome annotation with 
[gffcompare](#third-party-software-used).

- **Input**
  - Gene annotations (`.gtf`)
  - Merged structural definitions (`.gtf`)
- **Output**
  - Merged comparison (`.tmap`) used in [filter_by_coverage_and_relation](#filter_by_coverage_and_relation)

#### `get_transcripts_seqs`

Retrieve FASTA sequences with the reference genome from the merged transcript 
annotations file with [gffread](#third-party-software-used).

- **Input**
  - Reference genome (`.fasta`)
  - Merged structural definitions (`.gtf`)
- **Output**
  - Merged sequences (`.fasta`) used in 
  [create_blast_local_database](#create_blast_local_database)

#### `create_blast_local_database`

Create a local BLAST database [BLAST+](#third-party-software-used).

- **Input**
  - Merged sequences (`.fasta`)
- **Output**
  - Local database (`.fasta`) used in 
  [run_blastn_denovo_local_db](#run_blastn_denovo_local_db) && 
  [evaluate_blastn_results](#evaluate_blastn_results)

### _De novo_ assembly workflow

> This module goal is to generate some filtering steps from the _de novo_ 
assembly using [Trinity](#third-party-software-used) and 
[BLAST+](#third-party-software-used).

#### `finish_deNovo`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**
  - Trimmed reads (`.fasta`); from 
  [convert_reads_to_fasta](#convert_reads_to_fasta)
  - Blast evalue list (`.list`); from 
  [evaluate_blastn_results](#evaluate_blastn_results)
  - Filtered blast results (`.tsv`); from 
  [filter_by_coverage_and_relation](#filter_by_coverage_and_relation)
  - Filtered reads (`.fasta`); from 
  [filter_fasta_with_blastn_results](#filter_fasta_with_blastn_results)

#### `convert_reads_to_fasta`

Convert the trimmed paired-end reads from FASTQ to FASTA format and perform
_de novo_ transcriptome assembly with [Trinity](#third-party-software-used).

- **Input**
  - Paired-end trimmed read 1 (`.fq.gz`)
  - Paired-end trimmed read 2 (`.fq.gz`)
- **Output**
  - Trimmed reads (`.fasta`) used in 
  [obtain_trinity_stats](#obtain_trinity_stats) && [run_blastn_denovo_local_db](#run_blastn_denovo_local_db) && 
  [evaluate_blastn_results](#evaluate_blastn_results) && 
  [filter_fasta_with_blastn_results](#filter_fasta_with_blastn_results)

#### `obtain_trinity_stats`

Obtain statistics (N50, total base count, GC content, etc.) out of the 
_de novo_ transcriptome assembly with the _`TrinityStats.pl`_ script from 
[Trinity](#third-party-software-used).

> Optional

- **Input**
  - Trimmed reads (`.fasta`)
- **Output**
  - Trinitity statistics (`.stats`)

#### `run_blastn_denovo_local_db`

Search against the local sequence database using the query sequences generated 
from the _de novo_ assembly with [BLAST+](#third-party-software-used).

- **Input**
  - Local database (`.fasta`)
  - Trimmed reads (`.fasta`)
- **Output**
  - Assembly comparison (`.blastn`) used in 
  [evaluate_blastn_results](#evaluate_blastn_results)

#### `evaluate_blastn_results`

Evaluate the coverage of top hits obtained from the BLAST comparison results 
with the perl script _`analize_blastPlus_topHit_coverage.pl`_ from [Trinity](#third-party-software-used).

- **Input**
  - Local database (`.fasta`)
  - Assembly comparison (`.blastn`)
  - Trimmed reads (`.fasta`)
- **Output**
  - Blast evalue list (`.list`) used in [filter_by_coverage_and_relation](#filter_by_coverage_and_relation)

#### `filter_by_coverage_and_relation`

Keep all the sequences with a higher 80% similatiry from the list with a 
created bash script _`filter_blastn_results.sh`_ .

> Optional to obtain statistics

- **Input**
  - Blast evalue list (`.list`)
  - Merged comparison (`.tmap`)
  - Bash script (`.sh`)
- **Output**
  - Filtered blast results (`.tsv`) used in [filter_fasta_with_blastn_results](#filter_fasta_with_blastn_results)

#### `filter_fasta_with_blastn_results`

Extract subsequences from a larger sequence file based on specified coordinates
 or identifiers with [seqtk](#third-party-software-used).

- **Input**
  - Trimmed reads (`.fasta`)
  - Filtered blast results (`.tsv`)
- **Output**
  - Filtered reads (`.fasta`)

### Chimeric identification workflow

> This module goal is to detect Transposable elements and filter them for later
identification and processing.

#### `finish_chimeric_id`

Target rule as required by [Snakemake][docs-snakemake].

> Local rule

- **Input**


## Output file

> Brief description of the output files TEmerics provides.

### Main output

> The main output file is a `.tab` delimited file with the an overall of 
29 variables.

- **Sample_id**: identifier of the data provided by the user.
- **Transcript**: id of the transcript according to the reference genome data. 
- **Transcript Trinity**: ID according to Trinity.
- **TEfamily**: family of the TE.
- **TE coordinate**: coordinates where the TE was detected in the 
assembled region.
- **Transcript Stringtie**: id of the transcript according to the 
reference-guided assembly. Most genes have a FBtr ID, the ones starting 
with MSTRG are transcripts not found in the reference annotation.
- **Class**: transcript type code according to the reference-guided assembly.
- **Gene**: id of the gene in FBgn format.
- **Length transcript**: length of the assembled transcript in bp.
- **Total number of exons**: number of exons detected in the transcript after 
performing the alignment to the genome with [Minimap2][code-minimap].
- **Expon position**: position in which the TE was detected (number).
- **Exon position description**: position in which the TE is detected: first, 
last or middle exon (corresponds to the 3'/5' UTR, internal exon 
classification) or (TE inside) when is detected within.
- **Type**: exon_within_TE, exon_within_TE_3, exon_within_TE_5, TE_overlap_3, 
TE_overlap_5, TE_within_exon.
- **TE consensus**: name of the TE according to the library.
- **TE superfamily**: superfamily of the TE.
- **TE order**: order of the TE.
- **TE class**: class of the TE.
- **length TE incorporated**: length of the TE fragment detected in the 
chimeric transcipts in bp.
- **Length TE total**: length of the TE insertion in bp.
- **CP**: coding potential according to [CPAT][code-cpat]
- **status**: type of transcript according to reference genome annotations 
(RNA, ncRNA, pseudogene, pre_miRNA, tRNA, snoRNA).
- **avgExpr**: average level of expression (in TMM, averaging the expression 
level of three replicates).
- **SS**: a code specifing the meaning of FIMOenrichment column.
- **FIMOenrichment**: motif and p-value if a motif was found with FIMO.
- **group**: 1: overlap and AS insertions group, 2: internal insertions group.
- **score**: score of the TE according to [RepeatMasker][code-repeatmasker].
- **position of roo**: position in the roo in which there's a match 
(low complexity region).
- **roo type**: match of the roo fragment in the consensus
(low complexity region, or LTR).

### Secondary outputs

> Along with it there are some other outputs such as:

- Trinitity statistics (`.stats`)
- Fastp report (`.html`)
- Filtered blast results (`.tsv`)

[code-bedmap]: <https://github.com/bedops/bedops/blob/master/applications/bed/bedmap/src/Bedmap.cpp>
[code-bedtools]: <https://github.com/arq5x/bedtools2>
[code-blast]: <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source>
[code-cpat]: <https://github.com/liguowang/cpat?tab=readme-ov-file>
[code-fastp]: <https://github.com/OpenGene/fastp>
[code-fimo]: <https://github.com/PriceLab/fimoService>
[code-gffcompare]: <https://github.com/gpertea/gffcompare>
[code-gffread]: <https://github.com/gpertea/gffread>
[code-hisat2]: <https://github.com/DaehwanKimLab/hisat2>
[code-minimap]: <https://github.com/lh3/minimap2?tab=readme-ov-file>
[code-pfam]: <https://github.com/SMRUCC/GCModeller/blob/master/src/interops/scripts/PfamScan/PfamScan/pfam_scan.pl>
[code-repeatmasker]: <https://github.com/rmhubley/RepeatMasker>
[code-samtools]: <https://github.com/samtools/samtools?tab=readme-ov-file>
[code-seqtk]: <https://github.com/lh3/seqtk>
[code-stringtie]: <https://github.com/gpertea/stringtie/tree/master?tab=readme-ov-file>
[code-trf]: <https://github.com/Benson-Genomics-Lab/TRF?tab=readme-ov-file>
[code-trinity]: <https://github.com/trinityrnaseq/trinityrnaseq?tab=BSD-3-Clause-1-ov-file>
[docs-snakemake]: <https://snakemake.readthedocs.io/en/stable/>
[manual-bedmap]: <https://bedops.readthedocs.io/en/latest/content/reference/statistics/bedmap.html>
[manual-bedtools]: <https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html>
[manual-blast]: <https://www.ncbi.nlm.nih.gov/books/NBK279684/>
[manual-cpat]: <https://cpat.readthedocs.io/en/latest/#run-cpat-on-local-computer>
[manual-fastp]: <https://manpages.debian.org/testing/fastp/fastp.1.en.html>
[manual-fimo]: <https://meme-suite.org/meme/doc/fimo.html>
[manual-gffcompare]: <https://ccb.jhu.edu/software/stringtie/gffcompare.shtml>
[manual-gffread]: <https://manpages.ubuntu.com/manpages/focal/man1/gffread.1.html>
[manual-hisat2]: <https://daehwankimlab.github.io/hisat2/manual/>
[manual-minimap]: <https://lh3.github.io/minimap2/minimap2.html>
[manual-pfam]: <https://github.com/SMRUCC/GCModeller/blob/master/src/interops/scripts/PfamScan/PfamScan/README>
[manual-repeatmasker]: <https://www.animalgenome.org/bioinfo/resources/manuals/RepeatMasker.html>
[manual-samtools]: <https://www.htslib.org/doc/samtools.html>
[manual-seqtk]: <https://docs.csc.fi/apps/seqtk/>
[manual-stringtie]: <https://ccb.jhu.edu/software/stringtie/index.shtml?t=manual>
[manual-trf]: <https://tandem.bu.edu/trf/definitions>
[manual-trinity]: <https://github.com/trinityrnaseq/trinityrnaseq/wiki>
[publication-bedmap]: <http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract>
[publication-bedtools]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2832824/>
[publication-blast]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC441573/>
[publication-cpat]: <https://academic.oup.com/nar/article/41/6/e74/2902455>
[publication-fastp]: <https://academic.oup.com/bioinformatics/article/34/17/i884/5093234>
[publication-fimo]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3065696/>
[publication-gffcompare]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7222033/>
[publication-gffread]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7222033/>
[publication-hisat2]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4655817/>
[publication-minimap]: <https://academic.oup.com/bioinformatics/article/34/18/3094/4994778>
[publication-pfam]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3965110/>
[publication-repeatmasker]: <https://www.repeatmasker.org/webrepeatmaskerhelp.html>
[publication-samtools]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2723002/>
[publication-seqtk]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5051824/>
[publication-stringtie]: <https://www.nature.com/articles/nbt.3122>
[publication-trf]: <https://academic.oup.com/nar/article/27/2/573/1061099?login=true>
[publication-trinity]: <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3571712/>
[license-agpl3]: <https://opensource.org/license/agpl-v3>
[license-bsd3]: <https://opensource.org/license/bsd-3-clause>
[license-gpl2]: <https://opensource.org/license/gpl-2-0>
[license-gpl3]: <https://opensource.org/license/gpl-3-0/>
[license-mit]: <https://opensource.org/licenses/MIT>
[license-us]: <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/doc/public/LICENSE>
[license-osl2]: <https://opensource.org/license/osl-2-1>
