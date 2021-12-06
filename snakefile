import os
from os import path

from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

configfile: "config.yml"
workdir: config["workdir"]

WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]


# ----------------------------------------------------------------


rule all:
    input:
        # transcripts detected curves
        "Post_sqanti_qc/" + sample + "_transcriptsNE_curve.png",
        "Post_sqanti_qc/" + sample + "_transcriptsNFLR_perSample_curve.png",

        # transcript categories - pre filtering
        "Post_sqanti_qc/" + sample + "_preFilt_tc_count_withSd.png",
        "Post_sqanti_qc/" + sample + "_preFilt_tc_count.png",
        "Post_sqanti_qc/" + sample + "_preFilt_tc_NFLR.png",
        "Post_sqanti_qc/" + sample + "_preFilt_tc_NE.png",

        # transcript categories - post filtering
        "Post_sqanti_qc/" + sample + "_postFilt_tc_count_withSd.png",
        "Post_sqanti_qc/" + sample + "_postFilt_tc_count.png",
        "Post_sqanti_qc/" + sample + "_postFilt_tc_NFLR.png",
        "Post_sqanti_qc/" + sample + "_postFilt_tc_NE.png",

        # transcripts ranked
        "Post_sqanti_qc/" + sample + "_postFilt_transcriptsRanked_withSd.png",
        "Post_sqanti_qc/" + sample + "_postFilt_transcriptsRanked.png"

#########################################################################


rule process_sqanti_results:
    input:
        sqanti_class = config["sqanti_class"],
        sqanti_fasta = config["sqanti_fasta"]

    output:
        res = "Post_sqanti_qc/" + sample + "_classification_processed.txt"

    params:
        respath = "Post_sqanti_qc",
        prefix = sample,
        gene_ids = config.get("gene_ids", ""),
        script = SNAKEDIR + "/scripts/process_sqanti_results.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {input.sqanti_fasta} {params.prefix} {params.respath} {params.gene_ids}
        """



rule transcripts_detected_curve:
    input:
        sqanti_class = rules.process_sqanti_results.output.res

    output:
        p1 = "Post_sqanti_qc/" + sample + "_transcriptsNE_curve.png",
        p2 = "Post_sqanti_qc/" + sample + "_transcriptsNFLR_perSample_curve.png"

    params:
        respath = "Post_sqanti_qc",
        prefix = sample,
        script = SNAKEDIR + "/scripts/plot_transcripts_detected_curve.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {params.prefix} {params.respath}
        """



rule transcript_category_preFilt:
    input:
        sqanti_class = rules.process_sqanti_results.output.res

    output:
        p1 = "Post_sqanti_qc/" + sample + "_preFilt_tc_count_withSd.png",
        p2 = "Post_sqanti_qc/" + sample + "_preFilt_tc_count.png",
        p3 = "Post_sqanti_qc/" + sample + "_preFilt_tc_NFLR.png",
        p4 = "Post_sqanti_qc/" + sample + "_preFilt_tc_NE.png"

    params:
        respath = "Post_sqanti_qc",
        prefix = sample + "_preFilt",
        script = SNAKEDIR + "/scripts/plot_transcript_categories.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {params.prefix} {params.respath}
        """



rule filter_sqanti_results:
    input:
        sqanti_class = rules.process_sqanti_results.output.res

    output:
        res = "Post_sqanti_qc/" + sample + "_classification_filtered.txt"

    params:
        NE = config["NE"] if config.get("NE") else 0.3,
        NFLR = config["NFLR"] if config.get("NE") else 0,
        min_sample_perc = config["min_sample_perc"] if config.get("NE") else 0,
        respath = "Post_sqanti_qc",
        prefix = sample,
        script = SNAKEDIR + "/scripts/filter_sqanti_results.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {params.prefix} {params.respath} {params.NE} {params.NFLR} {params.min_sample_perc}
        """


rule transcript_category_postFilt:
    input:
        sqanti_class = rules.filter_sqanti_results.output.res

    output:
        p1 = "Post_sqanti_qc/" + sample + "_postFilt_tc_count_withSd.png",
        p2 = "Post_sqanti_qc/" + sample + "_postFilt_tc_count.png",
        p3 = "Post_sqanti_qc/" + sample + "_postFilt_tc_NFLR.png",
        p4 = "Post_sqanti_qc/" + sample + "_postFilt_tc_NE.png"

    params:
        respath = "Post_sqanti_qc",
        prefix = sample + "_postFilt",
        script = SNAKEDIR + "/scripts/plot_transcript_categories.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {params.prefix} {params.respath}
        """



rule transcripts_ranked:
    input:
        sqanti_class = rules.filter_sqanti_results.output.res

    output:
        p1 = "Post_sqanti_qc/" + sample + "_postFilt_transcriptsRanked_withSd.png",
        p2 = "Post_sqanti_qc/" + sample + "_postFilt_transcriptsRanked.png"

    params:
        respath = "Post_sqanti_qc",
        prefix = sample + "_postFilt",
        script = SNAKEDIR + "/scripts/plot_transcripts_ranked.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {params.prefix} {params.respath}
        """
