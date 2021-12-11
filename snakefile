import os
from os import path

from snakemake.utils import validate
from snakemake.utils import min_version
min_version("5.3")

# ----------------------------------------------------------------

SNAKEDIR = path.dirname(workflow.snakefile)

configfile: SNAKEDIR + "/config.yml"
validate(config, schema="schema/config_schema.yaml")

workdir: config["workdir"]

WORKDIR = config["workdir"]
sample = config["sample_name"]

report: "report/workflow.rst"

plot_group = ["main", "perSample"] if config["multisample"] == "yes" else ["main"]


# ----------------------------------------------------------------


rule all:
    input:
        # transcripts detected curves
        expand("Post_sqanti_qc/" + sample + "_NFLR_curve_{name}.png", name=plot_group),

        # transcript categories - pre filtering
        expand("Post_sqanti_qc/" + sample + "_preFilt_tc_count_{name}.png", name=plot_group),
        expand("Post_sqanti_qc/" + sample + "_preFilt_tc_NFLR_{name}.png", name=plot_group),

        # transcript categories - post filtering
        expand("Post_sqanti_qc/" + sample + "_postFilt_tc_count_{name}.png", name=plot_group),
        expand("Post_sqanti_qc/" + sample + "_postFilt_tc_NFLR_{name}.png", name=plot_group),

        # transcripts ranked
        expand("Post_sqanti_qc/" + sample + "_postFilt_transcriptsRanked_{name}.png", name=plot_group)

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
        gene_ids = config["gene_ids"],
        script = SNAKEDIR + "/scripts/process_sqanti_results.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {input.sqanti_fasta} {params.prefix} {params.respath} {params.gene_ids}
        """



rule transcripts_detected_curve:
    input:
        sqanti_class = rules.process_sqanti_results.output.res

    output:
        p1 = report(expand("Post_sqanti_qc/" + sample + "_NFLR_curve_{name}.png", name=plot_group), caption = "report/NFLR_curve.rst", category = "Pre-filtering", subcategory = "Transcripts detected vs. expression")


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
        p1 = report(expand("Post_sqanti_qc/" + sample + "_preFilt_tc_count_{name}.png", name=plot_group), caption = "report/tc_count.rst", category = "Pre-filtering", subcategory = "Transcript categories"),
        p3 = report(expand("Post_sqanti_qc/" + sample + "_preFilt_tc_NFLR_{name}.png", name=plot_group), caption = "report/tc_NFLR.rst", category = "Pre-filtering", subcategory = "Transcript categories")


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
        NE = config["NFLR_mean"],
        NFLR = config["NFLR_perSample"],
        min_sample_perc = config["min_sample_perc"],
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
        p1 = report(expand("Post_sqanti_qc/" + sample + "_postFilt_tc_count_{name}.png", name=plot_group), caption = "report/tc_count.rst", category = "Post-filtering", subcategory = "Transcript categories"),
        p3 = report(expand("Post_sqanti_qc/" + sample + "_postFilt_tc_NFLR_{name}.png", name=plot_group), caption = "report/tc_NFLR.rst", category = "Post-filtering", subcategory = "Transcript categories")


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
        p1 = report(expand("Post_sqanti_qc/" + sample + "_postFilt_transcriptsRanked_{name}.png", name=plot_group), caption = "report/transcriptsRanked.rst", category = "Post-filtering", subcategory = "Transcripts ranked")


    params:
        respath = "Post_sqanti_qc",
        prefix = sample + "_postFilt",
        script = SNAKEDIR + "/scripts/plot_transcripts_ranked.R"

    shell:
        """
        /software/R_v4.0.3/bin/Rscript {params.script} {input.sqanti_class} {params.prefix} {params.respath}
        """
