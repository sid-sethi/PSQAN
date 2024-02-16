# Configuration file for PSQAN pipeline

To configure this workflow, modify [config.yml](config.yml) according to your dataset and parameters.

> **_TIP:_** Run PSQAN with default values and use the [NFLR_curve plots](../example_output/ENSG00000TEST/NFLR_curve_main.png) to identify optimal values of `min_exp_perSample`, `min_exp_mean` and `min_sample_perc` for your dataset.

Parameter | Description
------ | -----------
workdir [STRING] | [**REQUIRED**] Path of working/project directory. This directory should exist. All PSQAN results would be saved within this directory
abundance [FILE] | [**REQUIRED**] Transcript abundance file. Abundance output file from `SQANTI` and `TALON` are accepted by PSQAN. For SQANTI, the file ending with "classification.txt" should be provided, while for TALON, read annotation file ending with "read_annot.tsv" should be used
abundance_type [STRING] | File type of transcript abundance input file. Accepted choices: ["SQANTI", "TALON"]; [default = "SQANTI"]
gene_ids [FILE] | [**REQUIRED**] A file containing Gene IDs of genes of interest to run PSQAN. Gene IDs column name should be **`gene_id`**. If the file contains more than one column, it should be a tab-seperated file - other columns in the file will be ignored
percentage_A_downstream_TTS [NUMERIC] [0-100] | Maximum percent of genomic "A"s allowed in the immediate downstream window of the read/isoform. This helps to filter out internal priming artifacts. If percent of genomic "A"s is high (say > 80), the 3' end site of the isoform is probably not reliable. If using `TALON` abundance file, this percentage is divided by 100 to convert into a fraction. [default = 80]
multisample [BOOLEAN] | Do you have transcript abundance of more than one sample/replicate in abundance file? [default = false]
min_exp_perSample [NUMERIC] [0-100] | Minimum value of normalised expression required for a transcript **PER** sample (or replicate). If $multisample=false$, transcripts with normalised expression $>=$ min_exp_perSample would be retained. If $multisample=true$, this parameter is combined with `min_sample_perc` and transcripts which pass the `min_exp_perSample` threshold in at least `min_sample_perc` % of samples are retained. [default = 0.3]
min_exp_mean [NUMERIC] [0-100] | Minimum value of mean normalised expression **ACROSS** all samples required for a transcript. If $multisample=false$, this parameter is ignored. [default = 0]
min_sample_perc [NUMERIC] [0-100] | Minimum % of samples which should pass the `min_exp_perSample` threshold. For instance, min_sample_perc=30 and min_exp_perSample=10 would mean that transcripts identified with a minimum normalised expression of 10 in at least 30% of the total samples would be retained. If $multisample=false$, this parameter is ignored. [default = 0]