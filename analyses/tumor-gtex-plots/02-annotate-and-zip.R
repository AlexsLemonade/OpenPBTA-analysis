# Script to annotate with MONDO, RMTL and EFO fields
# convert to JSONL and gzip

# load libraries
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(jsonlite))

option_list <- list(
  make_option(c("--input_file"), type = "character",
              help = "input file to annotate, convert to jsonl and gzip")
)
opt <- parse_args(OptionParser(option_list = option_list))
input_file <- opt$input_file

# root directory
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
annotator <- file.path(root_dir, 'analyses', 'long-format-table-utils', 'annotator', 'annotator-cli.R')

# run annotator
annotate_cmd <- paste('Rscript --vanilla', annotator, '-r -v -c MONDO,RMTL,EFO',
                      '-i', input_file,
                      '-o', input_file)
system(annotate_cmd)

# convert to JSON
annotated_file <- data.table::fread(input_file)
json_file <- gsub(".tsv", ".json", input_file)
jsonlite::write_json(x = annotated_file, path = json_file)

# convert to JSONL
jsonl_file <- gsub(".json", ".jsonl", json_file)
jq_cmd <- paste('jq --compact-output \'.[]\'', json_file, '>', jsonl_file)
system(jq_cmd)

# remove json file
rm_json_cmd <- paste('rm', json_file)
system(rm_json_cmd)

# gzip JSONL files
gzip_cmd <- paste('gzip --no-name', jsonl_file)
system(gzip_cmd)

