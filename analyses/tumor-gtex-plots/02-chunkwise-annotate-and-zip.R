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

# chunkwise conversion to jsonl
chunkSize <- 1000000
index <- 0
con <- file(description = input_file, open = "r")
repeat {
  index <- index + 1
  print(paste('Processing rows:', index * chunkSize))
  if(index == 1){
    dataChunk <- read.table(con, nrows = chunkSize, sep = "\t", header = T)
    cols <- colnames(dataChunk)
  } else {
    dataChunk <- read.table(con, nrows = chunkSize, skip = 0, sep="\t", col.names = cols)
  }
  
  if (nrow(dataChunk) == 0){
    print('Processed all files!')
    break
  }
  
  # convert to JSON
  json_file <- gsub(".tsv", replacement = paste0('_', index, ".json"), input_file)
  jsonlite::write_json(x = dataChunk, path = json_file)
  
  # convert to JSONL
  jsonl_file <- gsub(".json", ".jsonl", json_file)
  jq_cmd <- paste('jq --compact-output \'.[]\'', json_file, '>', jsonl_file)
  system(jq_cmd)
  
  # remove json file
  rm_json_cmd <- paste('rm', json_file)
  system(rm_json_cmd)
}
close(con)

# combine all jsonl files
jsonl_files <- gsub('.tsv', '*.jsonl', input_file)
final_jsonl_file <- gsub('.tsv', '.jsonl', input_file)
concat_jsonl_cmd <- paste('cat', jsonl_files, '>', final_jsonl_file)
system(concat_jsonl_cmd)

# gzip JSONL files
gzip_cmd <- paste('gzip --no-name', final_jsonl_file)
system(gzip_cmd)

# remove all intermediate jsonl files
rm_tmp_jsonl_cmd <- paste('rm', jsonl_files)
system(rm_tmp_jsonl_cmd)