#!/bin/bash
# This script runs the sex-prediction-from-RNASeq analysis
# Author's Name Bill Amadio 2019

set -e
set -o pipefail

# This script should always run as if it were being called from
# the directory it lives in.

script_directory="$(perl -e 'use File::Basename;
  use Cwd "abs_path";
  print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# option for only running a single tail percent in continuous integration
PERCENT=${OPENPBTA_PERCENT:-1}


#--------USER-SPECIFIED ARGUMENTS

# output directory of script 01, input directory of scripts 02 and 03
PROCESSED=processed_data

# outtput directory of script 02, input directory of script 03
MODELS=models

# output directory of script 03
RESULTS=results

# argument for script 01 processing and script 02 and 03 input file specification
SEED=36354

# argument for script 01 output file specification and script 02 and 03 input and output file specification
FILENAME_LEAD=kallisto_stranded

# argument for script 02 processing and script 03 input file specification
# if this is not in CI, we can use the full 10 values for the transcript tail %
# else, use 0.25 only
if [ "$PERCENT" -gt "0" ]; then
  TRANSCRIPT_TAIL_PERCENT_ARRAY=(0.25 0.125 0.0625 0.0313 0.0156 0.0078 0.0039 0.0020 0.0010 0.0005)
else
  TRANSCRIPT_TAIL_PERCENT_ARRAY=(0.25)
fi

# argument for script 01 processing
TRAIN_PERCENT=0.7

# argument for script 02 processing
# this specifies what column in the target data frame will be used as labels
# during training
TRAIN_TARGET_COLUMN=reported_gender

# we will evaluate the model's performance on both of these columns
# be sure to include the value of TRAIN_TARGET_COLUMN in this array as well
targetColumns=("reported_gender" "germline_sex_estimate")

#--------END USER-SPECIFIED ARGUMENTS

# output files for script 01, input files for script 02
# the training expression set and the labels
TRAIN_EXPRESSION_FILE_NAME=${FILENAME_LEAD}_${SEED}_train_expression.RDS
TRAIN_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_train_targets.tsv

# output files for script 01, input files for script 03
# the test expression set and the labels
TEST_EXPRESSION_FILE_NAME=${FILENAME_LEAD}_${SEED}_test_expression.RDS
TEST_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_test_targets.tsv

# output file for script 01
FULL_TARGETS_FILE_NAME=${FILENAME_LEAD}_${SEED}_full_targets.tsv

#concatenate targetColumns array to single string

targetColumns_to_pass=${targetColumns[0]}

for (( k = 1; k < ${#targetColumns[*]}; k++ )); do

  targetColumns_to_pass+=",${targetColumns[k]}"

done


# the first step is to split the data into a training and test set
Rscript --vanilla 01-clean_split_data.R \
  --expression ../../data/pbta-gene-expression-kallisto.stranded.rds \
  --metadata ../../data/pbta-histologies.tsv \
  --output_directory $PROCESSED \
  --train_expression_file_name $TRAIN_EXPRESSION_FILE_NAME \
  --test_expression_file_name $TEST_EXPRESSION_FILE_NAME \
  --train_targets_file_name $TRAIN_TARGETS_FILE_NAME  \
  --test_targets_file_name $TEST_TARGETS_FILE_NAME \
  --full_targets_file_name $FULL_TARGETS_FILE_NAME \
  --seed $SEED \
  --train_percent $TRAIN_PERCENT \
  --train_target_column $TRAIN_TARGET_COLUMN \
  --target_columns $targetColumns_to_pass



for (( i = 0; i < ${#TRANSCRIPT_TAIL_PERCENT_ARRAY[*]}; i++ )); do

  # output files for script 02
  # re: the elasticnet model
  MODEL_OBJECT_FILE_NAME=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_model_object.RDS
  MODEL_TRANSCRIPTS_FILE_NAME=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_model_transcripts.RDS
  MODEL_COEFS_FILE_NAME=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}_model_coefs.tsv


  Rscript --vanilla 02-train_elasticnet.R \
    --train_expression_file_name ${PROCESSED}/$TRAIN_EXPRESSION_FILE_NAME \
    --train_targets_file_name ${PROCESSED}/$TRAIN_TARGETS_FILE_NAME \
    --output_directory $MODELS \
    --model_object_file_name $MODEL_OBJECT_FILE_NAME \
    --model_transcripts_file_name $MODEL_TRANSCRIPTS_FILE_NAME \
    --model_coefs_file_name $MODEL_COEFS_FILE_NAME \
    --train_target_column $TRAIN_TARGET_COLUMN \
    --transcript_tail_percent ${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]} \
    --seed $SEED


  # if there is a test set, e.g., the entire dataset was not used for training
  if [ ! $TRAIN_PERCENT == 1 ]; then

    # the same filenames are used for evaluating predictions of reported_gender
    # and germline_sex_estimate
    RESULTS_FILENAME_LEAD=${FILENAME_LEAD}_${SEED}_${TRANSCRIPT_TAIL_PERCENT_ARRAY[i]}
    CM_SET_FILE=${RESULTS_FILENAME_LEAD}_prediction_details.tsv
    CM_SET=${RESULTS_FILENAME_LEAD}_confusion_matrix.RDS
    SUMMARY_FILE=${RESULTS_FILENAME_LEAD}_two_class_summary.RDS


    for t in ${targetColumns[@]}; do

      # set up results directory for the column being evaluated
      TEST_TARGET_COLUMN=${t}
      RESULTS_OUTPUT_DIRECTORY=${RESULTS}/${TEST_TARGET_COLUMN}

      # model evaluation
       Rscript --vanilla 03-evaluate_model.R \
         --test_expression_file_name ${PROCESSED}/$TEST_EXPRESSION_FILE_NAME \
         --test_targets_file_name ${PROCESSED}/$TEST_TARGETS_FILE_NAME \
         --model_object_file_name ${MODELS}/$MODEL_OBJECT_FILE_NAME \
         --model_transcripts_file_name ${MODELS}/$MODEL_TRANSCRIPTS_FILE_NAME \
         --output_directory $RESULTS_OUTPUT_DIRECTORY \
         --test_target_column $TEST_TARGET_COLUMN \
         --cm_set_file_name $CM_SET_FILE \
         --cm_file_name $CM_SET \
         --summary_file_name $SUMMARY_FILE
    done

  fi

done



Rscript -e "rmarkdown::render('04-present_results.Rmd', 'html_document', params = list(results_dir = '${RESULTS}', \
      model_dir = '${MODELS}', \
      train_target_column = '${TRAIN_TARGET_COLUMN}', \
      target_columns = '${targetColumns_to_pass}'))"



