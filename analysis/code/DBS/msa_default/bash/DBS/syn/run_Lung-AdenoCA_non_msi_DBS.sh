#!/bin/bash
PROJECT_DIR=/home/e0012078/packages/sig_attribution_paper_code
MSA_NF_FILE=common_code/MSA/run_auto_optimised_analysis.nf
DATA_DIR=analysis/raw_output/DBS/msa_default/syn/input/Lung-AdenoCA/non_msi
DATA_NAME=syn_data
OUTPUT_DIR=analysis/raw_output/DBS/msa_default/syn/output/Lung-AdenoCA/non_msi/raw
TOOL=conda
TYPE=DBS

ABS_MSA_NF_FILE="$PROJECT_DIR"/$MSA_NF_FILE
ABS_DATA_DIR="$PROJECT_DIR"/$DATA_DIR
ABS_OUTPUT_DIR="$PROJECT_DIR"/$OUTPUT_DIR

nice nextflow run $ABS_MSA_NF_FILE \
 -profile $TOOL --dataset $DATA_NAME --input_tables $ABS_DATA_DIR \
 --mutation_types $TYPE --output_path $ABS_OUTPUT_DIR \
 --signature_prefix $DATA_NAME --signature_tables $ABS_DATA_DIR/$DATA_NAME \

TAR_FILE_DIR="$ABS_OUTPUT_DIR"/output_tables/
TAR_FILE_NAME=MSA_output_"$DATA_NAME".tar.gz
cd $TAR_FILE_DIR
tar xf $TAR_FILE_NAME
EXPOSURE_FILE_DIR="$TAR_FILE_DIR"/$DATA_NAME
EXPOSURE_FILE_NAME=output_"$DATA_NAME"_"$TYPE"_mutations_table.csv
PRUNED_EXPOSURE_FILE_NAME=pruned_attribution_"$DATA_NAME"_"$TYPE"_abs_mutations.csv
cp $EXPOSURE_FILE_DIR/$EXPOSURE_FILE_NAME $ABS_OUTPUT_DIR/..
cp $EXPOSURE_FILE_DIR/$PRUNED_EXPOSURE_FILE_NAME $ABS_OUTPUT_DIR/..
cp $ABS_OUTPUT_DIR/nf-pipeline_info/MSA-nf_report.html $ABS_OUTPUT_DIR/..
