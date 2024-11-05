#!/bin/bash

if [ "$(basename "$PWD")" = "sig_attribution_paper_code" ]; then
  echo "Current directory is ok"
else 
  echo "Run this script from sig_attribution_paper_code" >&2
  exit 1
fi

PROJECT_DIR="."
DIR1=analysis/code/ID/msa_default
DIR2=analysis/code/ID/msa_thresholdx10
DIR3=analysis/code/ID/msa_thresholdx100
DIR4=analysis/code/ID/msa_thresholdx1000
BASH_FILE_NAME=run_bash_script_syn_ID.sh

BASH_FILE1="$PROJECT_DIR"/$DIR1/$BASH_FILE_NAME
BASH_FILE2="$PROJECT_DIR"/$DIR2/$BASH_FILE_NAME
BASH_FILE3="$PROJECT_DIR"/$DIR3/$BASH_FILE_NAME
BASH_FILE4="$PROJECT_DIR"/$DIR4/$BASH_FILE_NAME

bash $BASH_FILE1
bash $BASH_FILE2
bash $BASH_FILE3
bash $BASH_FILE4