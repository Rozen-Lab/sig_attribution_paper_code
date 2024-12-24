#!/bin/bash

if [ "$(basename "$PWD")" = "sig_attribution_paper_code" ]; then
  echo "Current directory is ok"
else 
  echo "Run this script from sig_attribution_paper_code" >&2
  exit 1
fi

PROJECT_DIR="."
DIR1=analysis/code/ID/msa_default
BASH_FILE_NAME=run_bash_script_syn_ID.sh

BASH_FILE1="$PROJECT_DIR"/$DIR1/$BASH_FILE_NAME

bash $BASH_FILE1