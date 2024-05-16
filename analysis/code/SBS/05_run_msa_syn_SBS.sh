#!/bin/bash
PROJECT_DIR=/home/e0012078/packages/raw/paa_paper_code
DIR1=analysis/code/SBS/msa
DIR2=analysis/code/SBS/msa_opt
BASH_FILE_NAME=run_bash_script_syn_SBS.sh

BASH_FILE1="$PROJECT_DIR"/$DIR1/$BASH_FILE_NAME
BASH_FILE2="$PROJECT_DIR"/$DIR2/$BASH_FILE_NAME

bash $BASH_FILE1
bash $BASH_FILE2