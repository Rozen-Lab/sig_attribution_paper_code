#!/bin/bash

echo "This might take a while ... grab a coffee/tea/water/wine"
mkdir -p logs

# default parameters
n=-1 # number of samples: -1 for all available ones (warning: PCAWG takes hours on a single machine)
dataset="SIM" # dataset name, change to "PCAWG" for real data
input_path="input_mutation_tables" # path to input data, change to "WGS_PCAWG" for the PCAWG dataset
signatures_path="signature_tables" # path to signatures, change to "signatures_random" for random sigantures
signature_prefix="sigProfiler" # change to sigRandom for random signatures (e.g. SIMrand dataset)

# read the arguments if changing defaults
while [[ "$#" > 0 ]]; do case $1 in
  -n|--number_of_samples) n="$2"; shift;;
  -d|--dataset) dataset="$2"; shift;;
  -i|--input_path) input_path="$2"; shift;;
  -s|--signatures_path) signatures_path="$2"; shift;;
  -p|--signature_prefix) signature_prefix="$2"; shift;;
  *) echo "Unknown parameter passed: $1"; exit 1;;
esac; shift; done

echo "Running on dataset: ${dataset}"
echo "Input data path: ${input_path}"
echo "Signatures path: ${signatures_path}"
echo "Signature prefix used: ${signature_prefix}"

# SBS: 96 context
echo "Running SBS attribution (96 context)..."
nohup python bin/run_NNLS.py -d ${dataset} -n ${n} -p ${signature_prefix} -s ${signatures_path} -i ${input_path} -t SBS -c 96 &> logs/${dataset}_NNLS_SBS_96.log &
nohup python bin/run_NNLS.py -d ${dataset} -n ${n} -p ${signature_prefix} -s ${signatures_path} -i ${input_path} -t SBS -c 96 -x &> logs/${dataset}_NNLS_SBS_96_optimised.log &
# SBS: 192 context
echo "Running SBS attribution (192 context)..."
nohup python bin/run_NNLS.py -d ${dataset} -n ${n} -p ${signature_prefix} -s ${signatures_path} -i ${input_path} -t SBS -c 192 &> logs/${dataset}_NNLS_SBS_192.log &
nohup python bin/run_NNLS.py -d ${dataset} -n ${n} -p ${signature_prefix} -s ${signatures_path} -i ${input_path} -t SBS -c 192 -x &> logs/${dataset}_NNLS_SBS_192_optimised.log &
wait;

echo "All done!"
