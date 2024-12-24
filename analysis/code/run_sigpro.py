""""
Run this with the top level directory as the working directory.

"""

def run_sigpro(input_dir, output_dir, seed_in_use, context_type):
  import os
  import random
  import glob
  from SigProfilerAssignment import Analyzer

  random.seed(int(seed_in_use))

  print(input_dir)
  sigs = input_dir+"/sigs.tsv"
  samples = input_dir+"/catalog.tsv"
  if not os.path.exists(output_dir):
    os.makedirs(output_dir, exist_ok=True)    
  Analyzer.cosmic_fit(samples=samples, output=output_dir, 
                      context_type=context_type, 
                      signature_database=sigs, genome_build="GRCh37",
                      collapse_to_SBS96=False)

import sys
# The for loop is just for sanity checking
for i, arg in enumerate(sys.argv):
  print(f"Argument {i:>6}: {arg}")


run_sigpro(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
