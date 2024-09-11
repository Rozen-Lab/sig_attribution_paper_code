
X = .... # Spectra
W = .... # Signatures

H, model = musical.refit.refit(X, W, method = 'likelihood_bidirectional', thresh = 0.001)

H, model = musical.refit.refit(X, W, method = 'thresh_naive', thresh = 0)

""""
Run this with the top level directory as the working directory.

"""

def run_musical(directory, seed_in_use):
  import numpy as np
  import scipy.stats as stats
  import matplotlib.pyplot as plt
  import seaborn as sns
  import matplotlib as mpl
  import pandas as pd
  import time
  import scipy as sp
  import pickle
  import musical

  
  import os
  import random
  import glob


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
