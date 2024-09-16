
# X = .... # Spectra
# W = .... # Signatures

# H, model = musical.refit.refit(X, W, method = 'likelihood_bidirectional', thresh = 0.001)

# H, model = musical.refit.refit(X, W, method = 'thresh_naive', thresh = 0)

""""
Run this with the top level directory as the working directory.

"""

def run_musical(in_directorytop, seed_in_use):
  import numpy as np
  import pandas as pd
  # import scipy.stats as stats
  # import matplotlib.pyplot as plt
  # import seaborn as sns
  # import matplotlib as mpl
  # import time
  # import scipy as sp
  # import pickle
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
  # We might be able to convert straight from R data structures
  X = pd.read_table(samples, sep = "\t") # read spectra from sigpro folder
  W = pd.read_table(sigs, sep = "\t") # read signatures from sigpro folder
  H, model = musical.refit.refit(X = X, W = W, methods = 'likelihood_bidirectional', thresh=0.001)
  return(H) 

