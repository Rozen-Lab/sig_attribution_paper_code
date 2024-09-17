
# X = .... # Spectra
# W = .... # Signatures

# H, model = musical.refit.refit(X, W, method = 'likelihood_bidirectional', thresh = 0.001)

# H, model = musical.refit.refit(X, W, method = 'thresh_naive', thresh = 0)

""""
Run this with the top level directory as the working directory.

"""
import sys

def py_run_musical(input_dir, output_dir, seed_in_use):
  import numpy as np
  import pandas as pd
  # import scipy.stats as stats
  # import matplotlib.pyplot as plt
  # yimport seaborn as sns
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

  X = pd.read_table(samples, sep = "\t", index_col=0) # read spectra from sigpro folder
  print(X)
  W = pd.read_table(sigs, sep = "\t", index_col=0) # read signatures from sigpro folder
  print(X)
  H, model = musical.refit.refit(X = X, W = W, method='likelihood_bidirectional', thresh=0.001)
  print(H)
  H.to_csv(output_dir+"/exposures.csv")


print(f"The argument passed is: {sys.argv[1]}")
py_run_musical(sys.argv[1], sys.argv[2], sys.argv[3])
print("done")

