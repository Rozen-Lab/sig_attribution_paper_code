This repository contains code for the manuscript:

Nanhai Jiang, Yang Wu, Steve Rozen, 
"Benchmarking 13 tools for mutational signature attribution, including a new and improved algorithm"

The R code for the "Presence Attribute Signature Activity" (PASA) algorithm is freely available at https://github.com/steverozen/mSigAct. The version used here is V3.0.1-branch, which can be installed with the R call 

`remotes::install_github(repo = "steverozen/mSigAct", ref = "v3.0.1-branch")`

The PASA algorithm is implemented in the function `mSigAct::PresenceAttributeSigActivity`.
