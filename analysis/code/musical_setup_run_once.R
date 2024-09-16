# Only run once

library(reticulate)


# At command line git clone https://github.com/parklab/MuSiCal.git in conda environment musical2
# setwd("~/github/") # adjust for the parent directory of the MuSiCal repo clone
# repl_python(input = "!pip install ./MuSiCal")

use_condaenv("musical2")
source_python("analysis/code/musical_analysis.py")
py_run_string('run_musical(".", 123)')
