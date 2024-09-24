library(tidyverse)
library(magrittr)
foo = data.table::fread("analysis/summary/SBS/assessment_each_sample_SBS.csv")

foo = dplyr::filter(foo, Tool == "siglasso") 
dim(foo)
foo = arrange(foo, Combined)
worst_c = foo$Sample.ID[c(1:3,5)]
foo = arrange(foo, prec)

foogt = mSigTools::read_exposure("synthetic_data/SBS/ground.truth.syn.exposures.csv")
foogt[ , worst_c]

fooinf = mSigTools::read_exposure("analysis/raw_output/SBS/siglasso/syn/inferred_exposures.csv")
fooinf[ , "Skin.Melanoma..S.14", drop = F]

fooinf[ , "Skin.Melanoma..S.21", drop = F]
