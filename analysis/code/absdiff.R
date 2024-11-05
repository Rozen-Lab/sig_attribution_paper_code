
absdiff = function(f1, f2) {
browser()
  x1 = read.csv(f1, row.names = 1)
  message(f1, " ", nrow(x1), " ", ncol(x1))
  x2 = read.csv(f2, row.names = 1)
  message(f2, " ", nrow(x2), " ", ncol(x2))
  
  if (!all.equal(dim(x1), dim(x2))) {
    return(Inf)
  }
  return(sum(abs(x1 - x2)))
  
}
