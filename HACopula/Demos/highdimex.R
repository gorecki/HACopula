# The high-dimensional example in R
#
# This script loads the data (U) generated for the high-dimensional example in MATLAB and computes
# a HAC estimate using the penalized Maximum likelihood method from the R package HAC.
#
# NOTE:
# 1) When loading the data, it is assumed that the folder Demos is the working directory.
# 2) The computation takes roughly 6.5 hours on Intel Core 2.83 GHz processor.
#    (to avoid the computation, load the estimate directly using load("highdimex_est.RData"))
#
# Copyright 2017 Jan Górecki

install.packages("HAC")
install.packages("R.matlab")

require("HAC")
require("R.matlab")

data <- readMat("highdimex_data.mat") # assuming Demos as working directory
U = data$U

start.time <- Sys.time()
est.obj = estimate.copula(U, type = 3, method = 4)
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

plot(est.obj)


