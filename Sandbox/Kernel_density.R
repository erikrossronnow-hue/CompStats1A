#Assignment 1 - A
library(readr)
library(data.table)
library(ggplot2)
setwd("~/Desktop/CompStat_folder/CompStat_kernel_density")
infrared = as.data.table(read.table("infrared.txt", header=TRUE))
infrared$RAh = as.integer(infrared$RAh)
?source

sourceDir <- function(path, trace = TRUE, ...) {
  op <- options(); on.exit(options(op)) # to reset after each
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
    options(op)
  }
}
sourceDir("R")

##############################################
set.seed(123)
n = 100
p = rgeom(n, 0.5)

x = rnorm(n)    
x[p == 1] = rnorm(sum(p == 1), mean = 7, sd = 1)
# x is now bimodal with two normal dist peaks. 
yy = seq(-4,15,0.001)
#ff_AMISE = lapply(yy, epa_kern_dens, x=x, h=h_AMISE(x))

# cv loglikelihood version with k=N (leave one out)
h_ll = cv_nll(epa_kern_dens, k=length(x), x)$par
h_ucv = UCV(epa_kern_dens, intsqr_epa, k=length(x), x)$par
#ff_CV = lapply(yy, epa_kern_dens, x=x, h = h_ll)

# scale to bw size to compare with R's implementation!
d_ll = density(x, kernel = "epanechnikov", bw = h_ll/ sqrt(5) )
d_ucv = density(x, kernel = "epanechnikov", bw = h_ucv/ sqrt(5) )
# USING GAUSSIAN KERNEL BTW. 
d_R = density(x, kernel = "gaussian", bw = bw.ucv(x))

# plugin estimator comparison
d_AMISE = density(x, kernel = "epanechnikov", bw = h_AMISE(x) / sqrt(5) )
d_R_plugin = density(x, kernel = "epanechnikov", bw = bw.nrd(x) ) #could use nrd0 also

hist(x, prob=TRUE, breaks = 15, main = "mixed normal")
lines(d_R, col="red")
lines(d_ll, col="blue")
lines(d_ucv, col="green")
legend(
  "topright",
  legend = c("gaussian - ucv", "epa - log lik", "epa - ucv"),
  col = c("red", "blue", "green"),
  lty = 1,
  lwd = 1.5,
  bty = "n"  # Optional: removes the outer border box
)

hist(x, prob=TRUE, breaks = 15, main = "mixed normal")
lines(d_R_plugin, col="red", lty=1)
lines(d_AMISE, col="blue", lty=2)
legend(
  "topright",
  legend = c("epa - nrd", "epa - AMISE"),
  col = c("red", "blue"),
  lty = 1,
  lwd = 1.5,
  bty = "n"  # Optional: removes the outer border box
)

# to compare we need to scale our estimates with / sqrt(5) to get bw.
# thats cause Var(K) = 1/5 and R calcs bw = h * sd(K) = h /sqrt(5). 
# see ?bw.ucv
h_AMISE(x) / sqrt(5)
bw.nrd(x)
bw.nrd0(x)
# plugins are pretty close, they dont use same formula exactly tho!

# cv are pretty close as well. Note that bw.ucv uses gaussian kernel!
h_ll / sqrt(5)
h_ucv / sqrt(5)
bw.ucv(x)


#real data
infrared <- read.table("infrared.txt", header = TRUE)
LF12 = log(infrared$F12)
LF12
d_AMISE = density(LF12, kernel = "epanechnikov", bw = h_AMISE(LF12) / sqrt(5) )
d_R_plugin = density(LF12, kernel = "epanechnikov", bw = bw.nrd(LF12) ) #could use nrd0 also

hist(LF12, prob=TRUE, breaks=30)
lines(d_AMISE, lt=1, col="black")
lines(d_R_plugin, col="green", lt=2)

# for the cv stuff:
x = LF12
h_ll = cv_nll(epa_kern_dens, k=length(x), x)$par
h_ucv = UCV(epa_kern_dens, intsqr_epa, k=length(x), x)$par
d_ll = density(x, kernel = "epanechnikov", bw = h_ll/ sqrt(5) )
d_ucv = density(x, kernel = "epanechnikov", bw = h_ucv/ sqrt(5) )
d_R = density(x, kernel = "epanechnikov") #THIS USES bw = nrd0!!!!!
hist(x, prob=TRUE, breaks = 20)
lines(d_R, col="red")
lines(d_ll, col="blue")
lines(d_ucv, col="green")

#more tests: