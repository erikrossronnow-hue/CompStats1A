setwd("~/Desktop/CompStat_folder/CompStat_kernel_density")
for (p in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(p)

set.seed(123)
n = 500
p = rgeom(n, 0.5)
x = rnorm(n)
x[p == 1] = rnorm(sum(p == 1), 10, 1)
k = n

f = list(AMISE_r = function() h_AMISE(x), AMISE_vec = function() h_AMISE_vec(x), AMISE_cpp = function() h_AMISE_cpp(x), cv_nll_r = function() cv_nll(epa_kern_dens, k, x), cv_nll_vec = function() cv_nll_vec(epa_kern_dens, k, x), cv_nll_cpp = function() cv_nll_cpp(epa_kern_dens, k, x), UCV_r = function() UCV(epa_kern_dens, intsqr_epa, k, x), UCV_vec = function() UCV_vec(epa_kern_dens, intsqr_epa, k, x), UCV_cpp = function() UCV_cpp(epa_kern_dens, intsqr_epa, k, x))
m = c(100000, 100000, 100000, 20, 20, 20, 10, 10, 10)
dir.create("work/profiles", recursive = TRUE, showWarnings = FALSE)
out = list()

# Repeat the short calls so Rprof has something to sample. Compile and warm up first.
for (i in seq_along(f)) {
  a = f[[i]]
  a()
  p = paste0("work/profiles/", names(f)[i], ".out")
  Rprof(p, interval = 0.01, memory.profiling = TRUE)
  for (j in seq_len(m[i])) a()
  Rprof(NULL)
  s = summaryRprof(p, memory = "both")
  out[[names(f)[i]]] = s
  print(names(f)[i])
  print(head(s$by.self, 8))
  print(head(s$by.total, 8))
}

saveRDS(list(profiles = out, session = sessionInfo()), "work/profiling.rds")

# Measured with n = k = 500 before adding the finite-start loop. Rerun for updated profiles.
# AMISE: IQR took about 83% in both R versions. The vec rewrite does the same work.
# cv_nll: the original spent about 68% in kern_dens and 32% in the rest of its objective.
# In vec, kern_dens took about 91%. Reusing the training sets removed much of that other work.
# In cpp, the R density callback took about 90%, so moving just the loop to C++ has limits.
# UCV: integrated_sq took about 54% in the original and 53% in vec, mostly matrix work.
# In cpp, the two R callbacks together took about 96%. The integral is still calculated in R.
# Rprof shows C++ as one call, not individual C++ lines. <Anonymous> contains the R callbacks.
