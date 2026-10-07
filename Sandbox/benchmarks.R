setwd("~/Desktop/CompStat_folder/CompStat_kernel_density")
for (p in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(p)

# Compilation is outside the timings. The same seed gives each version the same folds.
out = list()
fits = list()
for (n in c(100, 500, 1000)) {
  set.seed(123)
  p = rgeom(n, 0.5)
  x = rnorm(n)
  x[p == 1] = rnorm(sum(p == 1), 10, 1)

  for (i in seq_len(5)) for (f in list(h_AMISE, h_AMISE_vec, h_AMISE_cpp)) f(x)
  b = bench::mark(r = h_AMISE(x), vec = h_AMISE_vec(x), cpp = h_AMISE_cpp(x), iterations = 1000, check = FALSE, filter_gc = FALSE)
  out[[length(out) + 1]] = data.frame(n = n, k = NA, f = "AMISE", version = c("r", "vec", "cpp"), ms = as.numeric(b$median) * 1000, mb = as.numeric(b$mem_alloc) / 1024^2)

  for (k in c(5, n)) {
    for (j in seq_len(3)) {
      f = list(cv_nll, cv_nll_vec, cv_nll_cpp)[[j]]
      s = list(UCV, UCV_vec, UCV_cpp)[[j]]
      set.seed(1)
      a = f(epa_kern_dens, k, x)
      set.seed(1)
      u = s(epa_kern_dens, intsqr_epa, k, x)
      fits[[length(fits) + 1]] = data.frame(n = n, k = k, version = c("r", "vec", "cpp")[j], ll_h = a$par, ll_f = a$objective, ll_code = a$convergence, ucv_h = u$par, ucv_f = u$objective, ucv_code = u$convergence)
    }

    b = bench::mark(r = {set.seed(1); cv_nll(epa_kern_dens, k, x)}, vec = {set.seed(1); cv_nll_vec(epa_kern_dens, k, x)}, cpp = {set.seed(1); cv_nll_cpp(epa_kern_dens, k, x)}, iterations = 5, check = FALSE, filter_gc = FALSE)
    out[[length(out) + 1]] = data.frame(n = n, k = k, f = "cv_nll", version = c("r", "vec", "cpp"), ms = as.numeric(b$median) * 1000, mb = as.numeric(b$mem_alloc) / 1024^2)

    b = bench::mark(r = {set.seed(1); UCV(epa_kern_dens, intsqr_epa, k, x)}, vec = {set.seed(1); UCV_vec(epa_kern_dens, intsqr_epa, k, x)}, cpp = {set.seed(1); UCV_cpp(epa_kern_dens, intsqr_epa, k, x)}, iterations = 5, check = FALSE, filter_gc = FALSE)
    out[[length(out) + 1]] = data.frame(n = n, k = k, f = "UCV", version = c("r", "vec", "cpp"), ms = as.numeric(b$median) * 1000, mb = as.numeric(b$mem_alloc) / 1024^2)
    print(c(n = n, k = k))
  }
}

out = do.call(rbind, out)
fits = do.call(rbind, fits)
print(out, row.names = FALSE, digits = 4)
print(fits, row.names = FALSE, digits = 7)
dir.create("work", showWarnings = FALSE)
saveRDS(list(times = out, fits = fits, session = sessionInfo()), "work/benchmarks.rds")

# Timings below were measured before adding the finite-start loop; rerun this file for updated timings.
# At n = k = 500, cv_nll took about 240 ms in R, 165 ms in vec and 122 ms in cpp.
# UCV took about 326, 287 and 260 ms. Its matrix-based integral is still shared by all three.
# AMISE was about 0.043 ms in R, 0.041 ms in vec and 0.008 ms in cpp at n = 500.
# AMISE already uses vectorized sd/IQR; the shorter R version does not have much to gain.
# The CV versions save work by building each training set once. vapply still loops over observations.
# C++ still calls the R kernel, so these are modest gains. It is not a fully native density calculation.
# At n = 1000 and k = 5, vec UCV was slower than the original in this run. Faster is not guaranteed.
# mb is total R allocation per fit, not peak RAM. Keeping training sets uses O(n*k) storage.
# All timed fits had finite objectives and matched across versions. Seed 1 already had finite starts here.
