Rcpp::sourceCpp("src/bandwidth.cpp")

cv_nll_cpp = function(kern_dens, k, x) {
  n = length(x)
  g = sample(rep_len(seq_len(k), n))
  z = lapply(seq_len(k), function(j) x[g != j])
  f = function(h) nll_cpp(h, kern_dens, g, z, x)
  h = 1.06 * sd(x) * n^(-1/5)
  while (!is.finite(f(h))) { h = 2 * h }
  nlminb(start = h, objective = f, lower = 1e-8)
}

UCV_cpp = function(kern_dens, integrated_sq, k, x) {
  n = length(x)
  g = sample(rep_len(seq_len(k), n))
  z = lapply(seq_len(k), function(j) x[g != j])
  f = function(h) ucv_cpp(h, kern_dens, integrated_sq, g, z, x)
  h = 1.06 * sd(x) * n^(-1/5)
  nlminb(start = h, objective = f, lower = 1e-8)
}
