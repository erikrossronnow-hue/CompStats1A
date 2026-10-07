h_AMISE_vec = function(x) (40 * sqrt(pi))^(1/5) * min(sd(x), IQR(x)/1.34) * length(x)^(-1/5)

cv_nll_vec = function(kern_dens, k, x) {
  n = length(x)
  g = sample(rep_len(seq_len(k), n))
  z = lapply(seq_len(k), function(j) x[g != j])
  f = function(h) {
    v = vapply(seq_len(n), function(i) kern_dens(x[i], z[[g[i]]], h), numeric(1))
    -sum(log(v))
  }
  h = 1.06 * sd(x) * n^(-1/5)
  while (!is.finite(f(h))) { h = 2 * h }
  nlminb(start = h, objective = f, lower = 1e-8)
}

UCV_vec = function(kern_dens, integrated_sq, k, x) {
  n = length(x)
  g = sample(rep_len(seq_len(k), n))
  z = lapply(seq_len(k), function(j) x[g != j])
  f = function(h) {
    v = vapply(seq_len(n), function(i) kern_dens(x[i], z[[g[i]]], h), numeric(1))
    integrated_sq(x = x, h = h) - 2/n * sum(v)
  }
  h = 1.06 * sd(x) * n^(-1/5)
  nlminb(start = h, objective = f, lower = 1e-8)
}
