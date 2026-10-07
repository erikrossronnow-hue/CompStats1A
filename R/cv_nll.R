# this function implements k-fold cv using log-likelihood, NOT UCV!
# kern_dens is fhat
# k is cv folds
# x is data
cv_nll = function(kern_dens, k, x) {
  N = length(x)
  # cv
  labels = rep_len(1:k, N) # 1,2,...,k,1,2,...,k,1,2,3 (N long)
  groups = sample(labels)
  ind = seq_along(x)
  vals = numeric(N)
  
  nll = function(h) {
  for (i in ind) {
    #create the I^-i set:
    I_minus = ind[groups != groups[i]]
    #use the density kernel function to calculate each f^-i
    vals[i] = kern_dens(x[i], x[I_minus], h) 
    }
    -sum(log(vals))
  }
  
  h_start = 1.06*sd(x) * N^(-1/5)
  while (!is.finite(nll(h_start))) { h_start = 2 * h_start }
  #just return the opt object, then we can grab par / likelihood 
  nlminb(start = h_start, objective = nll, lower = 1e-8)
}


