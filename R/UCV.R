UCV = function(kern_dens, integrated_sq, k, x) {
  #k fold CV, first we set up a list of length N with the group numbers:
  N = length(x)
  ind = seq_along(x)
  # cv 
  labels = rep_len(1:k, N) # 1,2,...,k,1,2,...,k,1,2,3 (N long)
  groups = sample(labels)
  
  nll = function(h, x) {
    vals = numeric(N)
    for (i in ind) {
      #create the I^-i set:
      I_minus = ind[groups != groups[i]]
      #use the epa_kern function to calculate each f^-i
      vals[i] = kern_dens(x[i], x[I_minus], h)
    }
    
    #also need the ||f_h||_2^2 part:
    integrand = function(y) {
      sapply(y, function(z) {kern_dens(z, x = x, h = h)^2})
      #integrate needs a vectorized version...
    }
    
    integrated_sq(x = x, h = h) - 2/N * sum(vals)
  }
  
  #start our guess with something that looks like rule of thumb 
  opt = nlminb(start = 1.06*sd(x) * N^(-1/5), objective = nll,
               lower = 1e-8, x=x)
  opt
}
