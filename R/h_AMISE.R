h_AMISE = function(x) {
  N = length(x)
  s_tilde = min(sd(x), IQR(x)/1.34)
  C = ((8 * sqrt(pi) * 3/5)/(3 * (1/5)**2))**(1/5)
  #h
  C * s_tilde * N**(-1/5) 
}