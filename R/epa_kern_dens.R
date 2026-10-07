epa_kern_dens = function(y, x, h) {
  N = length(x)
  f = (y-x)/h
  x = x[(-1 < f) & (f < 1)]
  vals = 1-((y-x)/h)^2
  
  3/(4*N*h) * sum(vals)
}
