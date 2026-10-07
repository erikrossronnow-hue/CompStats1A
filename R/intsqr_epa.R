intsqr_epa = function(x, h) {
  # using computation trick from the book together with closed
  # form solution for norm of epa kernel. calcs ||f||²
  # o is the overlap 2-d thingy, cause if we're outside
  # this then it evals to zero. 
  
  d = abs(outer(x, x, "-")) / h
  o = pmax(2 - d, 0)
  
  vals = 3/160 * o^3 *
    (d^2 + 6 * d + 4)
  
  mean(vals) / h
}
