# A variant of Nelder-Mead
# Date: 28 Jan 2018
# Author: Saubhik Mukherjee
# Contact: saubhik[dot]mukherjee[at]gmail[dot]com

# Objective function to minimize
objective = function(x) {
    return (x[1]^2 + x[2]^2 )
}

# Lower bound for each coordinate
LB = c(-300, -300)
# Upper bound for each coordinate
UB = c(+300, +300)

# Create simplex from a single starting point
createSimplex = function (x) {
  # x is a single candidate point in the domain of f
  # dim is the dimension of Dom(f)
  dim = length(x)
  xtest_lo <- rep(list(x), dim)
  for (i in 1:dim) xtest_lo[[i]][i] = LB[i]
  xtest_hi <- rep(list(x), dim)
  for (i in 1:dim) xtest_hi[[i]][i] = UB[i]
  xtest = c(xtest_lo, list(x), xtest_hi)
  return (xtest)
}

# Make objective work on lists (simplices)
f = function(x) { unlist(lapply(x, FUN = objective)) }

# Measure the spread of the simplex
spread = function(x) {
  return (sd(unlist(lapply(x, FUN = function(y) {sqrt(sum((y-x[[1]]))^2)}))))
}

nelderMead = function(xtest) {
  # Minimize f(x), where x in Rn. Current test points: xtest = [x1, ..., x_m], where m = n+1
  
  # Standard values
  alpha = 1
  gamma = 2
  rho = 0.5
  sigma = 0.5
  threshold_f = 0.0001
  threshold_s = 0.0001
  
  # Order the test points
  xtest = xtest[order(f(xtest))]
  
  # check for termination
  # check sample standard deviation of values of f at 
  # the test points and also if the test points are close enough
  if (sd(f(xtest)) < threshold_f & spread(xtest) < threshold_s) {
    result = (c(xtest[1], f(xtest[1])))
    names(result) = c("optimal_Point", "minimized_Objective_Value")
    return (result)
  }
  
  # Calculate xo, centroid of all points except x_m
  xo = 0
  for(i in 1:(length(xtest) - 1)) xo = xo + xtest[[i]]
  xo = xo / (length(xtest) - 1)
  
  # Reflection 
  # alpha is a positive constant, determining extent of reflection, from xo
  xr = xo + alpha * (xo - unlist(xtest[length(xtest)]))
  if (any(xr > UB)) xr[xr > UB] = UB[xr > UB]
  if (any(xr < LB)) xr[xr < LB] = LB[xr < LB]
  
  # if reflected point, xr, is better than second worst, xn, but not better than best, then
  # replace worst point, xn = xtest[length(xtest)], by xr, and go to step 1
  if (f(list(xr)) < f(xtest[length(xtest) - 1]) & f(xtest[1]) <= f(list(xr))) 
    return (nelderMead(c(xtest[-length(xtest)], list(xr))))
  
  # Expansion
  # if reflected point, xr, is the best point so far, then compute the expanded point
  if (f(list(xr)) < f(xtest[1])) {
    # gamma is a constant greater than 1, determining the extent to go towards reflected point
    xe = xo + gamma * (xr - xo)
    if (any(xe > UB)) xe[xe > UB] = UB[xe > UB]
    if (any(xe < LB)) xe[xe < LB] = LB[xe < LB]
    if (f(list(xe)) < f(list(xr))) {
      # if expanded point is better than reflected point, 
      # replace the worst point with the expanded point
      return (nelderMead(c(xtest[-length(xtest)], list(xe))))
    } else {
      # replace the worst point with the reflected point
      return (nelderMead(c(xtest[-length(xtest)], list(xr))))
    }
  }
  
  # Contraction
  # Control will come here only if the reflected point is not better than
  # the second worst, i.e., f(xr) >= f(xn)
  # Compute the contracted point
  # rho is positive constant, but less than or equal to half
  # rho determines how much to move towards xm, from xo.
  xc = xo + rho  * (unlist(xtest[length(xtest)]) - xo)
  if (any(xc > UB)) xc[xc > UB] = UB[xc > UB]
  if (any(xc < LB)) xc[xc < LB] = LB[xc < LB]
  # if contracted point, xc, is better than worst point,
  # replace worst point with contracted point
  if (f(list(xc)) < f(xtest[length(xtest)]))
    return (nelderMead(c(xtest[-length(xtest)], list(xc))))
  
  # Shrink
  # Shrink all points except the best, xtest[1]
  # towards xtest[1]
  # sigma controls the how much the points should 
  # move towards the best point, xtest[1]
  xtest[-1] = lapply(xtest[-1], FUN = function(y) {xtest[[1]] + sigma * (y - xtest[[1]])})
  return (nelderMead(xtest))
  
}

# Test
result = nelderMead(createSimplex(c(1,3)))
print(result)

# Result will be:
# $optimal_Point
# [1] -1.417809e-05  3.421117e-05
# 
# $minimized_Objective_Value
# [1] 1.371423e-09