# Quantile treatment effect as a function of outcome value in control group 
#  using bootstrap.

# Matti Pirinen & Harri Hemila, version Sep 21, 2022.

# The source code in this file is published under 
# the MIT license listed at the end of this file.

# Contains two functions:
#
# qte_by_outcome <- function(Treatment, Control, at = NULL, 
#                           qte.conf = 0.95, B = 2000, 
#                           K = length(Control),
#                           verbose = TRUE)
#
# plot_qte_by_outcome <- function(x, plot.ci = TRUE, plot.rte = FALSE,
#                                pch = 19, col = "black",
#                                xlim = NULL, ylim = NULL, 
#                                xlab = NULL, ylab = NULL,
#                                xaxs = "r", yaxs = "r",
#                                xaxp = NULL, yaxp = NULL)
#
# For examples of usage, see file 'qte_by_outcome_examples.R'.  


qte_by_outcome <- function(Treatment, Control, at = NULL, 
                           qte.conf = 0.95, B = 2000, 
                           K = length(Control), 
                           verbose = TRUE){
  
  #Goal: 
  # We have observed survival/recovery times, or other continuous outcome values, 
  # in the treatment group and in the control group. 
  # Sample sizes can differ between the groups.
  # We want to estimate the treatment effect (TE) as a function of outcome value
  # at a given grid of values ('at') in the control group.
  
  #Approach:
  # Fix a grid of outcome values in control group at which
  #  TEs are computed (parameter 'at').
  #  Recommendation: 
  #   Do not try to estimate TEs at values that are outside
  #   empirical quantile levels (5/N, 1 - (5/N)), 
  #   where N is the sample size of controls because of low accuracy at the tails. 
  #
  # Estimate the quantile treatment effects as difference between the quantiles 
  #  of Treatment group and Control group at 'K' quantile levels 
  #  1/(K+1),...,K/(K+1). 
  #  Make a piecewise linear approximation of TE
  #  as a function of outcome value in controls, 
  #  based on the observed K quantiles. 
  #  Use that approximation
  #  to estimate the TE at the chosen grid of outcome values in controls.
  #
  # Use bootstrapping to estimate the uncertainty of TE at each grid point,
  #  by resampling with replacement a set of Treatment and Control group values that
  #  have the same sample size as the original Treatment and Control groups, respectively,
  #  and applying the procedure explained above to each bootstrap sample.
  #  Note that bootstrapping accounts for uncertainty of quantiles of both groups.
  # 
  # The final estimate of TE at a grid point is the mean over bootstrap samples
  #  and the confidence interval for TE is estimated from the quantiles of the bootstrap sample.
  
  #INPUT:  
  # 'Treatment', vector of outcome values (e.g. survival times) in Treatment group. 
  # 'Control', vector of outcome values in Control group.
  # 'at', vector of outcome values in controls at which TE is estimated.
  #       By default, between 10 and 20 points between the empirical quantile levels 
  #       of 5/N and 1 - (5/N) in controls where N = length(Control).
  # 'qte.conf' (default 0.95), confidence level for TE
  # 'B' (default 2000), number of bootstrap samples
  # 'K' (default 'length(Control)), number of quantiles where TE function is estimated.
  # 'verbose' (default TRUE) if TRUE, prints parameters and recommendations on console.
  
  #OUTPUT:
  #  data frame with 7 columns
  #1 'at', outcome value in control group
  #2 'qte', (quantile) treatment effect
  #3 'qte.low', lower CI end point for 'qte' 
  #4 'qte.up', upper CI end point for 'qte'
  #5 'rte', relative treatment effect with respect to control outcome value
  #6 'rte.low', lower CI end point for 'rte' 
  #7 'rte.up', upper CI end point for 'rte'
  
  if( K < 2 ) stop("K should be at least 2.")
  if( qte.conf <= 0 | qte.conf >= 1 ) stop("Value \'qte.conf\' not valid.")
  if( B < 1 ) stop("B is not positive integer.")
  C.n = length(Control)
  T.n = length(Treatment)
  
  #recommended to keep value in 'at' within the range 'x'
  x = as.numeric(quantile(Control, prob = c(5/C.n, 1 - 5/C.n)))
  if( is.null(at) ){ # Set 'at' to its default value 
    # that has between 10 and 20 values from interval between empirical quantiles
    # of Control vector of levels 5/C.n and 1-5/C.n.
    # The step size is set to the appropriate power of 2.
    d = 2^(floor(log2((x[2] - x[1])/10))) # step size
    at = d*seq(ceiling(x[1]/d), floor(x[2]/d), 1)
    # Check the endpoints for rounding errors and add possibly missing endpoints:
    tolerance = 1e-10 #tolerance for rounding error
    if(at[1] - d >= x[1] - tolerance){
      at = c(at[1] - d, at)} # add left endpoint
    if(at[length(at)] + d <= x[2] + tolerance) {
      at = c(at, at[length(at)] + d)} # add right end point
  }
  
  if( verbose ){
    cat("Running qte_by_outcome() with the following parameters.","\n")
    cat(paste("length(Treatment):",T.n),"\n")
    cat(paste("length(Control):",C.n),"\n")
    cat(paste("K:",K),"\n")
    cat(paste("B:",B),"\n")
    cat(paste("qte.conf:",qte.conf),"\n")
    cat(paste("at:",paste(signif(at,3),collapse=", ")),"\n")
    cat(paste0("If you choose to change 'at' values, a recommended interval is [",
               signif(x[1],3),", ",signif(x[2],3),"].\n"))
  }
  probs = (1:K)/(K+1) #quantile levels at which the TE-function is estimated
  C.b = matrix(NA, nrow = B, ncol = K) #bootstrapped control quantiles
  T.b = matrix(NA, nrow = B, ncol = K) #bootstrapped treatment quantiles
  TE.b = matrix(NA, nrow = B, ncol = length(at)) #bootstrapped approximations of TEs at grid 'at'
  for (ii in 1:B){
    C.b[ii,] = as.numeric(quantile(sample(Control, size = C.n, replace = T), prob = probs))
    T.b[ii,] = as.numeric(quantile(sample(Treatment, size = T.n, replace = T), prob = probs))
    TE.b[ii,] = as.numeric(approx(C.b[ii,], T.b[ii,] - C.b[ii,], xout = at, rule = 2, ties = mean)$y)
  }

  qte = apply(TE.b, 2, mean) # Mean as the final estimate 
  qte.low = as.numeric(apply(TE.b, 2, function(X){quantile(X, (1 - qte.conf)/2)}))
  qte.up = as.numeric(apply(TE.b, 2, function(X){quantile(X, (1 + qte.conf)/2)}))

  res = data.frame(
    at = at,
    qte = qte, 
    qte.low = qte.low,
    qte.up = qte.up,
    rte = qte/at,
    rte.low = qte.low/at,
    rte.up = qte.up/at)

  return(res)
}


plot_qte_by_outcome <- function(x, plot.ci = TRUE, plot.rte = FALSE,
                                pch = 19, col = "black",
                                xlim = NULL, ylim = NULL, 
                                xlab = NULL, ylab = NULL,
                                xaxs = "r", yaxs = "r",
                                xaxp = NULL, yaxp = NULL){
  
  #Print treatment effect on direct or relative scale
  #INPUT:
  # 'x' data.frame returned by function qte_by_outcome( )
  # plot.ci, if TRUE, plots confidence intervals around estimates
  # plot.rte, if TRUE, plots relative treatment effects, otherwise plots actual treatment effects
  # pch, col, xlim, ylim, xlab, ylab, xaxs, yaxs, xaxp, yaxp: 
  #      standard plotting parameters with sensible defaults.
  
  #Use percentages for relative effects
  if(plot.rte) x[,c("rte","rte.low","rte.up")] = 100*x[,c("rte","rte.low","rte.up")] 
  
  if(is.null(xlim)) {
    xlim = range(x$at) + c(-1,1) * 0.05 * diff(range(x$at))}
  if(is.null(ylim)) {
    if(!plot.rte){
      ylim = range(c(x$qte.low,x$qte.up)) + c(-1,1) * 0.05 * diff(range(c(x$qte.low,x$qte.up)))}
    else{
      ylim = range(c(x$rte.low,x$rte.up)) + c(-1,1) * 0.05 * diff(range(c(x$rte.low,x$rte.up)))}
  }
  if(is.null(xlab)) xlab = "Outcome in controls"
  if(is.null(ylab)){
    if(!plot.rte) ylab = "Treatment effect"
    else ylab = "Relative treatment effect (%)"
  }
  
  plot(NULL, 
       xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab,
       xaxs = xaxs, yaxs = yaxs,
       xaxp = xaxp, yaxp = yaxp)
  grid()
  abline(h = 0, lty = 2)
  
  #add intervals
  if(plot.ci){
    if(!plot.rte){y.low = x$qte.low; y.up = x$qte.up}
    else{y.low = x$rte.low; y.up = x$rte.up}
    arrows(x$at, y.low, 
           x$at, y.up, 
           col = col, code = 3, angle = 90, length = 0.0)
  }
  
  #add points
  if(!plot.rte) points(x$at, x$qte, pch = pch, col = col)
  else points(x$at, x$rte, pch = pch, col = col)
  
}


#MIT License

#Copyright (c) 2022 Matti Pirinen, Harri Hemila

#Permission is hereby granted, free of charge, to any person obtaining a copy
#of this software and associated documentation files (the "Software"), to deal
#in the Software without restriction, including without limitation the rights
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
#copies of the Software, and to permit persons to whom the Software is
#furnished to do so, subject to the following conditions:
  
#  The above copyright notice and this permission notice shall be included in all
#copies or substantial portions of the Software.

#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#SOFTWARE.
