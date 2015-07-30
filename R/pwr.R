#' Estimate the Power of a T-Test.
#'
#' @param n      The number of observations in each group.
#' @param delta  The effect size.
#' @param alpha  The level of the test.
#' @return  The simulation based power estimate.
#' @examples
#' pwr.ttest(10,0.3)
pwr.ttest <-
function(n, delta, alpha=0.05) mean(replicate(1000, t.test(rnorm(n), rnorm(n,delta))$p.value) < alpha)

#' Estimate the Power of a One-Way ANOVA.
#'
#' @param n        The number of observations in each group.
#' @param nlevels  The number of (levels of the) group(ing variable).
#' @param mu       The group (-level variable) means.
#' @param alpha    The level of the test.
#' @param nsim     The number of simulations on which the power estimate is based.
#' @return   A simulation based estimate of the power.
#' @examples
#' pwr.1way(10, 3, sqrt(0.3)*1:3, nsim=500)
#' power.anova.test(3, 10, 0.3, 1) # approximately the same
pwr.1way <-
  function(n, nlevels, mu=0.3*1:nlevels, alpha = 0.05, nsim = 100) {res = rowMeans(replicate(nsim, (fit <<- anova(lm(y~group, transform(data.frame(group=gl(nlevels, n)), y=mu[group] + rnorm(group)) )))$Pr) < alpha); names(res) = rownames(fit); res}

#' What is the power of a multiway ANOVA?
#'
#' @param n        The number(s) of observations of each row of design.
#' @param design   A data.frame that specifies the design of the multiway factorial structure.
#' @param delta    The effect size factor for group mean differences specified in mu.
#' @param mu       An array with dimensions equal to the levels of the columns of design,
#'                 specifying the additive effects (this should change to allow for
#'                 interaction effects; probably using an lm object).
#' @param alpha    Level of the tests.
#' @param nsim     The number of simulations performed for computing the power.
#' @return  A vector of simulation based power estimates for each of the (interaction) effects.
#' @examples
#' n = 1:12 * 5
#' power = sapply(n, pwr.factorial, delta = c(Ab=0.4, Ac=0.4))
#' ggplot2::qplot(n, power[1,], geom=c('point', 'smooth'))
pwr.factorial <-
function (n, delta = 0.2, design = expand.grid(A=letters[1:3], B=LETTERS[1:2]),
    alpha = 0.05, nsim = 200)
{
  if(missing(n)) stop("'n' not specified")
  df = design[rep(1:nrow(design), n), ]
  names(df) = names(design)
  formula = as.formula(paste("y ~ ", paste(names(df), collapse="*")))
  ngroups = prod(sapply(design, nlevels))
  df = transform(design, y = rnorm(n*ngroups)) # y should be a unique name
  lmfit = lm(formula, df);
  lmfit$coefficients[] = 0;
  if (length(names(delta)) > 0)
    lmfit$coefficients[names(delta)] = delta
  else
    lmfit$coefficients[] = delta
  lmfit$fitted.values = drop(model.matrix(lmfit) %*% lmfit$coefficients)
  lmfit$call[[2]] = formula
  lmfit$call$formula[[2]] = expression(as.matrix(simulate(lmfit, nsim)))[[1]] # change dependent to nsim simulated response vectors from model
  res = summary(aov(eval(lmfit$call))) # simulate and fit the model and return a summary anova table for each simulated response vector
  attr(res, "lm") = lmfit
  res
  structure(rowMeans(sapply(res, "$.data.frame", "Pr") < alpha), names=gsub("\\s","", rownames(res[[1]])), anovas = res, model = lmfit, design = design, alpha = alpha)
}

