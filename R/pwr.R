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
  if (length(names(delta)) > 1)
    lmfit$coefficients[names(delta)] = delta
  else
    lmfit$coefficients[] = delta
  lmfit$fitted.values = drop(model.matrix(lmfit) %*% lmfit$coefficients)
  lmfit$call[[2]] = formula
  lmfit$call$formula[[2]] = expression(as.matrix(simulate(lmfit, nsim)))[[1]] # change dependent to nsim simulated response vectors from model
  res = summary(aov(eval(lmfit$call))) # simulate and fit the model and return a summary anova table for each simulated response vector
  attr(res, "lm") = lmfit
  res
  structure(rowMeans(sapply(res, "$.data.frame", "Pr") < alpha), names=gsub("\\s","", rownames(res[[1]])), anovas = res, model = lmfit, design = design, alpha = alpha, formula = formula)
}

pwr.manova <-
  function (n, ndep=2, delta = 0, design = expand.grid(A=letters[1:3], B=factor(1:2)),
            M = ~ A, X = ~ 0, test = c("Pillai", "Wilks", "Hotelling-Lawley",
            "Roy", "Spherical"), Sigma= diag(ndep), alpha = 0.05, nsim = 200)
  {
    if (missing(n)) stop("'n' not specified")
    test = match.arg(test)
    formula = as.formula(paste("cbind(", paste("y",1:ndep, sep="", collapse=", "), ") ~ ", paste(names(df), collapse="*")))
    ngroups = prod(sapply(design, levels))
    df = transform(design, y = matrix(rnorm(n*ngroups*ndep),n*ngroups,ndep)); transform(design, y1 = rep(0,n*ngroups)); for(i in 2:ndep) df[[paste('y',i,sep='')]] = rep(0,n*ngroups);
    lmfit = lm(formula, df)
    lmfit$coefficients[] = 0;
    if (is.list(delta) && is.null(names(delta))) names(delta) = paste('y',1:length(delta),sep='') # make sure listed effects have names
    delta = unlist(delta)
    if (length(delta) > 1)
      lmfit$coefficients[do.call(rbind, strsplit(names(unlist(delta)), "\\."))[,2:1]] = delta
    else
      lmfit$coefficients[] = delta
    lmfit$fitted.values = drop(model.matrix(lmfit) %*% coef(lmfit))
    lmfit$call[[2]] = formula
    ssd   = SSD(lmfit)
    Resid = matrix(rnorm(nsim*length(lmfit$resid)), ncol=ncol(lmfit$resid)) %*% chol(ssd$SSD / ssd$df)
    sims  = array(fitted(lmfit), c(nrow(lmfit$resid), ncol(lmfit$resid), nsim)) + aperm(array(Resid, c(nrow(lmfit$resid), nsim, ncol(lmfit$resid))), c(1,3,2))
    res   = apply(sims, 3, function(y) anova(lm.fit(x=model.matrix(lmfit), y), M=M, X=X, test=test, idata=design))
    structure(rowMeans(sapply(res, "$.data.frame", "Pr") < alpha), names=gsub("\\s","", rownames(res[[1]])), anovas = res, model = lmfit, design = design, alpha = alpha, formula = formula)
  }
