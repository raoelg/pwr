#' What is the power of a t-test?
#'
#' @argument n     the number of observations in each group.
#' @argument delta the effect size.
#' @argument alpha the level of the test.
#' @return   the simulation based power estimate.
#' @example
#' pwr.ttest(10,0.3)
pwr.ttest <-
function(n, delta, alpha=0.05) mean(replicate(1000, t.test(rnorm(n), rnorm(n,delta))$p.value) < alpha)

#' What is the power of a one-way ANOVA?
#'
#' @argument n        the number of observations in each group.
#' @argument nlevels  the number of (levels of the) group(ing variable).
#' @argument mu       the group (-level variable) means.
#' @argument alpha    the level of the test.
#' @argument nsim     the number of simulations on which the power estimate is based.
#' @return   a simulation based estimate of the power.
#' @example
#' pwr.1way(10, 3, sqrt(0.3)*1:3, nsim=500) 
#' power.anova.test(3, 10, 0.3, 1) # approximately the same
pwr.1way <- 
function(n, nlevels, mu=0.3*1:nlevels, alpha = 0.05, nsim=100) {res = rowMeans(replicate(nsim, (fit <<- anova(lm(y~group, transform(data.frame(group=gl(nlevels, n)), y=mu[group] + rnorm(group)) )))$Pr) < alpha); names(res) = rownames(fit); res}

#' What is the power of a multiway ANOVA?
#'
#' @argument n        the number(s) of observations of each row of design.
#' @argument design   a data.frame that specifies the design of the multiway factorial structure.
#' @argument delta    the effect size factor for group mean differences specified in mu.
#' @argument mu       an array with dimensions equal to the levels of the columns of design,
#'                    specifying the additive effects (this should change to allow for 
#'                    interaction effects; probably using an lm object)
#' @argument alpha    level of the tests
#' @nsim              the number of simulations performed for computing the power
#' @return            vector of simulation based power estimates for each of the (interaction) effects
#' @example
#' n = 1:12 * 5
#' power = sapply(n, pwr.factorial, mu = 0.3 * (row(diag(3))+col(diag(3))))
#' require(ggplot2)
#' qplot(n, power[2,], geom=c('point', 'smooth'))
pwr.factorial <- 
function (n, design = expand.grid(A=letters[1:3], B=LETTERS[1:3]), 
    delta = 0.2, mu = delta * slice.index(array(0, dim = sapply(design, 
        nlevels)),1), alpha = 0.05, nsim = 100) 
{
	if(missing(n)) stop("'n' not specified")
    df = design[rep(1:nrow(design), n), ]
    names(df) = names(design)
    formula = as.formula(paste("y ~ ", paste(names(df), collapse="*")))
    # lmfit = lm(y ~ V1*V2, transform(expand.grid(V1=letters[1:3],V2=LETTERS[1:3]), y= rnorm(10*9))); 
    # lmfit$coefficients[] = c(0,0,1,2,0,0,0,0,0)
    # kan zo: lm(simulate(lmfit,10) ~ V1*V2, lmfit$model)
    # of zo: lmfit$call$formula[[2]] = expression(as.matrix(simulate(lmfit,100)))[[1]]; eval(lmfit$call); summary(aov(.Last.value))
    res = rowMeans(replicate(nsim, {
        y = rep(0, nrow(df))
        for (i in ncol(df)) y = y + mu[df[, i], i]
        df$y = y + rnorm(y)
        fit <<- anova(lm(formula, df))
        fit$Pr
    }) < alpha)
    names(res) = rownames(fit);
    res
}

pwr.1way(100,3)