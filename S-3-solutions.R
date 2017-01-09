# Binomial test --------------------------------------------------------------------

get_binom_pvalue <- function(np, nq, p0 = 0.5) {
  1 - pbinom(np - 1, np + nq, p0)  ## p(x >= np)
}

get_binom_pvalue(10, 3)
binom.test(c(10, 3), p = 0.5, alternative = 'g')










# Wald test --------------------------------------------------------------------

# 1. 
x1 <- c(3, 5, 6, 6, 7, 10, 13, 15, 18, 22)
x2 <- c(9, 12)
lambda_mle1 <- mean(x1)
lambda_mle2 <- mean(x2)
print(c(lambda_mle1, lambda_mle2))

# log(L) = -n * lambda + log(lambda) * sum(x_i) - sum(log(x_i !))
poisson.loglikelihood <- function(lambda, x) {
  n <- length(x)
  -n * lambda + log(lambda) * sum(x) - sum(lfactorial(x))
}

# 2. 

get_lr1 <- function(lambda) {
  poisson.loglikelihood(lambda, x1) - poisson.loglikelihood(lambda_mle1, x1)
}

get_lr2 <- function(lambda) {
  poisson.loglikelihood(lambda, x2) - poisson.loglikelihood(lambda_mle2, x2)
}

curve(get_lr1, from=3, to=25, xlab=expression(lambda))
curve(get_lr2, add=TRUE, lty="dashed")
legend(x=12, y=-40, legend=c("Sample 1", "Sample 2", "MLE"), lty=c("solid", "dashed", "dotted"), bty="n")

abline(v=lambda_mle1, lty="dotted")  # MLE for both samples

# 3. 

# L''(lambda) = sum(x) / n^2
poisson.mle_variance <- function(x) {
  n <- length(x)
  sum(x) /(n ^ 2)
}
poisson.mle_variance(x1)
# [1] 1.05
poisson.mle_variance(x2)
# [1] 5.25

# For the second sample the variance is higher

# 4. 
lambda_to_test <- 9
lambda_mle <- mean(x1)
var <- poisson.mle_variance(x1)
Z0 <- (lambda_to_test - lambda_mle) / sqrt(var)

Z0
pnorm(Z0)
# 0.07161745 < 0.05 -> don't reject













# Likelihood ratio test --------------------------------------------------------------------

# 1.
curve(poisson.loglikelihood(x, data), from = 2, to = 30,
      xlab = expression(lambda), ylab = "Log-likelihood")
maxlik.mle <- poisson.loglikelihood(lambda_mle, x1)

# Paint interval
segments(x0 = lambda_mle, y0 = -140, x1 = mle, y1 = maxlik.mle, lty="dotted") 
#segments(x0 = 0, y0 = maxlik.mle, x1 = mle, y1 = maxlik.mle, lty="dotted")

# 2. 
abline(h = maxlik.mle - qchisq(p = 0.95, df = 1) / 2, lty = "dashed")


# 3. 
op.func<-function(lambda, data, alpha, maxlikMLE) {
  -2 * (poisson.loglikelihood(lambda, data) - maxlikMLE) - qchisq(p = 1 - alpha, df = 1)
}

# Intersect the dashed line with the curve with uniroot
lower <- uniroot(f = op.func, interval = c(0, mle), data = data, alpha = 0.05, maxlikMLE = maxlik.mle)
upper <- uniroot(f = op.func, interval = c(mle, 30), data = data, alpha = 0.05, maxlikMLE = maxlik.mle)
conf_int <- c(lower$root, upper$root)
print(conf_int)
abline(v = conf_int)

#check if we have correctly found mle
optimize(f = op.func, interval = c(0, 28), data=data, alpha = 0.05, maxlikMLE = maxlik.mle)















# Score test --------------------------------------------------------------------

# 1.

# log(L) = -n * lambda + log(lambda) * sum(x_i) - sum(log(x_i !))
# dlog(L) = -n  + 1 / lambda * sum(x_i)
poisson.score <- function(lambda, data) {
  n <- length(data)
  -n + 1 / lambda * sum(data)
}

poisson.fisherI <- function(lambda, data) {
  1 / lambda ^ 2 * sum(data)
} 

lambda_mle <- mean(x1)
lambda_to_test <- 9.0
test_statistic_score <- poisson.score(lambda_to_test, data) / sqrt(poisson.fisherI(lambda_to_test, data))
test_statistic_score
1 - pnorm(test_statistic_score)
# 0.07161745 < 0.05 -> don't reject

# 2. 
test_statistic_wald <- (lambda_mle - lambda_to_test) / sqrt(poisson.mle_variance(data))
1 - pnorm(test_statistic_wald)
