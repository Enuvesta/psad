# Kate --------------------------------------------------------------------


data <- c(8, 12, 7, 6, 12)
theta = 1 / mean(c(8, 12, 7, 6, 12))
# 1 / theta = N
N = mean()

d_loglikelihood_kate <- function(theta, d) {
  n = length(d)
  return (n / theta  - (sum(d)  - n) / (1-theta))
}

curve(function(x) sapply(x, d_loglikelihood_kate(x, d=data)), from=0.001, to=10)
f <- function(x) sapply(x, function(x) d_loglikelihood_kate(x, d=data))
curve(f , from=0.001, to=2)
(solution <- uniroot(d_loglikelihood_kate, c(0.01, 0.99), data))
theta <- solution$root

plot(dgeom(seq(1,60), theta), type="l")



# Wind speed - MLE --------------------------------------------------------



data <- read.csv("./S-2/Turbine.csv")
head(data)
windspeed <- data$AveSpeed
windspeed

d_loglikelihood_wind <- function(k, x) {
  alpha <- sum(x^k)
  1 / k  + mean(log(x)) - 1 / alpha * sum( x^k * log(x))
}

(k <- uniroot(d_loglikelihood_wind, c(0, 10), x=windspeed)$root)
lambda <- mean(windspeed ^ k) ^ (1/k)

hist(windspeed, breaks=9, probability = T)
curve(dweibull(x, k, lambda), from=0, to=100,  add = TRUE, col=2, lw=5)


# Boostrap - wind speed probability ---------------------------------------

1 - mean(windspeed < 5)
mean(windspeed > 5)
1 - mean(windspeed <= 5)
1 - pweibull(5, k, lambda)


get_p_simple <- function(data, indicies, p=5) {
  mean(data[indicies] > p)
}


get_p_weibull <- function(data, indicies, p=5) {
  k <- uniroot(d_loglikelihood_wind, c(0, 20), x=data[indicies])$root
  lambda <-  mean(data[indicies] ^ k) ^ (1/k)
  1 - pweibull(p, k, lambda)
}


library(boot)
N = 1000
(bt.simple <- boot(windspeed, get_p_simple, R = N))
(bt.weibull <- boot(windspeed, get_p_weibull, R = N))
(bias.simple <- mean(bt.simple$t) - bt.simple$t0)
(bias.simple <- mean(bt.weibull$t) - bt.weibull$t0)


boot.ci(bt.simple)
boot.ci(bt.weibull)

plot(bt.simple)
plot(bt.weibull)


# Delta method ------------------------------------------------------------

df <- read.csv("./S-2/Delta method/aspirin.csv", row.names = 1)

GetOdds <- function(x, indicies) {
  x <- x[indicies, ]
  p1 <- mean(x[x$aspirin==1, "stroke"])
  p2 <- mean(x[x$aspirin==0, "stroke"])
  p1 / (1-p1)  / p2 * (1-p2)
}

GetLogOdds <- function(x, indicies) {
  log(GetOdds(x, indicies))
}

N = 1000
bt.odds <- boot(df, GetOdds,  N, stype = "i")
bt.logodds <- boot(df, GetLogOdds,  N, stype = "i")


odds <- bt.odds$t0
var.odds <- as.numeric(var(bt.odds$t))

logodds <- bt.logodds$t0
var.logodds <- as.numeric(var(bt.logodds$t))

var.logodds
var.odds * (1 /(odds))^2
