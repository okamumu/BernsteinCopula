library(BernsteinCopula)
library(MASS)
library(ggplot2)


set.seed(123)
Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
data <- MASS::mvrnorm(n = 200, mu = c(0, 0), Sigma = Sigma)
x <- exp(data[, 1])
y <- exp(data[, 2])


Fx <- normaldist(mean = mean(x), sd = sd(x))
Gy <- normaldist(mean = mean(y), sd = sd(y))
# Fx <- exponentialdist(rate = 1/mean(x))
# Gy <- exponentialdist(rate = 1/mean(y))

m <- 4
n <- 4
R0 <- matrix(1 / (m * n), m, n)

result <- emloop_with_F(
  x = x,
  y = y,
  R_init = R0,
  F = Fx,
  G = Gy,
  options = list(verbose=TRUE)
)

print("Estimated parameters only from marginals:")
print("Fx params:")
print(Fx$params)
print("Gy params:")
print(Gy$params)

print("Estimated parameters with joint copula:")
print("Fx params:")
print(result$params1)
print("Gy params:")
print(result$params2)

