library(BernsteinCopula)
library(mapfit)
library(MASS)
library(ggplot2)

set.seed(123)
Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
data <- MASS::mvrnorm(n = 20, mu = c(5, 5), Sigma = Sigma)
x <- data[, 1]
y <- data[, 2]


xres <- phfit.point(ph=cf1(10), x=x)
Fx <- phdist(alpha = xres$model$alpha(), rate = xres$model$rate())
yres <- phfit.point(ph=cf1(10), x=y)
Gy <- phdist(alpha = yres$model$alpha(), rate = yres$model$rate())

m <- 4
n <- 4
R0 <- matrix(1 / (m * n), m, n)

result <- fit.copula.point.joint(
  x = x,
  y = y,
  R_init = R0,
  Fx = Fx,
  Gy = Gy,
  options=list(verbose=FALSE)
)

cat("\nEstimated parameters with joint copula:\n")
cat("Fx params:\n")
print(result$params1)
cat("Gy params:\n")
print(result$params2)

Rhat <- result$R
Fx$params <- result$params1
Gy$params <- result$params2

xg <- seq(min(x), max(x), length.out = 50)
yg <- seq(min(y), max(y), length.out = 50)
Z <- bernstein_joint_pdf_grid(xg, yg, Rhat, Fx, Gy)

df_data <- data.frame(x = x, y = y)

plot_contour(xg, yg, Z,
             title = "Estimated Joint Density (filled)",
             xlabel = "x", ylabel = "y")

plot_contour(xg, yg, Z, filled = FALSE,
             title = "Estimated Joint Density with Data",
             xlabel = "x", ylabel = "y",
             points = df_data)
