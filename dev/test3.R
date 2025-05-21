library(BernsteinCopula)
library(MASS)
library(copulaData)
library(ggplot2)


# set.seed(123)
# Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
# data <- MASS::mvrnorm(n = 200, mu = c(0, 0), Sigma = Sigma)
# x <- data[, 1]
# y <- data[, 2]

data("NELS88")
x <- NELS88$Math
y <- NELS88$Reading

data <- bincount2d(x, y, xbreaks = seq(0, 100, 10), ybreaks = seq(0, 100, 10))

Fx <- normaldist(mean=mean(x), sd=sd(x))
Gy <- normaldist(mean=mean(y), sd=sd(y))

m <- 3
n <- 3
R0 <- matrix(1/(m*n), m, n)

opts <- list(
  verbose = TRUE,
  steps = 20
)

result <- fit.copula.group.joint(data$xbreaks, data$ybreaks, data$counts, R_init=R0, Fx=Fx, Gy=Gy, options=opts)

Rhat <- result$R
Fx$params <- result$params1
Gy$params <- result$params2

xg <- seq(min(data$xbreaks), max(data$xbreaks), length.out = 50)
yg <- seq(min(data$ybreaks), max(data$ybreaks), length.out = 50)
Z <- bernstein_joint_pdf_grid(xg, yg, Rhat, Fx, Gy)

df_data <- data.frame(x = x, y = y)

print(
  plot_contour(xg, yg, Z,
             title = "Estimated Joint Density (filled)",
             xlabel = "x", ylabel = "y")
)

print(
  plot_contour(xg, yg, Z, filled = FALSE,
             title = "Estimated Joint Density with Data",
             xlabel = "x", ylabel = "y",
             points = df_data)
)
