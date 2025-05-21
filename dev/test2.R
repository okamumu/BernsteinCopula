library(BernsteinCopula)
library(mapfit)
library(ggplot2)
library(copulaData)


data("NELS88")

x <- NELS88$Math
y <- NELS88$Reading

# Fx <- normaldist(mean=mean(x), sd=sd(x))
# Gy <- normaldist(mean=mean(y), sd=sd(y))

xres <- phfit.point(ph=cf1(10), x=x)
yres <- phfit.point(ph=cf1(10), x=y)
Fx <- phdist(alpha = xres$model$alpha(), rate = xres$model$rate())
Gy <- phdist(alpha = yres$model$alpha(), rate = yres$model$rate())
cat("Estimated parameters only from marginals:\n")
cat("Fx params:\n")
print(Fx$params)
cat("Gy params:\n")
print(Gy$params)

m <- 10
n <- 10
R0 <- matrix(1 / (m * n), m, n)

result <- emloop_with_F(
  x = x,
  y = y,
  R_init = R0,
  Fx = Fx,
  Gy = Gy,
  options=list(verbose=FALSE)
)

cat("Estimated parameters from copula:\n")
cat("R params:\n")
print(result$R)
cat("Fx params:\n")
print(result$params1)
cat("Gy params:\n")
print(result$params2)

xgrid <- seq(0, 100, length.out = 200)

df_fx <- data.frame(x = xgrid, density = Fx$pdf(xgrid))
print(
  ggplot(df_fx, aes(x = x, y = density)) +
  geom_line(color = "steelblue") +
  labs(title = "Updated Marginal of X", x = "X", y = "Density") +
  theme_minimal()
)
df_gy <- data.frame(y = xgrid, density = Gy$pdf(xgrid))

print(
  ggplot(df_gy, aes(x = y, y = density)) +
  geom_line(color = "darkgreen") +
  labs(title = "Updated Marginal of Y", x = "Y", y = "Density") +
  theme_minimal()
)
Rhat <- result$R
Fx$params <- result$params1
Gy$params <- result$params2

xg <- seq(0, 100, length.out = 50)
yg <- seq(0, 100, length.out = 50)
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
