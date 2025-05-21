library(copulaData)
library(BernsteinCopula)

x <- NELS88$Math
y <- NELS88$Reading
u <- pnorm(x, mean=mean(x), sd=sd(x))
v <- pnorm(y, mean=mean(y), sd=sd(y))
# ビンの境界を指定（[0,1] を10分割）
u_breaks <- seq(0, 1, length.out = 11)
v_breaks <- seq(0, 1, length.out = 11)
# 2次元ヒストグラムを作成
counts <- table(
  cut(u, breaks = u_breaks, include.lowest = TRUE, right = FALSE),
  cut(v, breaks = v_breaks, include.lowest = TRUE, right = FALSE)
)
counts <- as.matrix(counts)

image(
  x = u_breaks[-length(u_breaks)],
  y = v_breaks[-length(v_breaks)],
  z = t(counts),  # 転置が必要
  col = heat.colors(20),
  xlab = "x", ylab = "y", main = "2D Binned Data"
)

m <- 4
n <- 4
R_init <- matrix(1 / (m*n), m, n)
result <- fit.copula.group(u_breaks, v_breaks, counts, R_init)

