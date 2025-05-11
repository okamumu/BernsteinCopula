# BernsteinCopula

BernsteinCopula is an R package for statistical modeling and copula functions. It provides tools for estimating Bernstein copulas using the Expectation-Maximization (EM) algorithm and includes utilities for working with marginal distributions, joint densities, and copula transformations.

## Features

- **Bernstein Copula Estimation**: Estimate copulas using the EM algorithm.
- **Marginal Distribution Utilities**: Tools for estimating and visualizing marginal distributions.
- **Joint Density Visualization**: Plot joint density estimates and contours.
- **Integration with Rcpp**: High-performance C++ implementations for computationally intensive tasks.

## Installation

To install the package, clone the repository and build it using R:

```R
# Install devtools if not already installed
install.packages("devtools")

devtools::install_github("okamumu/BernsteinCopula")
```

## Usage

### Simulating Data

Generate bivariate normal data:

```R
library(MASS)
set.seed(123)
Sigma <- matrix(c(1, 0.7, 0.7, 1), 2, 2)
data <- MASS::mvrnorm(n = 200, mu = c(0, 0), Sigma = Sigma)
x <- data[,1]
y <- data[,2]
```

### Estimating a Bernstein Copula

Transform data to uniform margins and estimate the copula:

```R
library(BernsteinCopula)
Fx <- normaldist(mean = mean(x), sd = sd(x))
Gy <- normaldist(mean = mean(y), sd = sd(y))

u <- Fx$cdf(x)
v <- Gy$cdf(y)

m <- 4
n <- 4
R0 <- matrix(1/(m * n), m, n)
result <- emloop(u, v, R_init = R0)
```

### Visualizing the Copula Density

Plot the estimated joint density:

```R
xg <- seq(min(x), max(x), length.out = 50)
yg <- seq(min(y), max(y), length.out = 50)
Z <- bernstein_joint_pdf_grid(xg, yg, result$R, Fx, Gy)

plot_contour(xg, yg, Z, title = "Estimated Joint Density", xlabel = "x", ylabel = "y")
```

## Documentation

Detailed documentation is available in the `man/` directory and the vignette `vignettes/BernsteinCopula.Rmd`. To view the vignette, use:

```R
vignette("BernsteinCopula")
```

## Contributing

Contributions are welcome! Please submit issues or pull requests on the GitHub repository.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.