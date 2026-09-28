###############################################################################
# MAT5182 Computational Statistics — Assignment 1: Simulation
# Yann-Emmanuel Obrou, Winter 2025
#
# This script reconstructs the R code from the submitted report
# (reports/assignement-1-computational-stat.pdf) after the original
# .Rmd file was lost. It is written to reproduce the analyses, figures
# and printed results of that report as closely as possible.
#
# Required packages: none (base R only).
#
# NOTE ON SEEDS: the original report only called set.seed() in Question 4.
# A seed is set at the top here so that the plots in Questions 1-3 are
# reproducible; it does not affect Question 4 (which re-seeds itself).
###############################################################################

set.seed(2025)

## ----------------------------------------------------------------------------
## Question 1 — Simulating from an arbitrary Beta distribution
## ----------------------------------------------------------------------------
# If Y ~ Gamma(a,1), Z ~ Gamma(b,1) independent, then Y/(Y+Z) ~ Beta(a,b).
# Sample of size n = 1000 from Beta(1/2, 1) via the Gamma-ratio method,
# then (as in the report) a direct rbeta() sample, which is what feeds
# the plots below.

n <- 1000
a <- 0.5
b <- 1

# Simulate Y ~ Gamma(a, 1) and Z ~ Gamma(b, 1)
Y <- rgamma(n, shape = a, rate = 1)
Z <- rgamma(n, shape = b, rate = 1)

# X = Y / (Y + Z)  ~  Beta(a, b)
X <- Y / (Y + Z)

# Simulate X ~ Beta(0.5, 1) directly (this draw is used for the figures)
X <- rbeta(n, shape1 = 0.5, shape2 = 1)

# Sample of size 1000 with Y = U^2, U standard uniform
n <- 1000
U <- runif(n, 0, 1)
Y <- U^2

# Histograms for X and Y
par(mfrow = c(2, 1))
hist(X, breaks = 40, main = "Histogram of X ~ Beta(0.5, 1)",
     xlab = "X", col = "lightblue")
hist(Y, breaks = 40, main = "Histogram of Y = U^2",
     xlab = "Y", col = "pink")
par(mfrow = c(1, 1))

# The two histograms look the same: U^2 has the Beta(1/2, 1) distribution.

# QQ-plot of X vs Y
qqplot(X, Y, main = "QQ Plot of X vs. Y",
       xlab = "Quantiles of X", ylab = "Quantiles of Y")
abline(0, 1, col = "red")
# The points follow the line y = x: U^2 is equivalent in distribution
# to Beta(1/2, 1).

## ----------------------------------------------------------------------------
## Question 2 — Simulating from a Pareto distribution
## ----------------------------------------------------------------------------
# Pareto(alpha) density: f(x) = alpha * x^(-alpha-1), x > 1.
# Inverse transform: X = (1 - U)^(-1/alpha), U ~ Uniform(0,1).

alpha <- 2
n <- 1000

u <- runif(n, 0, 1)
X <- (1 - u)^(-1 / alpha)

# Histogram with the theoretical density overlaid
hist(X, freq = FALSE)
curve(alpha * x^(-alpha - 1), xlim = c(0, 10), add = TRUE, col = "blue")

# QQ-plot against the theoretical Pareto quantiles
orderX <- sort(X)
n <- length(X)
grid <- (1:n) / (n + 1)
quantiles <- (1 - grid)^(-1 / alpha)
plot(quantiles, orderX,
     main = "QQ-plot: simulated vs theoretical Pareto quantiles",
     xlab = "Theoretical quantiles", ylab = "Sample quantiles")
# The histogram follows the theoretical density and the QQ-plot follows
# the y = x line: the inverse-transform sample is Pareto(alpha = 2).

## ----------------------------------------------------------------------------
## Question 3 — Gamma with non-integer shape by rejection sampling
## ----------------------------------------------------------------------------
# Target: X ~ Gamma(3/2, 1). Proposal g: Exponential with the same mean
# as X (mean = 3/2, so rate lambda = 2/3). The envelope constant
# alpha = 1.5 was chosen manually in the report for a high acceptance rate.

n <- 1000
alpha <- 1.5   # envelope constant, manually adjusted (report)

# Gamma(3/2, 1) density (target)
f <- function(x) {
  (x^(3/2 - 1) * exp(-x)) / gamma(3/2)
}

# Exponential proposal density with the same mean as the target
lambda <- 2 / 3
g <- function(x) {
  lambda * exp(-lambda * x)
}

# Rejection sampling
samples <- numeric(n)
count <- 0
while (count < n) {
  Y <- rexp(1, rate = lambda)   # sample from g
  U <- runif(1)                 # uniform draw
  if (U <= f(Y) / (alpha * g(Y))) {
    count <- count + 1
    samples[count] <- Y
  }
}

# Comparison with R's built-in rgamma()
rgamma_samples <- rgamma(n, shape = 3/2, rate = 1)

par(mfrow = c(1, 2))
hist(samples, prob = TRUE, main = "Rejection Sampling", xlim = c(0, 10))
curve(f, add = TRUE, col = "red")
hist(rgamma_samples, prob = TRUE, main = "rgamma function", xlim = c(0, 10))
curve(f, add = TRUE, col = "blue")
par(mfrow = c(1, 1))

## ----------------------------------------------------------------------------
## Question 4 — Estimating a censored exponential mean (consulting problem)
## ----------------------------------------------------------------------------
# Clients observe laptop lifetimes only when the laptop breaks within
# 4 years. They propose lambda_hat' = (sum of observed X_i) / (sum of B_i),
# where B_i = 1 if X_i <= 4. We simulate the experiment to study the bias.

set.seed(123)  # for reproducibility (as in the report)

lambda_values <- c(1/4, 1/5, 1/3)  # true rate parameters (means 4, 5, 3 yrs)
n <- 200                           # laptops per experiment
num_simulations <- 200             # simulations per lambda value

# Function simulating the experiment for a given lambda
simulate_experiment <- function(lambda, n, time_limit, num_simulations = 1000) {
  estimators <- numeric(num_simulations)
  for (sim in 1:num_simulations) {
    # Generate exponential lifetimes
    u <- runif(n)
    EXPO <- -log(u) / lambda
    # Indicator: laptop failed within the study period
    B <- ifelse(EXPO <= time_limit, 1, 0)
    # Estimate lambda using only observed failures (proposed estimator)
    observed_failures <- EXPO[B == 1]
    lambda_hat_prime <- mean(observed_failures) / mean(B)
    estimators[sim] <- lambda_hat_prime
  }
  return(estimators)
}

# Run the experiment for each lambda value
results <- lapply(lambda_values,
                  function(l) simulate_experiment(l, n, time_limit = 4,
                                                  num_simulations))

# Distribution of the estimator for each lambda
par(mfrow = c(1, length(lambda_values)))
for (i in 1:length(lambda_values)) {
  hist(results[[i]], breaks = 30,
       main = paste("Estimator Distribution for lambda =", lambda_values[i]),
       xlab = expression(hat(lambda) ~ "Estimator"),
       col = "lightblue", border = "black")
  abline(v = lambda_values[i], col = "red", lwd = 2, lty = 2)
  legend("topright", legend = paste("True lambda =", lambda_values[i]),
         col = "red", lwd = 2, lty = 2, bty = "n")
}
par(mfrow = c(1, 1))

# Bias of the estimator (computed as in the report)
biases <- sapply(results, function(est) mean(est) - mean(lambda_values))
bias_table <- data.frame(Lambda = lambda_values, Bias = biases)
print(bias_table)
# Report output:
#      Lambda     Bias
# 1 0.2500000 2.428617
# 2 0.2000000 2.899972
# 3 0.3333333 1.857672

# Conclusion (report): the proposed estimator is biased because it only
# uses lifetimes observed within 4 years, ignoring longer-lived laptops.
# Following participants for more than 4 years would reduce the bias.

## ----------------------------------------------------------------------------
## Optional Question 1 (chosen) — Simulating from Wigner's semicircle law
## ----------------------------------------------------------------------------
# Wigner (semicircle) density on [-1, 1]: f(x) = sqrt(1 - x^2) / (2*pi).
# Rejection sampling with a Uniform(-1, 1) proposal, g(x) = 1/2,
# and envelope constant c = 1/pi.

n <- 1000
c <- 1 / pi

# Wigner density
f <- function(x) {
  (1 / (2 * pi) * sqrt(1 - x^2)) * (abs(x) <= 1)
}

# Uniform(-1, 1) proposal density
g <- function(x) {
  0.5 * (abs(x) <= 1)
}

# Rejection sampling
samples1 <- numeric(n)
count <- 0
while (count < n) {
  Y <- runif(1, -1, 1)   # sample from g
  U <- runif(1)
  if (U <= f(Y) / (c * g(Y))) {
    count <- count + 1
    samples1[count] <- Y
  }
}

# Histogram with the theoretical density overlaid
hist(samples1, breaks = 40, probability = TRUE,
     main = "Histogram of Wigner Law Samples", xlab = "x")
curve(f, from = -1, to = 1, add = TRUE, col = "red", lwd = 2)
