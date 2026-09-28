###############################################################################
# MAT5182 Computational Statistics — Assignment 4: Bootstrap
# Yann-Emmanuel Obrou, Winter 2025
#
# Reconstructs the R code from the submitted report
# (reports/assignement-4-bootstrap.docx) after the original .Rmd file
# was lost. RNG seeds and draw order match the report.
#
# Required packages: none (base R only).
###############################################################################

## ----------------------------------------------------------------------------
## Question 1 — Bootstrap confidence intervals for the median of |X|, X ~ N(0,1)
## ----------------------------------------------------------------------------
# The true median of |X| is Phi^{-1}(0.75).

true_theta <- qnorm(0.75)
cat("the true theta is:", true_theta, "\n")
# Report: the true theta is: 0.6744898

set.seed(123)
n <- 1000

# Sample X ~ N(0,1), set Y = |X|
X <- rnorm(n, mean = 0, sd = 1)
Y <- abs(X)

# Crude Monte Carlo estimate: sample median of Y
theta_hat <- median(Y)
cat("Crude MC estimated median of |X|:", theta_hat, "\n")
# Report: Crude MC estimated median of |X|: 0.6486526

# --- Bootstrap confidence intervals ---
B <- 1000
boot_medians <- numeric(B)

# Resample Y with replacement; store the median of each replicate
for (i in 1:B) {
  boot_sample <- sample(Y, n, replace = TRUE)
  boot_medians[i] <- median(boot_sample)
}

# Bootstrap standard error
se_boot <- sd(boot_medians)

# Normal-approximation CI: theta_hat +/- z_{1-alpha/2} * se_boot
alpha <- 0.05
z_alpha <- qnorm(1 - alpha / 2)
IC_normal <- c(theta_hat - z_alpha * se_boot, theta_hat + z_alpha * se_boot)
cat("Normal bootstrap CI:", IC_normal, "\n")
# Report: Normal bootstrap CI: 0.6009459 0.6963594

# Percentile bootstrap CI
IC_percentile <- quantile(boot_medians, probs = c(alpha / 2, 1 - alpha / 2))
cat("Percentile bootstrap CI:", IC_percentile, "\n")
# Report: Percentile bootstrap CI: 0.6030025 0.7024453

# Justification (report): n = 1000 is large (CLT for the median), and the
# normal CI (0.6009, 0.6964) agrees closely with the percentile CI
# (0.6030, 0.7024), so the normal approximation is reasonable here.

## ----------------------------------------------------------------------------
## Question 2 — Estimating the endpoint theta of Uniform(0, theta)
## ----------------------------------------------------------------------------
# The MLE is theta_hat = max(X_1, ..., X_n); it is biased downward.
# (Theoretical arguments were done by hand; code below.)

set.seed(1234)
theta_true <- 10
n <- 100
M <- 1000

# (a) Monte Carlo distribution of the MLE
theta_estimates <- numeric(M)
for (i in 1:M) {
  U <- runif(n, min = 0, max = theta_true)
  theta_estimates[i] <- max(U)
}

theta_MC_mean <- mean(theta_estimates)
bias_MC <- theta_MC_mean - theta_true
cat("Monte Carlo (Part a) estimate of theta: =", theta_MC_mean,
    "and Bias =", bias_MC, "\n")
# Report: estimate = 9.901523, Bias = -0.09847736

hist(theta_estimates, col = "lightblue",
     main = "MC: Distribution of MLE", xlab = expression(hat(theta)))

# One reference sample for the bootstrap/subsampling parts
set.seed(5678)
U_sample <- runif(n, min = 0, max = theta_true)
theta_hat <- max(U_sample)

# (b) "Parametric" bootstrap (B = 1000)
# NOTE (faithful to the report): the code below actually resamples the
# original sample with replacement instead of drawing from the fitted
# Uniform(0, theta_hat) model.
B <- 1000
boot_thetas <- numeric(B)
for (i in 1:B) {
  boot_sample <- sample(U_sample, n, replace = TRUE)
  boot_thetas[i] <- max(boot_sample)
}
parametric_estimate <- mean(boot_thetas)
cat("Parametric Bootstrap (Part b) estimate of theta:", parametric_estimate, "\n")
# Report: Parametric Bootstrap (Part b) estimate of theta: 9.900374

hist(boot_thetas, col = "lightgreen",
     main = "Bootstrap: Distribution of MLE", xlab = expression(hat(theta)^"*"))

# (c) Subsampling with m = floor(sqrt(n)) = 10 (B = 1000)
m <- floor(sqrt(n))
B <- 1000
subsample_thetas <- numeric(B)
for (i in 1:B) {
  subsample <- sample(U_sample, m, replace = FALSE)
  subsample_thetas[i] <- max(subsample)
}
subsampling_estimate <- mean(subsample_thetas)
cat("Subsampling (Part c) estimate of theta:", subsampling_estimate, "\n")
# Report: Subsampling (Part c) estimate of theta: 9.204225

hist(subsample_thetas, col = "lightcoral",
     main = "Subsampling: Distribution of MLE", xlab = "Max (subsample)")

# (d) Nonparametric bootstrap (B = 1000)
B <- 1000
boot_thetas_np <- numeric(B)
for (i in 1:B) {
  boot_sample_np <- sample(U_sample, n, replace = TRUE)
  boot_thetas_np[i] <- max(boot_sample_np)
}
nonparametric_estimate <- mean(boot_thetas_np)
cat("Nonparametric Bootstrap (Part d) estimate of theta:",
    nonparametric_estimate, "\n")
# Report: Nonparametric Bootstrap (Part d) estimate of theta: 9.902998

hist(boot_thetas_np, col = "lightyellow",
     main = "Nonparametric Bootstrap: MLE", xlab = expression(hat(theta)^"*"))

# (e) Compare the four histograms
par(mfrow = c(2, 2))
hist(theta_estimates, col = "lightblue", main = "MC: MLE",
     xlab = expression(hat(theta)), breaks = 20)
hist(boot_thetas, col = "lightgreen", main = "Parametric Bootstrap: MLE",
     xlab = expression(hat(theta)^"*"), breaks = 20)
hist(subsample_thetas, col = "lightcoral", main = "Subsampling: MLE",
     xlab = "Max (subsample, m=10)", breaks = 20)
hist(boot_thetas_np, col = "lightyellow", main = "Nonparametric Bootstrap: MLE",
     xlab = expression(hat(theta)^"*"), breaks = 20)
par(mfrow = c(1, 1))

# Conclusion (report): MC (9.9015), "parametric" (9.9004) and nonparametric
# (9.9030) bootstrap agree; subsampling with m = 10 is much lower (9.2042)
# because smaller subsamples produce smaller maxima.

## ----------------------------------------------------------------------------
## Question 3 — Bootstrap for kurtosis of Exp(1): parametric works,
##              nonparametric struggles
## ----------------------------------------------------------------------------
# Empirical kurtosis: theta_hat = m4 / m2^2, with m2, m4 the 2nd and 4th
# central sample moments.

set.seed(4321)
n <- 1000
B <- 1000

calc_kurtosis <- function(x) {
  xbar <- mean(x)
  m2 <- mean((x - xbar)^2)
  m4 <- mean((x - xbar)^4)
  return(m4 / m2^2)
}

# (a) Empirical distribution: B samples from Exp(1), kurtosis of each
kurtosis_empirical <- numeric(B)
for (i in 1:B) {
  sample_exp <- rexp(n, rate = 1)
  kurtosis_empirical[i] <- calc_kurtosis(sample_exp)
}

cat("Empirical(part a) mean kurtosis:", mean(kurtosis_empirical), "\n")
# Report: Empirical(part a) mean kurtosis: 8.772293

hist(kurtosis_empirical, col = "lightblue",
     main = "Empirical Distribution of Kurtosis", xlab = "Kurtosis")

# (b) Parametric bootstrap: one sample from Exp(1), estimate the rate,
# then B samples from Exp(lambda_hat)
set.seed(9876)
X <- rexp(n, rate = 1)
lambda_hat <- 1 / mean(X)

kurtosis_param_boot <- numeric(B)
for (i in 1:B) {
  X_star <- rexp(n, rate = lambda_hat)
  kurtosis_param_boot[i] <- calc_kurtosis(X_star)
}

cat("Parametric Bootstrap (Part b) mean kurtosis:",
    mean(kurtosis_param_boot), "\n")
# Report: Parametric Bootstrap (Part b) mean kurtosis: 8.804002

hist(kurtosis_param_boot, col = "lightgreen",
     main = "Parametric Bootstrap: Kurtosis", xlab = "Kurtosis")

# (c) Nonparametric bootstrap: resample the original sample X
kurtosis_nonparam_boot <- numeric(B)
for (i in 1:B) {
  X_star_np <- sample(X, n, replace = TRUE)
  kurtosis_nonparam_boot[i] <- calc_kurtosis(X_star_np)
}

cat("Nonparametric Bootstrap (Part c) mean kurtosis:",
    mean(kurtosis_nonparam_boot), "\n")
# Report: Nonparametric Bootstrap (Part c) mean kurtosis: 7.968489

hist(kurtosis_nonparam_boot, col = "lightcoral",
     main = "Nonparametric Bootstrap: Kurtosis", xlab = "Kurtosis")

# (d) The three histograms side by side
par(mfrow = c(1, 3))
hist(kurtosis_empirical, col = "lightblue",
     main = "Empirical Kurtosis", xlab = "Kurtosis")
hist(kurtosis_param_boot, col = "lightgreen",
     main = "Parametric Bootstrap", xlab = "Kurtosis")
hist(kurtosis_nonparam_boot, col = "lightcoral",
     main = "Nonparametric Bootstrap", xlab = "Kurtosis")
par(mfrow = c(1, 1))

# Conclusion (report): the parametric bootstrap reproduces the tail
# behaviour of the exponential and matches the empirical distribution
# (means ~8.77 vs ~8.80), while the nonparametric bootstrap is noisier
# (mean ~7.97, wider spread) because kurtosis is sensitive to extremes
# that the finite observed sample may under-represent.
