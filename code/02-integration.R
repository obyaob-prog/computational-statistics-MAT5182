###############################################################################
# MAT5182 Computational Statistics — Assignment 2: Integration
# Yann-Emmanuel Obrou, Winter 2025
#
# Reconstructs the R code from the submitted report
# (reports/assignement-2-computational-stat.docx) after the original
# .Rmd file was lost. Order of the random draws is kept identical to the
# report so the printed numbers are reproduced.
#
# Required packages: ggplot2 (only used for the Asian-option plots).
#   install.packages("ggplot2")
#
# True value of the integral (Riemann sums): I = 0.6278151.
###############################################################################

## ----------------------------------------------------------------------------
## Question 1 — Monte Carlo integration of I = integral_0^1 x^2 * exp(x^2) dx
## ----------------------------------------------------------------------------
# Write I = E[g(U)], U ~ Uniform(0,1), with g(x) = x^2 * exp(x^2).

set.seed(123)

n_values <- c(1e3, 1e4, 1e5)   # sample sizes
n <- c(10^3, 10^4, 10^5)

g <- function(x) x^2 * exp(x^2)

# --- Crude Monte Carlo ---
mc_estimates <- sapply(n, function(n) mean(g(runif(n))))
mc_estimates
# Report: 0.6187158 0.6173007 0.6274209

# --- Control variate (h(x) = x^2) ---
cv_function <- function(x) x^2

# Estimate E[h(U)] with a large pilot sample
mu_cv <- mean(cv_function(runif(100000)))

# Optimal coefficient gamma_cv = -Cov(g, h) / Var(h), from pilot samples
# (named gamma_cv so it does not shadow base::gamma(), used in Question 3)
covariance <- cov(g(runif(10000)), cv_function(runif(10000)))
variance <- var(cv_function(runif(10000)))
gamma_cv <- -covariance / variance

cv_estimates <- sapply(n, function(n)
  mean(g(runif(n)) + gamma_cv * (cv_function(runif(n)) - mu_cv))
)
cv_estimates
# Report: 0.6646270 0.6356376 0.6272932

# --- Antithetic variates ---
antithetic_estimate <- function(n) {
  u <- runif(n / 2)          # half the draws; pair u with 1 - u
  (mean(g(u)) + mean(g(1 - u))) / 2
}
antithetic_estimates <- sapply(n, antithetic_estimate)
antithetic_estimates
# Report: 0.6444586 0.6226092 0.6265897

# --- Variance comparison (as computed in the report) ---
var_table <- data.frame(
  n,
  Crude = var(mc_estimates),
  ControlVar = var(cv_estimates),
  Antithetic = var(antithetic_estimates)
)
print(var_table)
# Report:
#       n        Crude   ControlVar   Antithetic
# 1 1e+03 3.003344e-05 0.0003839702 0.0001354233
# 2 1e+04 3.003344e-05 0.0003839702 0.0001354233
# 3 1e+05 3.003344e-05 0.0003839702 0.0001354233

# --- Beta(3,1) representation: I = (1/3) * E[exp(V^2)], V ~ Beta(3,1) ---

# Crude Monte Carlo with Beta samples
beta_crude_mc <- function(n) {
  V <- rbeta(n, 3, 1)
  estimates <- exp(V^2)
  I_hat <- mean(estimates) / 3
  I_var <- var(estimates) / (n * 9)
  return(c(I_hat, I_var))
}

# Control variate with Beta samples: h(V) = V^2, E[V^2] = 3*4/(4*5) = 0.6
beta_control_mc <- function(n) {
  V <- rbeta(n, 3, 1)
  f_V <- exp(V^2)
  h_V <- V^2
  beta_hat <- -cov(f_V, h_V) / var(h_V)
  adjusted_estimates <- f_V + beta_hat * (h_V - 0.6)
  I_hat <- mean(adjusted_estimates) / 3
  I_var <- var(adjusted_estimates) / (n * 9)
  return(c(I_hat, I_var))
}

# Antithetic variates with Beta samples (inverse transform: V = U^(1/3))
beta_antithetic_mc <- function(n) {
  n <- ceiling(n / 2) * 2   # ensure n is even
  U <- runif(n / 2)
  V1 <- U^(1/3)
  V2 <- (1 - U)^(1/3)
  estimates <- (exp(V1^2) + exp(V2^2)) / 2
  I_hat <- mean(estimates) / 3
  I_var <- var(estimates) / (n * 9)
  return(c(I_hat, I_var))
}

n_values <- c(1e3, 1e4, 1e5)
beta_crude_results <- t(sapply(n_values, beta_crude_mc))
beta_control_results <- t(sapply(n_values, beta_control_mc))
beta_antithetic_results <- t(sapply(n_values, beta_antithetic_mc))

results_estimators <- data.frame(
  "Sample Size" = n_values,
  "Crude MC (Beta)" = round(beta_crude_results[, 1], 6),
  "Antithetic MC (Beta)" = round(beta_antithetic_results[, 1], 6),
  "CV (Beta)" = round(beta_control_results[, 1], 6),
  "True Value" = 0.6278151
)
print(results_estimators)
# Report:
#   Sample.Size Crude.MC..Beta. Antithetic.MC..Beta. CV..Beta. True.Value
# 1       1e+03        0.619362             0.627835  0.628361   0.627815
# 2       1e+04        0.626626             0.627803  0.627811   0.627815
# 3       1e+05        0.627424             0.627813  0.627843   0.627815

## ----------------------------------------------------------------------------
## Question 2 — Covariance inequality (graduate students)
## ----------------------------------------------------------------------------
# Theoretical proof done by hand (paper); no R code was submitted.

## ----------------------------------------------------------------------------
## Question 3 — Weibull tail probabilities
## ----------------------------------------------------------------------------
# F(x) = 1 - exp(-(x/c)^beta). Estimate P(X > q) with
# q = c * (-log(p))^(1/beta), for p = 0.001 and p = 0.005.

set.seed(123)
n <- 100
c_weib <- 1
beta <- 1.5

# Weibull sample of size 100 via the quantile (inverse-CDF) method
U <- runif(n)
X <- c_weib * (-log(U))^(1 / beta)

# --- Crude Monte Carlo, beta > 1 ---
n_mc <- 10000
beta <- 1.5
p_values <- c(0.001, 0.005)

for (p in p_values) {
  q <- c_weib * (-log(p))^(1 / beta)
  X_mc <- c_weib * (-log(runif(n_mc)))^(1 / beta)
  indicator <- X_mc > q

  p_hat <- mean(indicator)
  var_hat <- var(indicator) / n_mc
  ci_lower <- p_hat - 1.96 * sqrt(var_hat)
  ci_upper <- p_hat + 1.96 * sqrt(var_hat)

  cat("p =", p, ": p_hat =", p_hat, "Variance =", var_hat,
      "CI = [", ci_lower, ",", ci_upper, "]\n")
}
# Report:
# p = 0.001 : p_hat = 8e-04 Variance = 7.994399e-08 CI = [ 0.0002458224 , 0.001354178 ]
# p = 0.005 : p_hat = 0.005 Variance = 4.975498e-07 CI = [ 0.003617471 , 0.006382529 ]

# --- Crude Monte Carlo, beta < 1 ---
beta <- 0.5

for (p in p_values) {
  q <- c_weib * (-log(p))^(1 / beta)
  X_mc <- c_weib * (-log(runif(n_mc)))^(1 / beta)
  indicator <- X_mc > q

  p_hat <- mean(indicator)
  var_hat <- var(indicator) / n_mc
  ci_lower <- p_hat - 1.96 * sqrt(var_hat)
  ci_upper <- p_hat + 1.96 * sqrt(var_hat)

  cat("p =", p, ": p_hat =", p_hat, "Variance =", var_hat,
      "CI = [", ci_lower, ",", ci_upper, "]\n")
}
# Report:
# p = 0.001 : p_hat = 0.0012 Variance = 1.19868e-07 CI = [ 0.0005214097 , 0.00187859 ]
# p = 0.005 : p_hat = 0.0043 Variance = 4.281938e-07 CI = [ 0.003017444 , 0.005582556 ]

# --- Importance sampling, beta > 2, with Gamma(beta, 1/c) proposal ---
beta <- 3   # beta > 2

for (p in p_values) {
  q <- c_weib * (-log(p))^(1 / beta)
  Y <- rgamma(n_mc, shape = beta, scale = c_weib)
  w <- (beta * gamma(beta) / c_weib) * exp(-(Y / c_weib)^beta + Y / c_weib)
  indicator <- Y > q

  p_hat_is <- mean(indicator * w)
  var_hat_is <- var(indicator * w) / n_mc
  ci_lower_is <- p_hat_is - 1.96 * sqrt(var_hat_is)
  ci_upper_is <- p_hat_is + 1.96 * sqrt(var_hat_is)

  # Crude MC for comparison
  X_mc <- c_weib * (-log(runif(n_mc)))^(1 / beta)
  p_hat_mc <- mean(X_mc > q)
  var_hat_mc <- var(X_mc > q) / n_mc
  rel_eff <- ifelse(var_hat_mc > 0, var_hat_mc / var_hat_is, NA)

  cat("p =", p, ": IS p_hat =", p_hat_is, "Var =", var_hat_is,
      "CI = [", ci_lower_is, ",", ci_upper_is, "], Rel Eff =", rel_eff, "\n")
}
# Report:
# p = 0.001 : IS p_hat = 0.001035038 Var = 2.021812e-09 CI = [ 0.0009469078 , 0.001123169 ], Rel Eff = 29.6615
# p = 0.005 : IS p_hat = 0.004734837 Var = 4.114014e-08 CI = [ 0.004337289 , 0.005132384 ], Rel Eff = 15.69853

# --- Importance sampling, beta < 1, with a log-normal proposal ---
# (For beta < 1 the Gamma proposal is no longer suitable, so a
#  LogNormal(mu = log(q), sigma = 1) proposal centred at q is used.)
beta <- 0.5   # beta < 1

for (p in p_values) {
  q <- c_weib * (-log(p))^(1 / beta)
  mu <- log(q)
  sigma <- 1
  Y <- rlnorm(n_mc, meanlog = mu, sdlog = sigma)
  f_X <- (beta / c_weib) * (Y / c_weib)^(beta - 1) * exp(-(Y / c_weib)^beta)  # Weibull pdf
  f_Y <- dlnorm(Y, meanlog = mu, sdlog = sigma)                # LogNormal pdf
  w <- f_X / f_Y
  indicator <- Y > q

  p_hat_is <- mean(indicator * w)
  var_hat_is <- var(indicator * w) / n_mc
  ci_lower_is <- p_hat_is - 1.96 * sqrt(var_hat_is)
  ci_upper_is <- p_hat_is + 1.96 * sqrt(var_hat_is)

  cat("p =", p, ": IS p_hat =", p_hat_is, "Var =", var_hat_is,
      "CI = [", ci_lower_is, ",", ci_upper_is, "]\n")
}
# Report:
# p = 0.001 : IS p_hat = 0.0009906009 Var = 3.692509e-10 CI = [ 0.0009529377 , 0.001028264 ]
# p = 0.005 : IS p_hat = 0.005054143 Var = 6.994995e-09 CI = [ 0.004890216 , 0.00521807 ]

## ----------------------------------------------------------------------------
## Question 4 — Asian options: B0 (max-based) and D0 (average-based)
## ----------------------------------------------------------------------------
# B0 = E[exp(-rT) * max(max_i S(t_i) - K, 0)]
# D0 = E[exp(-rT) * max(mean_i S(t_i) - K, 0)]
# Methods: crude MC, antithetic sampling, control variates with S_T,
# control variates with (S_T - K)+ (European call, Black-Scholes price).

library(ggplot2)
set.seed(123)

asian_option_estimator <- function(sigma, r, T, S0, K, N, n, m) {
  dt <- T / N

  # European call price (Black-Scholes) used as the control-variate mean
  d1 <- (log(S0 / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  european_price <- S0 * pnorm(d1) - K * exp(-r * T) * pnorm(d2)

  # Simulate geometric Brownian motion paths.
  # NOTE (as in the submitted report): antithetic = TRUE draws fresh
  # normals and negates them, rather than reusing the original path's
  # normals.
  simulate_gbm <- function(n_paths, antithetic = FALSE) {
    dW <- matrix(rnorm(n_paths * N, 0, sqrt(dt)), nrow = n_paths, ncol = N)
    if (antithetic) dW <- -dW
    W <- t(apply(dW, 1, cumsum))
    t_grid <- matrix(seq(dt, T, by = dt), nrow = n_paths, ncol = N, byrow = TRUE)
    S <- S0 * exp((r - sigma^2 / 2) * t_grid + sigma * W)
    return(S)
  }

  # One row per replication (m rows), one column per method
  results_B <- matrix(0, nrow = m, ncol = 4)
  results_D <- matrix(0, nrow = m, ncol = 4)
  colnames(results_B) <- colnames(results_D) <-
    c("Crude MC", "Antithetic", "CV_S", "CV_Option")

  for (k in 1:m) {
    S <- simulate_gbm(n, antithetic = FALSE)
    S_antithetic <- simulate_gbm(n, antithetic = TRUE)

    # B0: payoff on the path maximum
    max_S <- apply(S, 1, max)
    payoff_B <- exp(-r * T) * pmax(max_S - K, 0)
    max_S_anti <- apply(S_antithetic, 1, max)
    payoff_B_anti <- exp(-r * T) * pmax(max_S_anti - K, 0)
    B_antithetic <- (payoff_B + payoff_B_anti) / 2

    # D0: payoff on the path average
    avg_S <- rowMeans(S)
    payoff_D <- exp(-r * T) * pmax(avg_S - K, 0)
    avg_S_anti <- rowMeans(S_antithetic)
    payoff_D_anti <- exp(-r * T) * pmax(avg_S_anti - K, 0)
    D_antithetic <- (payoff_D + payoff_D_anti) / 2

    # Control variate with S_T (E[S_T] = S0 * exp(r*T))
    S_T <- S[, N]
    gamma_B <- -cov(payoff_B, S_T) / var(S_T)
    B_cv_S <- mean(payoff_B) + gamma_B * (mean(S_T) - S0 * exp(r * T))
    gamma_D <- -cov(payoff_D, S_T) / var(S_T)
    D_cv_S <- mean(payoff_D) + gamma_D * (mean(S_T) - S0 * exp(r * T))

    # Control variate with (S_T - K)+ (European call payoff)
    payoff_EU <- exp(-r * T) * pmax(S_T - K, 0)
    gamma_B_EU <- -cov(payoff_B, payoff_EU) / var(payoff_EU)
    B_cv_EU <- mean(payoff_B) +
      gamma_B_EU * (mean(payoff_EU) - european_price)
    gamma_D_EU <- -cov(payoff_D, payoff_EU) / var(payoff_EU)
    D_cv_EU <- mean(payoff_D) +
      gamma_D_EU * (mean(payoff_EU) - european_price)

    results_B[k, ] <- c(mean(payoff_B), mean(B_antithetic), B_cv_S, B_cv_EU)
    results_D[k, ] <- c(mean(payoff_D), mean(D_antithetic), D_cv_S, D_cv_EU)
  }

  return(list(B0 = results_B, D0 = results_D))
}

# Parameters from the assignment
sigma <- 1
r <- 0.1
T <- 1
S0 <- 1
K <- 1
N <- 100    # grid size
n <- 100    # MC paths per estimate
m <- 100    # number of estimates per method

results <- asian_option_estimator(sigma, r, T, S0, K, N, n, m)

# Plot the m estimates for B0
df_B <- data.frame(
  Estimate = c(results$B0),
  Method = factor(rep(colnames(results$B0), each = m)),
  Index = rep(1:m, times = 4)
)
ggplot(df_B, aes(x = Index, y = Estimate, color = Method)) +
  geom_line() +
  labs(title = "Estimates for B0 (Maximum-Based Asian Option)", y = "Estimate") +
  theme_minimal()

# Variance comparison for B0
var_B <- apply(results$B0, 2, var)
cat("Variances for B0:\nCrude MC:", var_B[1], "\nAntithetic:", var_B[2],
    "\nCV with S_T:", var_B[3], "\nCV with EU Option:", var_B[4], "\n")
# Report:
# Crude MC: 0.01732895 / Antithetic: 0.008139353 /
# CV with S_T: 0.004177641 / CV with EU Option: 0.004498833

# Plot the m estimates for D0
df_D <- data.frame(
  Estimate = c(results$D0),
  Method = factor(rep(colnames(results$D0), each = m)),
  Index = rep(1:m, times = 4)
)
ggplot(df_D, aes(x = Index, y = Estimate, color = Method)) +
  geom_line() +
  labs(title = "Estimates for D0 (Average-Based Asian Option)", y = "Estimate") +
  theme_minimal()

# Variance comparison for D0
var_D <- apply(results$D0, 2, var)
cat("\nVariances for D0:\nCrude MC:", var_D[1], "\nAntithetic:", var_D[2],
    "\nCV with S_T:", var_D[3], "\nCV with EU Option:", var_D[4], "\n")
# Report:
# Crude MC: 0.002471226 / Antithetic: 0.001216017 /
# CV with S_T: 0.0009418998 / CV with EU Option: 0.00103126

## ----------------------------------------------------------------------------
## Optional Question 1 (chosen) — Expected Shortfall of N(0,1) at p = 0.01
## ----------------------------------------------------------------------------
# Theoretical ES: ES_Z(p) = phi(Phi^{-1}(1 - p)) / p.

set.seed(123)
p <- 0.01
n <- 10000

# Sample from N(0,1)
Z <- rnorm(n)

# Theoretical VaR at p = 0.01 (99th percentile)
q <- qnorm(1 - p)

# Tail sample and crude MC estimate of ES
tail_Z <- Z[Z > q]
ES_hat <- mean(tail_Z)

# Standard error for the confidence interval
se <- sd(tail_Z) / sqrt(length(tail_Z))

# 95% confidence interval
CI_lower <- ES_hat - qnorm(0.975) * se
CI_upper <- ES_hat + qnorm(0.975) * se

# Theoretical ES
true_ES <- dnorm(qnorm(1 - p)) / p

cat("Estimated ES:", ES_hat, "\n")
# Report: Estimated ES: 2.672539
cat("95% Confidence Interval:", CI_lower, "to", CI_upper, "\n")
# Report: 95% Confidence Interval: 2.609354 to 2.735725
cat("Theoretical ES:", true_ES, "\n")
# Report: Theoretical ES: 2.665214

# Conclusion (report): with 10,000 samples the Monte Carlo estimate
# (2.672539) closely matches the theoretical value (2.665214).
