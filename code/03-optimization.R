###############################################################################
# MAT5182 Computational Statistics — Assignment 3: Optimization
# Yann-Emmanuel Obrou, Winter 2025
#
# Reconstructs the R code from the submitted report
# (reports/assignement-3.docx) after the original .Rmd file was lost.
# Questions 1 (graduate) was a hand-written proof; the code below covers
# Questions 2-5 exactly as submitted.
#
# Required packages: none (base R only).
###############################################################################

## ----------------------------------------------------------------------------
## Question 1 — EM expressions (graduate students only)
## ----------------------------------------------------------------------------
# Theoretical derivation done by hand (paper); no R code was submitted.

## ----------------------------------------------------------------------------
## Question 2 — Maximizing f(x) = log(x) / (1 + x)
## ----------------------------------------------------------------------------
# A curve plot over [0, 10] shows the maximum lies between 1 and 10,
# so the grid search uses a = 1, b = 10.

curve(log(x) / (1 + x), from = 0, to = 10)

# --- Grid search ---
gridS <- function(f, n, a, b) {
  x_grid <- a + ((0:n) / n) * (b - a)   # n+1 points between a and b
  res <- rep(-1, n + 1)
  for (i in 1:(n + 1)) {
    res[i] <- f(x_grid[i])
  }
  return(list(M = max(res), argmax = x_grid[which.max(res)]))
}

# f(x) = ln(x) / (1 + x)
f <- function(x) {
  return(log(x) / (1 + x))
}

grid_result <- gridS(f, 10000, 1, 10)
cat("Grid Search:\n")
cat("Maximum found at x =", grid_result$argmax,
    "with f(x) =", grid_result$M, "\n")
# Report: Maximum found at x = 3.5911 with f(x) = 0.2784645

# --- Newton-Raphson ---

# First derivative of f
fprime <- function(x) {
  return(((1 + x) / x - log(x)) / (1 + x)^2)
}

# Second derivative of f
f2prime <- function(x) {
  A <- (1 + x) / x - log(x)
  Aprime <- -1 / (x^2) - 1 / x
  return((Aprime * (1 + x)^2 - 2 * A * (1 + x)) / (1 + x)^4)
}

NR <- function(x0, tol = 1e-6, maxit = 1000) {
  x <- x0
  for (i in 1:maxit) {
    x_new <- x - fprime(x) / f2prime(x)
    if (abs(x_new - x) < tol) break
    x <- x_new
  }
  return(x)
}

x_nr <- NR(3)   # initial guess x0 = 3
cat("Newton-Raphson:\n")
cat("Maximum found at x =", x_nr, "with f(x) =", f(x_nr), "\n\n")
# Report: Maximum found at x = 3.591121 with f(x) = 0.2784645

# --- Gradient ascent ---

# Simple gradient ascent: builds the full path of iterates
SimpleGrad <- function(grad, n, gamma, x0) {
  res <- rep(x0, n + 1)
  for (i in 1:n) {
    res[i + 1] <- res[i] + gamma[i] * grad(res[i])
  }
  return(res)
}

n_iter <- 100000
gamma_vec <- rep(0.01, n_iter)
x0 <- 3
x_values <- SimpleGrad(fprime, n_iter, gamma_vec, x0)
x_ga <- tail(x_values, 1)   # last iterate is the maximizer
cat("Gradient Ascent:\n")
cat("Maximum found at x =", x_ga, "with f(x) =", f(x_ga), "\n")
# Report: Maximum found at x = 3.591121 with f(x) = 0.2784645

## ----------------------------------------------------------------------------
## Question 3 — MLE of the Pareto shape parameter: analytical vs Newton-Raphson
## ----------------------------------------------------------------------------

# Pareto(alpha) sampler via the inverse transform (CDF inversion)
rpareto <- function(n, alpha, x_min = 1) {
  U <- runif(n)
  return(x_min / (1 - U)^(1 / alpha))
}

# Analytical MLE: alpha_hat = n / sum(log(x_i))
analytical_MLE <- function(x) {
  n <- length(x)
  return(n / sum(log(x)))
}

# Newton-Raphson for the root of l'(alpha) = 0,
# with l'(theta) = n/theta - L, l''(theta) = -n/theta^2, L = sum(log(x))
NR <- function(x, theta0, tol = 1e-6, maxit = 1000) {
  theta <- theta0
  n <- length(x)
  L <- sum(log(x))
  for (i in 1:maxit) {
    lprime <- n / theta - L
    lprime2 <- -n / theta^2
    theta_new <- theta - lprime / lprime2
    if (abs(theta_new - theta) < tol) break
    theta <- theta_new
  }
  return(theta)
}

set.seed(123)

alpha_true <- 2
n <- 1000     # sample size
m <- 100      # number of replications

analytical_estimates <- numeric(m)
NR_estimates <- numeric(m)

for (i in 1:m) {
  x <- rpareto(n, alpha_true)
  analytical_estimates[i] <- analytical_MLE(x)
  NR_estimates[i] <- NR(x, theta0 = 1)   # initial guess theta0 = 1
}

cat("Analytical MLE: Mean =", mean(analytical_estimates),
    "SD =", sd(analytical_estimates), "\n")
# Report: Analytical MLE: Mean = 2.006097 SD = 0.05570272
cat("Newton-Raphson MLE: Mean =", mean(NR_estimates),
    "SD =", sd(NR_estimates), "\n")
# Report: Newton-Raphson MLE: Mean = 2.006097 SD = 0.05570272

# Conclusion (report): both methods agree; the means are close to the true
# alpha = 2 and the low SD shows the estimates are stable across samples.

## ----------------------------------------------------------------------------
## Question 4 — Two-component Gaussian mixture fitted by EM
## ----------------------------------------------------------------------------
# p(y_i) = q * N(y_i | mu1, sigma1^2) + (1 - q) * N(y_i | mu2, sigma2^2)
# Latent component labels z_i are unobserved -> EM algorithm.

EM_gauss_mix <- function(y, theta0, tol = 1e-6, maxit = 1000) {
  # theta0 = (q, mu1, mu2, sigma1, sigma2)
  theta <- theta0
  iter <- 0
  n <- length(y)

  while (iter < maxit) {
    q      <- theta[1]
    mu1    <- theta[2]
    mu2    <- theta[3]
    sigma1 <- theta[4]
    sigma2 <- theta[5]

    # E-step: responsibilities (log-densities for numerical stability)
    log_tau1 <- log(q) + dnorm(y, mean = mu1, sd = sigma1, log = TRUE)
    log_tau2 <- log(1 - q) + dnorm(y, mean = mu2, sd = sigma2, log = TRUE)
    gamma1 <- exp(log_tau1) / (exp(log_tau1) + exp(log_tau2))
    gamma2 <- 1 - gamma1

    # M-step: weighted MLE updates
    q_new      <- mean(gamma1)
    mu1_new    <- sum(gamma1 * y) / sum(gamma1)
    mu2_new    <- sum(gamma2 * y) / sum(gamma2)
    sigma1_new <- sqrt(sum(gamma1 * (y - mu1_new)^2) / sum(gamma1))
    sigma2_new <- sqrt(sum(gamma2 * (y - mu2_new)^2) / sum(gamma2))

    theta_new <- c(q_new, mu1_new, mu2_new, sigma1_new, sigma2_new)

    if (max(abs(theta_new - theta)) < tol) break

    theta <- theta_new
    iter <- iter + 1
  }

  return(list(theta_hat = theta, iterations = iter))
}

# Dataset from Houtot et al. 2002 (Geophysics)
y_data <- c(8.1, 8.2, 8.1, 8.2, 8.2, 7.4, 7.3, 7.4, 8.1, 8.1,
            7.9, 7.8, 8.2, 7.9, 7.9, 8.1, 8.1)

# Three different initializations: theta0 = (q, mu1, mu2, sigma1, sigma2)
init1 <- c(0.5, 7.5, 8.5, 0.5, 0.5)
init2 <- c(0.3, 7.8, 8.0, 0.3, 0.3)
init3 <- c(0.7, 7.2, 8.4, 0.7, 0.7)

result1 <- EM_gauss_mix(y_data, init1)
result2 <- EM_gauss_mix(y_data, init2)
result3 <- EM_gauss_mix(y_data, init3)

cat("EM Algorithm Results:\n")
cat("Initialization 1:\n")
print(result1$theta_hat)
# Report: 0.44742083 7.69530202 8.14025989 0.27766802 0.04905319
cat("Number of iterations:", result1$iterations, "\n\n")
# Report: 18
cat("Initialization 2:\n")
print(result2$theta_hat)
# Report: 0.17647015 7.36666659 8.06428536 0.04714048 0.12877056
cat("Number of iterations:", result2$iterations, "\n\n")
# Report: 16
cat("Initialization 3:\n")
print(result3$theta_hat)
# Report: 0.44742325 7.69530445 8.14025988 0.27766925 0.04905319
cat("Number of iterations:", result3$iterations, "\n\n")
# Report: 14

# Conclusion (report): EM converges, but to different local maxima
# depending on the starting point (inits 1 and 3 agree, init 2 differs),
# so multiple initializations should always be tried.

## ----------------------------------------------------------------------------
## Question 5 — Multinomial model with a single parameter theta
## ----------------------------------------------------------------------------
# pi1 = 1/2 + theta/4, pi2 = 1/4 - theta/4, pi3 = 1/4 - theta/4,
# pi4 = 1 - (pi1 + pi2 + pi3) = theta/4,  theta in (0, 1).
# Data: n = 197, x = (125, 18, 20, 34).

x1 <- 125; x2 <- 18; x3 <- 20; x4 <- 34
n <- x1 + x2 + x3 + x4   # 197

pi1 <- function(th) 0.5 + th / 4
pi2 <- function(th) 0.25 - th / 4
pi3 <- function(th) 0.25 - th / 4
pi4 <- function(th) 1 - (pi1(th) + pi2(th) + pi3(th))

# Log-likelihood
ell <- function(th) {
  return(x1 * log(pi1(th)) + x2 * log(pi2(th)) +
         x3 * log(pi3(th)) + x4 * log(pi4(th)))
}

# Derivative of the log-likelihood (note: 38 = x2 + x3)
ell_prime <- function(th) {
  return((125 / (4 * (0.5 + th / 4))) -
         (38 / (4 * (0.25 - th / 4))) + (34 / th))
}

# Plot the log-likelihood and its derivative
theta_seq <- seq(0.001, 0.999, length.out = 500)
ell_values <- unlist(sapply(theta_seq, ell))
plot(theta_seq, ell_values, type = "l",
     main = "Log-likelihood", xlab = "theta", ylab = "ell(theta)")

ellp_values <- unlist(sapply(theta_seq, ell_prime))
plot(theta_seq, ellp_values, type = "l",
     main = "Derivative of log-likelihood", xlab = "theta", ylab = "ell'(theta)")

# The derivative crosses zero near theta = 0.626, inside (0, 1): the MLE.

# --- Newton's method, starting at theta0 = 0.5 ---

# Second derivative of the log-likelihood
ell_second <- function(th) {
  return((-125 / (16 * (0.5 + th / 4)^2)) +
         (-38 / (16 * (0.25 - th / 4)^2)) + (-34 / (th^2)))
}

NR <- function(theta0, tol = 1e-6, maxit = 1000) {
  th <- theta0
  for (i in 1:maxit) {
    d1 <- ell_prime(th)
    d2 <- ell_second(th)
    th_new <- th - d1 / d2
    if (abs(th_new - th) < tol) break
    th <- th_new
  }
  return(th)
}

theta_NR <- NR(0.5)
cat("Newton's method result:", theta_NR, "\n")
# Report: Newton's method result: 0.6268215

# --- EM algorithm (complete-data MLE exists in closed form) ---
# E-step/M-step update derived in the report (paper derivation).

EM_multinomial <- function(theta0, tol = 1e-6, maxit = 1000) {
  th <- theta0
  iter <- 0
  while (iter < maxit) {
    p <- (0.5) / (0.5 + th / 4)
    E_Y <- 125 * p

    # M-step: closed-form update
    th_new <- (68 + 159 * th) / (144 + 197 * th)

    if (abs(th_new - th) < tol) break
    th <- th_new
    iter <- iter + 1
  }
  return(list(theta_hat = th, iterations = iter))
}

theta_EM <- EM_multinomial(0.5)
cat("EM Algorithm:\n")
cat("Solution by EM (theta) =", theta_EM$theta_hat,
    "after", theta_EM$iterations, "iterations.\n")
# Report: Solution by EM (theta) = 0.6268207 after 6 iterations.

# --- Gradient ascent, starting at theta0 = 0.5 ---
# NOTE: a learning rate of 0.01 overshoots (the report obtained ~ -63.95);
# 0.0001 converges steadily to the MLE.

SimpleGrad <- function(grad, n, gamma, x0) {
  res <- rep(x0, n + 1)
  for (i in 1:n) {
    res[i + 1] <- res[i] + gamma[i] * grad(res[i])
  }
  return(res)
}

n_iter <- 1000
gamma_vec <- rep(0.0001, n_iter)
theta_values <- SimpleGrad(ell_prime, n_iter, gamma_vec, 0.5)
theta_GA <- tail(theta_values, 1)
cat("Gradient Ascent:\n")
cat("Solution by Gradient Ascent (theta) =", theta_GA, "\n")
# Report: Solution by Gradient Ascent (theta) = 0.6268215
