###############################################################################
# MAT5182 Computational Statistics — Assignment 5: Advanced Monte Carlo
# Yann-Emmanuel Obrou, Winter 2025
#
# Reconstructs the R code from the submitted report
# (reports/Assignement-5.pdf) after the original .Rmd file was lost.
# Question 3 (posterior density) was a hand-written derivation; no R code
# was submitted for it.
#
# Required packages: none (base R only).
###############################################################################

## ----------------------------------------------------------------------------
## Question 1 — Finite Markov chain: stationary distribution
## ----------------------------------------------------------------------------
# Transition matrix on {1, 2, 3}. The stationary distribution is
# pi = (1/3, 1/3, 1/3) (shown by hand in the report).

# (b) 10th power of Q
Q <- matrix(c(0.12, 0.24, 0.64,
              0.34, 0.40, 0.26,
              0.54, 0.36, 0.10), nrow = 3, byrow = TRUE)

# Function computing the power of a matrix
matrix_power <- function(matrix, power) {
  result <- diag(nrow(matrix))   # start with the identity matrix
  for (i in 1:power) {
    result <- result %*% matrix   # multiply by Q at each iteration
  }
  return(result)
}

Q10 <- matrix_power(Q, 10)
print(Q10)
# Report:
#           [,1]      [,2]      [,3]
# [1,] 0.3335765 0.3333908 0.3330327
# [2,] 0.3333013 0.3333258 0.3333729
# [3,] 0.3331222 0.3332834 0.3335944
# As expected, every row is close to (1/3, 1/3, 1/3).

# (c) Simulate the chain and compare visit proportions with pi

# Function simulating the Markov chain
simulate_markov <- function(transition_matrix, initial_state, n) {
  states <- numeric(n)
  states[1] <- initial_state
  for (i in 2:n) {
    current_state <- states[i - 1]   # 1, 2 or 3: used directly as row index
    states[i] <- sample(1:3, size = 1,
                        prob = transition_matrix[current_state, ])
  }
  return(states)
}

set.seed(123)

# NOTE: the assignment asked for a sample of size 100, while the report's
# code table shows sample_size <- 100000. Both runs are reproduced here.
sample_size <- 100
chain_sample <- simulate_markov(Q, initial_state = 1, sample_size)
prop_table <- table(chain_sample) / sample_size
print(prop_table)
# Report (n = 100, quoted in the text): approximately (0.36, 0.39, 0.25)

sample_size <- 100000
chain_sample <- simulate_markov(Q, initial_state = 1, sample_size)
prop_table <- table(chain_sample) / sample_size
print(prop_table)
# Report:
# chain_sample
#         1         2         3
# 0.33260   0.33548   0.33192

# By the law of large numbers for (irreducible, aperiodic) Markov chains,
# the visit proportions converge almost surely to pi = (1/3, 1/3, 1/3).

## ----------------------------------------------------------------------------
## Question 2 — Metropolis-Hastings for a non-standard density
## ----------------------------------------------------------------------------
# Target (unnormalized) on (-2, 2): f*(x) = -exp(-x) * sin(x) + 1.
# Normalizing constant C = integral_{-2}^{2} f*(x) dx (unknown).

# (a) Plot the unnormalized density
f_unnormalized <- function(x) {
  -exp(-x) * sin(x) + 1
}

x_vals <- seq(-2, 2, length.out = 1000)
y_vals <- f_unnormalized(x_vals)

plot(x_vals, y_vals, type = "l", col = "blue", lwd = 2,
     main = "Unnormalized Density f*(x)", xlab = "x", ylab = "f*(x)")

# (b) Estimate C by Monte Carlo: C = 4 * E[f*(X)], X ~ Uniform(-2, 2)
set.seed(123)

N <- 100000
x_samples <- runif(N, min = -2, max = 2)
C_estimate <- 4 * mean(f_unnormalized(x_samples))   # 4 = length of (-2, 2)
print(C_estimate)
# Report: [1] 8.936518

# (c) Metropolis-Hastings with an exponential proposal kernel
# Proposal: p(y) = lambda * exp(-lambda * y), y > 0. It is made symmetric
# by assigning a random sign to each exponential step, so the proposal
# ratio cancels in the acceptance probability. Proposed values outside
# (-2, 2) are resampled.

f_target <- function(x) {
  -exp(-x) * sin(x) + 1
}

MH_exponential <- function(x, lambda) {
  repeat {
    # Exponential step with a random sign
    delta <- rexp(1, rate = lambda)
    sign <- sample(c(-1, 1), 1)
    y <- x + sign * delta
    # Accept only proposals inside the support (-2, 2)
    if (y > -2 && y < 2) break
  }
  # Acceptance probability (proposal symmetric -> ratio simplifies)
  a <- min(1, f_target(y) / f_target(x))
  if (runif(1) < a) {
    return(y)   # accept the move
  } else {
    return(x)   # reject and stay
  }
}

# Simulate a Metropolis-Hastings chain
simulate_MH_chain <- function(lambda, n_iter = 5000, start = 0) {
  chain <- numeric(n_iter)
  chain[1] <- start
  for (i in 2:n_iter) {
    chain[i] <- MH_exponential(chain[i - 1], lambda)
  }
  return(chain)
}

# Normalized density used to overlay the true density on the histograms
f_normalized <- function(x) {
  f_target(x) / C_estimate
}

lambda_values <- c(0.1, 1, 10)
par(mfrow = c(1, 2))

for (lambda in lambda_values) {
  set.seed(123)
  chain <- simulate_MH_chain(lambda, n_iter = 5000, start = 0)

  # Histogram of the chain (after burn-in) vs the true normalized density
  hist(chain[1001:5000], breaks = 30, probability = TRUE,
       main = paste("Histogram, lambda =", lambda),
       xlab = "x", col = "lightblue", border = "black")
  x_curve <- seq(-2, 2, length.out = 1000)
  lines(x_curve, f_normalized(x_curve), col = "red", lwd = 2)

  # Trace (path) plot of the chain
  plot(chain, type = "l", col = "blue",
       main = paste("Trace plot, lambda =", lambda),
       xlab = "Iteration", ylab = "x")
}
par(mfrow = c(1, 1))

# Commentary (report):
# - lambda = 0.1: large jumps (mean step ~ 10); the chain bounces across
#   (-2, 2) but many proposals fall outside the support and are rejected,
#   so the histogram is rougher.
# - lambda = 1: moderate steps; good balance of acceptance and exploration;
#   the histogram tracks the true density most smoothly.
# - lambda = 10: small steps (mean ~ 0.1); high acceptance but slow,
#   local exploration; the chain takes longer to traverse (-2, 2).

## ----------------------------------------------------------------------------
## Question 3 — Posterior density (theory)
## ----------------------------------------------------------------------------
# X | Lambda = lambda ~ Exp(lambda), prior Lambda ~ Exp(beta):
# the posterior density was derived by hand (paper); no R code submitted.
