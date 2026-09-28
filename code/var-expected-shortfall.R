###############################################################################
# MAT5182 Computational Statistics — Project:
# Estimating VaR and Expected Shortfall Using Monte Carlo Methods
# Yann-Emmanuel Obrou, Winter 2025
#
# Reconstructs the R code from the appendix of the submitted report
# (reports/project-Yann-OBROU.pdf, graded 20/20) after the original
# .Rmd file was lost.
#
# Setup:  tail probability p = 0.01, n = 10,000 samples per model.
#   Model 1: Z ~ N(0, 1)            (loss variable)
#   Model 2: X ~ Exp(rate = 10)     (loss variable)
#
# Methods: crude Monte Carlo (two implementations for the normal model)
# and importance sampling, each with 95% confidence intervals, compared
# against the analytical VaR/ES benchmarks:
#   Normal:      VaR = Phi^{-1}(1-p) ~ 2.326348, ES = phi(q)/p ~ 2.665214
#   Exponential: VaR = -log(p)/10    ~ 0.460517, ES = VaR + 1/10 ~ 0.560517
#
# Required packages: none (base R only).
###############################################################################

p <- 0.01
n <- 10000

###############################################################################
# Part 1: Normal Model (N(0,1))
###############################################################################

set.seed(10)

# Analytical benchmark
true_ES_norm <- dnorm(qnorm(1 - p)) / p
cat("True ES =", true_ES_norm, "\n\n")
# Report: True ES = 2.665214

Z <- rnorm(n)
threshold <- qnorm(1 - p)   # theoretical VaR at p = 0.01 (~2.326348)

# --- Crude Monte Carlo: Complete Vector Approach ---
# A length-n vector of zeros; tail values are filled in by a loop.
MCvalues <- rep(0, n)
for (i in 1:n) {
  if (Z[i] > threshold) {
    MCvalues[i] <- Z[i]
  }
}

ES_MC <- mean(MCvalues) / p
var_MC <- var(MCvalues) / n
lower_CI <- ES_MC - 1.96 * sqrt(var_MC)
upper_CI <- ES_MC + 1.96 * sqrt(var_MC)

cat("Crude MC (normal, complete vector) : ES =", ES_MC, "\n")
# Report: Crude MC (normal, complete vector) : ES = 3.075687
cat("95% CI : [", lower_CI, ",", upper_CI, "] \n\n")
# Report: 95% CI : [ 3.070056 , 3.081318 ]

# --- Crude Monte Carlo: Direct Extraction ---
# Extract the tail values directly and average them.
Y <- Z[Z > threshold]
ES_MC2 <- mean(Y)
var_MC2 <- var(Y) / n
lower_CI2 <- ES_MC2 - 1.96 * sqrt(var_MC2)
upper_CI2 <- ES_MC2 + 1.96 * sqrt(var_MC2)

cat("Crude MC (normal, direct extraction) : ES =", ES_MC2, "\n")
# Report: Crude MC (normal, direct extraction) : ES = 2.67451
cat("95% CI : [", lower_CI2, ",", upper_CI2, "] \n\n")
# Report: 95% CI : [ 2.668107 , 2.680914 ]

# --- Importance Sampling for the Normal Model ---
# Proposal: N(2.33, 1), i.e. a normal shifted so its mean is ~ the VaR,
# putting many more samples in the tail. Weights w(z) = phi(z)/phi(z; 2.33, 1).
set.seed(10)

mu_IS <- 2.33
Y_IS <- rnorm(n, mean = mu_IS, sd = 1)
w <- dnorm(Y_IS) / dnorm(Y_IS, mean = mu_IS, sd = 1)

ind <- (Y_IS > threshold)
p_IS <- mean(ind * w)
num_IS <- mean(Y_IS * ind * w)
ES_IS <- num_IS / p_IS

var_IS <- var(Y_IS[ind] * w[ind]) / n
lower_CI_IS <- ES_IS - 1.96 * sqrt(var_IS)
upper_CI_IS <- ES_IS + 1.96 * sqrt(var_IS)

cat("Importance Sampling (normal) : ES =", ES_IS, "\n")
# Report (results table & discussion): ES = 2.668681
# NOTE: the console output pasted in the report's appendix shows 2.66597
# instead. The two runs differ slightly; the table/discussion value is the
# one the report's conclusions are based on.
cat("95% CI : [", lower_CI_IS, ",", upper_CI_IS, "] \n\n")
# Report (table): 95% CI : [ 2.667839 , 2.669523 ]

###############################################################################
# Part 2: Exponential Model (Exp(10))
###############################################################################

set.seed(20)

lambda <- 10
VaR_exp <- -log(p) / lambda        # ~0.460517
true_ES_exp <- VaR_exp + 1 / lambda  # memoryless property -> ~0.560517
cat("True_exponential ES =", true_ES_exp, "\n\n")
# Report: True_exponential ES = 0.560517

# --- Crude Monte Carlo for Exp(10): complete-vector approach ---
X <- rexp(n, rate = lambda)
MCvalues_exp <- rep(0, n)
for (i in 1:n) {
  if (X[i] > VaR_exp) {
    MCvalues_exp[i] <- X[i]
  }
}

ES_MC_exp <- mean(MCvalues_exp) / p
var_MC_exp <- var(MCvalues_exp) / n
lower_CI_exp <- ES_MC_exp - 1.96 * sqrt(var_MC_exp)
upper_CI_exp <- ES_MC_exp + 1.96 * sqrt(var_MC_exp)

cat("Crude MC (exponential) : ES =", ES_MC_exp, "\n")
# Report: Crude MC (exponential) : ES = 0.4823476 (table: 0.482348)
cat("95% CI : [", lower_CI_exp, ",", upper_CI_exp, "] \n\n")
# Report: 95% CI : [ 0.4813239 , 0.4833713 ]

# --- Importance Sampling for Exp(10) ---
# Proposal: Exp(lambda_IS) with lambda_IS = 10 (same rate as the target,
# as in the report). Weights: w(x) = (lambda/lambda_IS) *
# exp(-(lambda - lambda_IS) * x).
lambda_IS <- 10
X_IS <- rexp(n, rate = lambda_IS)
w_exp <- (lambda / lambda_IS) * exp(-(lambda - lambda_IS) * X_IS)

ind_exp <- (X_IS > VaR_exp)
p_IS_exp <- mean(ind_exp * w_exp)
num_IS_exp <- mean(X_IS * ind_exp * w_exp)
ES_IS_exp <- num_IS_exp / p_IS_exp

var_IS_exp <- var(X_IS[ind_exp] * w_exp[ind_exp]) / n
lower_CI_IS_exp <- ES_IS_exp - 1.96 * sqrt(var_IS_exp)
upper_CI_IS_exp <- ES_IS_exp + 1.96 * sqrt(var_IS_exp)

cat("Importance Sampling (exponential) : ES =", ES_IS_exp, "\n")
# Report: table gives 0.5562955; the discussion text quotes 0.5545475.
# (The table's CI upper bound equals the estimate, which looks like a
# copy/paste typo in the original report.)
cat("95% CI : [", lower_CI_IS_exp, ",", upper_CI_IS_exp, "] \n")
# Report: 95% CI : [ 0.5527996 , 0.5562955 ]

###############################################################################
# Summary of the report's results
# ---------------------------------------------------------------------------
# Model        Method                        True_ES    Estimate   95% CI
# Normal       Crude MC (complete vector)    2.665214   3.075687   [3.070056, 3.081318]
# Normal       Crude MC (direct extraction)  2.665214   2.67451    [2.668107, 2.680914]
# Normal       Importance Sampling           2.665214   2.668681   [2.667839, 2.669523]
# Exponential  Crude MC                      0.560517   0.482348   [0.4813239, 0.4833713]
# Exponential  Importance Sampling           0.560517   0.5562955  [0.5527996, 0.5562955]
#
# Takeaway (report): importance sampling matches the theoretical ES closely
# when the proposal is well calibrated to the tail; the complete-vector
# crude MC overestimates ES for the normal model and underestimates it
# for the exponential model, while direct extraction of tail samples is
# more accurate.
###############################################################################
