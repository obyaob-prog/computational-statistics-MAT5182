# Computational Statistics — MAT5182 (uOttawa, Winter 2025)

Coursework portfolio for **MAT5182 Computational Statistics** (uOttawa), covering
simulation, Monte Carlo integration, numerical optimization, bootstrap, advanced
Monte Carlo methods, and a final project on financial risk measurement.

> **How this repo was built:** the original `.Rmd` files for this course were
> lost. All R code in `code/` was **reconstructed from the graded submitted
> reports** (copied unchanged in `reports/`), which contained the code chunks,
> parameter values, seeds and printed outputs. The reconstruction is faithful to
> those reports — including their quirks (see *Limitations* below) — but the
> scripts have **not been re-run in this environment**, so please run them in
> R/RStudio and confirm the outputs match the values quoted in the comments.

## Research questions

- **Simulation:** can we generate samples from non-standard distributions
  (Beta via Gamma ratios, Pareto by inversion, Gamma by rejection, Wigner's
  semicircle law) and validate them against theory? How biased is a
  mean-lifetime estimator when lifetimes are censored at 4 years?
- **Integration:** how do crude Monte Carlo, control variates and antithetic
  sampling compare when estimating an integral — with and without a
  variance-reducing Beta(3,1) reparameterization? Same question for rare-event
  tail probabilities of a Weibull and for Asian option prices.
- **Optimization:** grid search vs Newton–Raphson vs gradient ascent on
  `log(x)/(1+x)`; analytical vs Newton–Raphson MLE for Pareto(α); EM for a
  two-component Gaussian mixture (multiple starts, local maxima) and for a
  one-parameter multinomial model.
- **Bootstrap:** normal vs percentile CIs for the median of |N(0,1)|; MLE of a
  Uniform(0,θ) endpoint under MC / bootstrap / subsampling; parametric vs
  nonparametric bootstrap for the kurtosis of Exp(1).
- **Advanced Monte Carlo:** stationary distribution of a finite Markov chain;
  Metropolis–Hastings with an exponential proposal for a non-standard density
  (effect of the proposal scale λ ∈ {0.1, 1, 10}).
- **Final project — Value-at-Risk & Expected Shortfall:** estimate 1%-tail VaR
  and ES for N(0,1) and Exp(10) losses by crude Monte Carlo vs importance
  sampling, with 95% CIs benchmarked against analytical values. Written report
  graded **20/20**; presentation **7/10** (instructor feedback: it ran over the
  time limit and tried to cover too much — the deck should stay focused).

## Repository layout

```
code/
  01-simulations.R            Assignment 1: simulation (base R)
  02-integration.R            Assignment 2: integration  (base R + ggplot2)
  03-optimization.R           Assignment 3: optimization (base R)
  04-bootstrap.R              Assignment 4: bootstrap    (base R)
  05-advanced-monte-carlo.R   Assignment 5: MCMC         (base R)
  var-expected-shortfall.R    Final project: VaR & ES    (base R)
reports/
  assignement-1-computational-stat.pdf
  assignement-2-computational-stat.docx
  assignement-3.docx
  assignement-4-bootstrap.docx
  Assignement-5.pdf
  project-Yann-OBROU.pdf                      (final project report, 20/20)
  project-Yann-OBROU-annotated-by-instructor.pdf
```

## Key findings

- **Simulation:** `U²` with U ~ Uniform(0,1) is Beta(1/2, 1) (QQ-plot on the
  y = x line); inverse-transform Pareto(2) and rejection-sampled Gamma(3/2, 1)
  samples match their theoretical densities. The censored-lifetime estimator is
  heavily biased (reported biases 2.43 / 2.90 / 1.86 for true rates
  1/4, 1/5, 1/3) — extending the 4-year follow-up reduces it.
- **Integration:** with the Beta(3,1) reparameterization all three estimators
  land within ~0.001 of the true 0.6278151 at n = 10⁵; importance sampling
  with a Γ(3, 1) proposal beats crude MC by ~15–30× relative efficiency on
  Weibull tail probabilities; a LogNormal proposal works for β < 1 where Gamma
  does not. For Asian options, control variates with S_T cut the variance of
  crude MC by ~4×.
- **Optimization:** all three optimizers agree on x* ≈ 3.5911, f(x*) ≈ 0.2785
  (gradient ascent needs a small step: 0.01 works for the log-ratio, but the
  multinomial problem required 0.0001). Pareto MLE: analytical and
  Newton–Raphson agree (mean 2.0061, SD 0.0557 over 100 replications). EM on
  the Gaussian mixture converges to *different* local maxima from different
  starts — always try multiple initializations.
- **Bootstrap:** for the median of |N(0,1)| the normal CI [0.6009, 0.6964] and
  percentile CI [0.6030, 0.7024] agree, justifying the normal approximation.
  For kurtosis of Exp(1), the parametric bootstrap (mean 8.80) reproduces the
  empirical sampling distribution (mean 8.77); the nonparametric bootstrap is
  noisier (mean 7.97) because kurtosis is sensitive to tail extremes.
- **VaR/ES project:** importance sampling with a well-calibrated proposal
  (N(2.33,1) for the normal model) matches the theoretical ES (2.6687 vs
  2.6652; 0.5563 vs 0.5605) with very tight CIs. The "complete vector" crude
  MC implementation overestimates ES for the normal model (3.0757) and
  underestimates it for the exponential (0.4823); directly extracting tail
  samples is more accurate.

## How to run

```r
# Everything is base R except Assignment 2, which needs ggplot2:
install.packages("ggplot2")

source("code/01-simulations.R")
source("code/02-integration.R")
source("code/03-optimization.R")
source("code/04-bootstrap.R")
source("code/05-advanced-monte-carlo.R")
source("code/var-expected-shortfall.R")
```

Each script is self-contained: it sets its own seeds (the same seeds recorded
in the reports) and every printed number is annotated with the value from the
corresponding report, so you can check the reproduction line by line.

## Limitations (read before citing numbers)

1. **Reconstructed, not re-run.** The code was rewritten from the reports and
   could not be executed here. Treat every figure as "expected output" until
   you run it in R and confirm.
2. **Report quirks preserved.** A few oddities from the original submissions
   were kept verbatim (with `NOTE` comments in the code): in Assignment 1 Q1
   the Gamma-ratio sample is overwritten by `rbeta()` before plotting; in
   Assignment 1 Q4 the bias column subtracts the *mean* of the lambda grid;
   in Assignment 4 Q2(b) the so-called "parametric" bootstrap actually
   resamples the original data; in Assignment 2 Q4 the "antithetic" path draws
   fresh normals and negates them instead of reusing the original path.
3. **Unresolved ambiguities.** (a) Assignment 5 Q1(c): the instructions asked
   for a chain of length 100 but the report's code uses 100,000 — both runs
   are in the script. (b) VaR project: the appendix console print for normal
   importance sampling shows ES = 2.66597 while the results table/discussion
   says 2.668681; for exponential IS the table says 0.5562955 while the text
   says 0.5545475 (and the table's CI upper bound equals the estimate, likely
   a copy/paste typo). The scripts follow the tables; re-running will settle it.
4. **Missing page.** Page 5 of the VaR project PDF is image-only and could not
   be recovered as text (likely a figure); it does not affect the code.
5. **Hand-written parts.** Proofs/derivations (A1 Q1 proof, A2 Q2, A3 Q1, A4
   theory parts, A5 Q3 posterior) were submitted on paper/photos and contain no
   R code; they are marked as such in the scripts.
