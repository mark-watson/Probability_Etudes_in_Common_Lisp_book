# Probability Theory Experiments

A series of self-contained Common Lisp programs for studying probability theory.
Each file is independent: load it on its own and it runs a demonstration.
The code comments explain the underlying **probability theory** — not the Lisp
mechanics.

## Running the Examples

Examples 01 through 10 are self-contained batch programs: load one and it runs
its demonstration. Using LispWorks (invoked here as `lw`):

```
lw -eval "(progn (load \"01_basic_probability.lisp\") (quit))"
```

or SBCL:

```
sbcl --script 01_basic_probability.lisp
```

Replace `01_basic_probability` with any of the example names below.

Example 11 is different: it is an interactive REPL, not a batch script. Load it
and type `demo` at the `ppl>` prompt:

```
rlwrap sbcl --load 11_probabilistic_dsl.lisp
```

The source and its comments assume a UTF-8 terminal, since a few comments use
symbols such as Σ and π.

## Examples

The examples are ordered to build concepts progressively.

### 01 — Basic Probability
**File:** `01_basic_probability.lisp`

Sample spaces, events, and the Kolmogorov axioms of probability.
Demonstrates the classical definition of probability
(`P(A) = |A| / |Omega|`) using the sample space of two dice, and verifies
the complement rule `P(not A) = 1 - P(A)`.

### 02 — Conditional Probability
**File:** `02_conditional_probability.lisp`

Conditional probability `P(A|B) = P(A ∩ B) / P(B)`, independence of events,
the law of total probability, and Bayes' theorem. Includes the classic
medical-screening example showing how a positive test for a rare disease
yields a surprisingly low posterior probability.

### 03 — Discrete Random Variables
**File:** `03_discrete_random_variables.lisp`

Probability mass functions (PMF), cumulative distribution functions (CDF),
expected value `E[X] = Σ x·p(x)`, and variance `Var(X) = E[X²] - (E[X])²`.
Applies these to a loaded die, a Bernoulli distribution, and the sum of two
fair dice.

### 04 — Binomial and Geometric Distributions
**File:** `04_binomial_geometric.lisp`

The binomial distribution `Binomial(n, p)` as the count of successes in `n`
independent Bernoulli trials, and the geometric distribution as the number of
trials until the first success. Verifies the memoryless property of the
geometric distribution: `P(Y > m+n | Y > m) = P(Y > n)`.

### 05 — Continuous Distributions
**File:** `05_continuous_distributions.lisp`

Probability density functions (PDF) and the relationship `P(a ≤ X ≤ b) = ∫ f(x) dx`.
Covers the Uniform, Exponential, and Normal distributions. Computes means and
variances both by closed-form formulas and by numerical integration, and
verifies the 68-95-99.7 empirical rule for the normal distribution.

### 06 — Law of Large Numbers
**File:** `06_law_of_large_numbers.lisp`

Empirical demonstration of the (Weak) Law of Large Numbers: the sample mean
`M_n = (X₁ + ... + Xₙ)/n` converges to the true mean `μ` as `n` grows.
Tracks `M_n` for a fair die (μ = 3.5) and a Bernoulli process across
increasing sample sizes from 10 to 1,000,000.

### 07 — Central Limit Theorem
**File:** `07_central_limit_theorem.lisp`

The CLT states that the sample mean of i.i.d. variables (with finite variance)
is approximately normal for large `n`, regardless of the original distribution.
Demonstrates by drawing sample means from a Bernoulli(0.5) source — a very
non-normal, two-point distribution — and showing the histogram of means is
bell-shaped, with the empirical mean and variance matching the CLT predictions.

### 08 — Monte Carlo Estimation
**File:** `08_monte_carlo.lisp`

Monte Carlo methods estimate quantities as expectations and rely on the LLN
for convergence. Estimates π by throwing random points into a unit square and
measuring the fraction that land in the quarter disk (probability = π/4).
Reports estimates and standard errors across sample sizes, illustrating the
slow `1/√n` convergence rate.

### 09 — Bayesian Inference
**File:** `09_bayesian_inference.lisp`

Bayesian belief updating using a conjugate Beta prior with a Bernoulli
likelihood. Starting from a uniform `Beta(1,1)` prior, updates the posterior
after observing simulated coin flips (true bias 0.7). Shows how the posterior
mean converges to the true parameter and the variance shrinks as data
accumulates (Bernoulli's theorem).

### 10 — Markov Chains
**File:** `10_markov_chains.lisp`

Discrete-time Markov chains, the Markov property, transition matrices, and
stationary distributions. Uses a two-state weather model (Sunny/Rainy) and
shows the distribution converging to the unique stationary distribution from
any starting state — both by long-run iteration and by solving `π = πP`
exactly.

### 11 — A Probabilistic Programming Language
**File:** `11_probabilistic_dsl.lisp`

A small but complete probabilistic programming language. It provides a macro
DSL (`defmodel`, `sample`, `observe`) for declaring Bayesian models and three
inference engines written from scratch: random-walk Metropolis-Hastings,
Hamiltonian Monte Carlo with forward-mode automatic differentiation, and
mean-field variational inference. A terminal REPL, text plots, and the ESS and
R-hat diagnostics let you fit and inspect the built-in coin, mean-variance, and
linear-regression models. Unlike the other examples, this one is interactive.

## Concepts Covered

| Example | Key Concepts |
|---------|-------------|
| 01 | Sample space, events, Kolmogorov axioms, complement rule |
| 02 | Conditional probability, independence, total probability, Bayes' theorem |
| 03 | PMF, CDF, expectation, variance, linearity of expectation |
| 04 | Bernoulli trials, binomial distribution, geometric distribution, memorylessness |
| 05 | PDF, CDF, numerical integration, uniform/exponential/normal distributions |
| 06 | Law of Large Numbers, sample mean convergence |
| 07 | Central Limit Theorem, standardization, normal approximation |
| 08 | Monte Carlo estimation, standard error, convergence rate |
| 09 | Bayesian inference, conjugate priors, posterior updating |
| 10 | Markov chains, transition matrices, stationary distribution, ergodicity |
| 11 | Probabilistic programming, MCMC, Hamiltonian Monte Carlo, variational inference, automatic differentiation |
