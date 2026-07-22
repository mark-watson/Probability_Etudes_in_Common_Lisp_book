# Central Limit Theorem {#clt}

The Law of Large Numbers tells us that the sample mean converges to the true mean. But it does not tell us how the sample mean is **distributed** around that mean. The **Central Limit Theorem** (CLT) answers this question, and its answer is one of the most remarkable results in all of mathematics: the distribution is approximately normal, no matter what the original distribution looks like.

The example program for this chapter is in the file **07_central_limit_theorem.lisp**.

## A Brief History

The Central Limit Theorem has a long and distinguished history. The earliest form was proved by Abraham de Moivre in 1733 for the special case of Bernoulli trials: he showed that the binomial distribution is well approximated by a normal distribution for large n. Pierre-Simon Laplace generalized De Moivre's result in his monumental *Theorie Analytique des Probabilites* (1812) and applied it to problems in astronomy and geodesy, where he needed to reason about the distribution of averaged measurement errors.

The modern general form of the CLT, allowing arbitrary distributions with finite variance, was proved by the Russian mathematician Aleksandr Lyapunov in 1901. Later, the theorem was refined to handle non-identically-distributed variables under a technical condition known as the Lindeberg condition. Even today, generalizations of the Central Limit Theorem are an active area of research, especially for weakly dependent variables and high-dimensional problems.

## The Theorem

Let `X_1, X_2, \ldots`$ be i.i.d. random variables with mean `\mu`$ and finite variance `\sigma^2`$. Let `S_n = X_1 + X_2 + \cdots + X_n`$ be their sum. The CLT says that for large `n`$, the **standardized sum**:

```$
\frac{S_n - n \mu}{\sigma \sqrt{n}}
```

converges in distribution to the **standard normal** distribution `\text{Normal}(0, 1)`$.

Equivalently, the sample mean `M_n = S_n / n`$ is approximately `\text{Normal}(\mu,\, \sigma^2/n)`$ for large `n`$. Notice that the variance of `M_n`$ shrinks as `1/n`$: more data means a tighter distribution around the mean. But the CLT tells us more than the LLN did: it tells us the *shape* of the distribution of `M_n`$, not just that it concentrates.

Written as a limit of CDFs, the CLT says:

```$
P\!\left( \frac{S_n - n \mu}{\sigma \sqrt{n}} \leq x \right) \to \Phi(x) \quad \text{as } n \to \infty,
```

where `\Phi`$ is the CDF of the standard normal. This mode of convergence is called **convergence in distribution** or **weak convergence**. It is the weakest of the modes we surveyed in the previous chapter, and even so, it delivers enormous practical value.

## Why This Is Remarkable

The CLT says that the original `X_i`$ can have **any** distribution with finite variance, and the sum will still become approximately normal. The `X_i`$ could be Bernoulli, uniform, exponential, or some weird custom distribution. It does not matter. Sums of enough i.i.d. random variables always tend toward the normal distribution.

This is why the normal distribution appears everywhere in nature. Many natural quantities are sums or averages of many small independent effects. Heights, measurement errors, blood pressure, and countless other quantities are approximately normal because they are built from the sum of many independent contributions.

The result is worth savoring. A theorem that starts with almost no assumptions about the `X_i`$ except finite variance ends with a very specific conclusion: the standardized sum is a standard normal. It is almost as if the normal distribution is the mathematical shadow of independence and averaging.

## Why It Is True: A Sketch via Generating Functions

The moment generating function from Chapter 3 gives a clean way to see *why* the normal appears, not just that it does. Standardize the individual variables by writing `Y_i = (X_i - \mu)/\sigma`$, so each `Y_i`$ has mean `0`$ and variance `1`$. The standardized sum is

```$
Z_n = \frac{S_n - n\mu}{\sigma\sqrt{n}} = \frac{1}{\sqrt{n}}\sum_{i=1}^{n} Y_i.
```

Because the `Y_i`$ are independent, the MGF of `Z_n`$ is the `n`$-th power of one scaled MGF:

```$
M_{Z_n}(t) = \left[ M_Y\!\left( \frac{t}{\sqrt{n}} \right) \right]^{n}.
```

Now expand `M_Y`$ near `0`$. Since `M_Y(0) = 1`$, `M_Y'(0) = E[Y] = 0`$, and `M_Y''(0) = E[Y^2] = 1`$, the Taylor expansion is `M_Y(s) = 1 + \tfrac{1}{2}s^2 + o(s^2)`$. Substituting `s = t/\sqrt{n}`$ and raising to the `n`$-th power,

```$
M_{Z_n}(t) = \left( 1 + \frac{t^2}{2n} + o\!\left(\tfrac{1}{n}\right) \right)^{n} \longrightarrow e^{t^2/2}.
```

The limit `e^{t^2/2}`$ is exactly the MGF of the standard normal (set `\mu = 0`$, `\sigma = 1`$ in the normal MGF of Chapter 5). Since the MGF determines the distribution, `Z_n`$ converges to `\text{Normal}(0, 1)`$. The whole theorem comes down to one fact: only the first two moments survive the `1/\sqrt{n}`$ scaling. Every distribution with finite variance has the same quadratic leading behaviour, so they all flow to the same bell curve. The `\sqrt{n}`$ in the denominator is forced, being the unique scaling that holds the variance of `Z_n`$ at `1`$; divide by `n`$ and you get the constant `\mu`$ of the Law of Large Numbers, divide by less and the variance blows up. A fully rigorous proof replaces the MGF with the characteristic function `\varphi(t) = E[e^{itX}]`$, which exists for every distribution, but the algebra is the same.

## Rate of Convergence: Berry-Esseen

The plain CLT is a limit statement: it does not tell us how large `n`$ has to be for the normal approximation to be accurate. The **Berry-Esseen theorem** provides a quantitative bound. Under the additional assumption that the `X_i`$ have finite third absolute moment `E[|X_i|^3]`$, the maximum difference between the CDF of the standardized sum and the standard normal CDF satisfies:

```$
\sup_{x} |F_n(x) - \Phi(x)| \leq \frac{C \, \rho}{\sigma^3 \sqrt{n}},
```

for some absolute constant `C`$ (the best known value is roughly `0.4748`$) and where `\rho = E[|X_i - \mu|^3]`$. The error shrinks as `1/\sqrt{n}`$, the same rate we saw in the Law of Large Numbers. Doubling `n`$ improves the CLT approximation by a factor of `\sqrt{2}`$, roughly `1.41`$.

The key qualitative message is that the approximation works better when:

- The original distribution is nearly symmetric (small third moment).
- The variance is well-behaved (not too small relative to the third moment).
- The sample size `n`$ is large.

For heavy-tailed distributions or highly skewed distributions, the convergence can be quite slow, and much larger samples may be needed before the normal approximation is trustworthy.

## The Simulation

The program demonstrates the CLT with a deliberately non-normal source: a `\text{Bernoulli}(0.5)`$ distribution, which takes only the values `0`$ and `1`$. This is about as far from a bell curve as you can get. We draw sample means of size `n = 50`$ and collect `10{,}000`$ of them:

{lang="lisp",linenos=off}
~~~~~~~~
(defun sample-mean (p n)
  "M_n = average of n i.i.d. Bernoulli(p)."
  (/ (loop for i below n sum (bernoulli p)) n 1.0d0))

(defun collect-sample-means (p n num-samples)
  "Collect NUM-SAMPLES independent sample means, each from n Bernoulli(p)
   trials."
  (loop for s below num-samples collect (sample-mean p n)))
~~~~~~~~

The CLT predicts that these sample means should have:

- Mean `= \mu = p = 0.5`$
- Variance `= \sigma^2/n = p(1-p)/n = 0.25/50 = 0.005`$

The program computes both the theoretical predictions and the empirical values from the simulation, then prints a text histogram.

## Running the Example

```
=== Central Limit Theorem ===
Source distribution: Bernoulli(p=0.5), very non-normal.
Sample size n = 50, number of sample means = 10000

CLT predicts for the sample mean M_n:
  E[M_n]        = mu        =  0.5
  Var(M_n)      = sigma^2/n = 0.0050
Empirical from simulation:
  mean of M_n's = 0.5011
  var  of M_n's = 0.0050
  (These should be close to the CLT predictions.)

Histogram of sample means (bell-shaped = CLT at work):
   0.325 | **
   0.375 | *******
   0.425 | ****************************
   0.475 | ********************************
   0.525 | **************************************************
   0.575 | **********************
   0.625 | **************
   0.675 | **
```

Look at that histogram. It is bell-shaped, even though the source distribution has only two possible values. The CLT has transformed a highly non-normal distribution into a nearly normal one, just by averaging `50`$ independent draws.

The empirical mean (`0.5011`$) and variance (`0.0050`$) are very close to the CLT predictions of `0.5`$ and `0.005`$. The small deviation in the mean is normal sampling variability.

## The Theoretical Basis for Statistical Practice

The CLT is the foundation for much of classical statistics. Confidence intervals, hypothesis tests (z-tests, t-tests), and regression analysis all rely on the normal approximation that the CLT provides.

### Confidence Intervals

Suppose we sample `n`$ i.i.d. values with unknown mean `\mu`$ and known standard deviation `\sigma`$. By the CLT, the sample mean `M_n`$ is approximately `\text{Normal}(\mu,\, \sigma^2/n)`$. A **95% confidence interval** for `\mu`$ is:

```$
M_n \pm 1.96 \cdot \frac{\sigma}{\sqrt{n}}.
```

The number `1.96`$ is the 97.5th percentile of the standard normal (the two-sided critical value for 95% coverage). This interval is constructed so that in 95% of hypothetical repetitions of the experiment, the true mean would lie inside. The width of the interval shrinks as `1/\sqrt{n}`$, the familiar Monte Carlo scaling.

### Margin of Error in Polling

When a pollster reports a "margin of error" of plus or minus 3%, that number comes from the CLT. The sample proportion is approximately normal (by the CLT), and the margin of error is about `2`$ standard deviations of that normal distribution. For a sample of `n`$ voters, the margin of error is approximately `1/\sqrt{n}`$ in percentage points, ignoring some multiplicative constants.

This is why pollsters usually collect around `1000`$ respondents. With `n = 1000`$, the margin of error is roughly `3`$ percentage points. To halve it, they would need `n = 4000`$. To reach a margin of error of `0.3`$ percentage points, they would need `n = 100{,}000`$, which is prohibitively expensive for most polls.

### Quality Control

When a quality control engineer tests whether a manufacturing process is within spec, the CLT justifies using normal-based control charts even when the underlying process distribution is not normal. As long as the daily quantity being tracked is an average of many small independent effects, the CLT guarantees it will be approximately normal.

## How Large Does n Need to Be?

A common rule of thumb is that `n \geq 30`$ is sufficient for the normal approximation to be reasonable. But this depends on the shape of the original distribution. For symmetric distributions, `n = 10`$ or even `n = 5`$ may be enough. For highly skewed distributions, you may need `n = 50`$ or more.

For `\text{Binomial}(n, p)`$ situations, an alternative rule of thumb is that the normal approximation is good when both `np`$ and `n(1 - p)`$ exceed `5`$ (some sources use `10`$). This ensures that the distribution is far from being pushed against `0`$ or `n`$ and has enough room in both tails.

In our simulation, `n = 50`$ works beautifully for the symmetric `\text{Bernoulli}(0.5)`$ case. If we used a `\text{Bernoulli}(0.01)`$ instead, which is highly skewed, we would need a much larger `n`$ to see the bell shape emerge.

### The Continuity Correction

When the CLT approximates a *discrete* distribution, such as a binomial count, a small adjustment sharpens the result. A discrete variable puts probability on the integers, while the approximating normal spreads probability continuously, so the event `\{X \geq 60\}`$ is matched better by the normal event `\{X \geq 59.5\}`$: we widen each integer to the half-unit interval around it. To approximate `P(X \geq k)`$ for an integer-valued `X`$, apply the normal approximation to `P(X \geq k - \tfrac12)`$; to approximate `P(X \leq k)`$, use `P(X \leq k + \tfrac12)`$. This **continuity correction** typically cuts the approximation error several-fold for moderate `n`$, and it is worth applying whenever the underlying variable is a count. Problem 7.7 asks you to measure its effect against an exact binomial tail.

## Beyond i.i.d.: Extensions of the CLT

The CLT we have stated assumes i.i.d. random variables. Several extensions relax this:

- **Lindeberg-Feller CLT**: allows non-identically-distributed independent variables as long as no single term dominates in the variance.
- **Martingale CLT**: allows a specific kind of dependence (martingale differences), useful in stochastic processes.
- **Multivariate CLT**: extends to random vectors; the limit is a multivariate normal.
- **Functional CLT (Donsker's theorem)**: the entire path of a random walk converges to Brownian motion.

These extensions show how central the CLT is: even in complicated settings, sums of many small independent effects tend to look normal in some appropriate sense.

## Why This Matters

The Central Limit Theorem is the reason the normal distribution deserves its special status. It is not just one distribution among many. It is the **attractor** distribution for sums of independent random variables. Wherever independent effects accumulate, the normal distribution appears. This makes it the most important distribution in probability and statistics, and the CLT tells us why.

## Problem Set

**Problem 7.1.** For the example simulation (`n = 50, p = 0.5`$), the CLT predicts `M_n`$ is approximately `\text{Normal}(0.5,\, 0.005)`$. Compute the standard deviation of that normal (should be about `0.0707`$) and identify the range `M_n \pm 2\sigma`$. Approximately what fraction of the `10{,}000`$ simulated sample means should fall inside this range?

**Problem 7.2.** Modify the simulation to use `p = 0.05`$ (a highly skewed Bernoulli). With `n = 50`$, is the histogram of sample means still bell-shaped? Increase `n`$ to `500, 5000`$. At what point does the bell shape become clear? Explain what you observe using the Berry-Esseen intuition.

**Problem 7.3 (Uniform source).** Replace the Bernoulli source with a `\text{Uniform}(0, 1)`$ source. Draw `10{,}000`$ sample means with `n = 12`$. Confirm that the sample means look approximately `\text{Normal}(0.5,\, 1/144)`$. Because `\text{Uniform}(0, 1)`$ has mean `1/2`$ and variance `1/12`$, and the mean of `12`$ samples has variance `1/144`$, this is the CLT.

**Problem 7.4 (Exponential source).** Replace the source with an exponential of rate `1`$ (mean `1`$, variance `1`$). Draw `10{,}000`$ sample means with `n = 30, 100`$, and `300`$. Because the exponential is right-skewed, the bell shape emerges more slowly than with the uniform. Compare the empirical histograms.

**Problem 7.5 (Sample size for margin of error).** How large a random sample `n`$ do you need so that a 95% confidence interval for a Bernoulli proportion has half-width at most `0.02`$? Use the worst case `p = 0.5`$. What if the half-width must be `0.005`$?

**Problem 7.6 (CLT vs. Chebyshev).** For a `\text{Bernoulli}(0.5)`$ source with `n = 100`$, use Chebyshev's inequality to bound `P(|M_n - 0.5| \geq 0.1)`$. Then use the CLT to compute the same probability under the normal approximation. Compare the two bounds. Which is tighter?

**Problem 7.7 (Approximating a binomial tail).** Use the CLT (with continuity correction if you know it) to approximate `P(X \geq 60)`$ where `X`$ is `\text{Binomial}(100, 0.5)`$. Then compute the exact binomial probability using the example program from Chapter 4. How close are the two answers?

**Problem 7.8 (Coding exercise).** Add a function that computes the empirical CDF of the sample means. Compare it against the CLT-predicted normal CDF and print the maximum absolute difference. This is a numerical estimate of the Berry-Esseen bound in action. Verify that this maximum difference shrinks as you increase `n`$.

**Problem 7.9 (When the CLT fails).** The CLT requires finite variance. Recall the Cauchy distribution from Chapter 6. Simulate `10{,}000`$ sample means of `n = 50`$ Cauchy draws. Plot the histogram. Is it bell-shaped, or does it look like the Cauchy itself? Explain your observation.
