# Central Limit Theorem {#clt}

The Law of Large Numbers tells us that the sample mean converges to the true mean. But it does not tell us how the sample mean is **distributed** around that mean. The **Central Limit Theorem** (CLT) answers this question, and its answer is one of the most remarkable results in all of mathematics: the distribution is approximately normal, no matter what the original distribution looks like.

The example program for this chapter is in the file **07_central_limit_theorem.lisp**.

## The Theorem

Let X1, X2, ... be i.i.d. random variables with mean mu and finite variance sigma^2. Let S_n = X1 + X2 + ... + Xn be their sum. The CLT says that for large n, the **standardized sum**:

    (S_n - n * mu) / (sigma * sqrt(n))

converges in distribution to the **standard normal** distribution Normal(0, 1).

Equivalently, the sample mean M_n = S_n / n is approximately Normal(mu, sigma^2 / n) for large n. Notice that the variance of M_n shrinks as 1/n: more data means a tighter distribution around the mean.

## Why This Is Remarkable

The CLT says that the original Xi can have **any** distribution with finite variance, and the sum will still become approximately normal. The Xi could be Bernoulli, uniform, exponential, or some weird custom distribution. It does not matter. Sums of enough i.i.d. random variables always tend toward the normal distribution.

This is why the normal distribution appears everywhere in nature. Many natural quantities are sums or averages of many small independent effects. Heights, measurement errors, blood pressure, and countless other quantities are approximately normal because they are built from the sum of many independent contributions.

## The Simulation

The program demonstrates the CLT with a deliberately non-normal source: a Bernoulli(0.5) distribution, which takes only the values 0 and 1. This is about as far from a bell curve as you can get. We draw sample means of size n=50 and collect 10,000 of them:

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

- Mean = mu = p = 0.5
- Variance = sigma^2 / n = p(1-p) / n = 0.25 / 50 = 0.005

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

Look at that histogram. It is bell-shaped, even though the source distribution has only two possible values. The CLT has transformed a highly non-normal distribution into a nearly normal one, just by averaging 50 independent draws.

The empirical mean (0.5011) and variance (0.0050) are very close to the CLT predictions of 0.5 and 0.005. The small deviation in the mean is normal sampling variability.

## The Theoretical Basis for Statistical Practice

The CLT is the foundation for much of classical statistics. Confidence intervals, hypothesis tests (z-tests, t-tests), and regression analysis all rely on the normal approximation that the CLT provides.

When a pollster reports a "margin of error" of plus or minus 3%, that number comes from the CLT. The sample proportion is approximately normal (by the CLT), and the margin of error is about 2 standard deviations of that normal distribution.

When a quality control engineer tests whether a manufacturing process is within spec, the CLT justifies using normal-based control charts even when the underlying process distribution is not normal.

## How Large Does n Need to Be?

A common rule of thumb is that n >= 30 is sufficient for the normal approximation to be reasonable. But this depends on the shape of the original distribution. For symmetric distributions, n=10 or even n=5 may be enough. For highly skewed distributions, you may need n=50 or more.

In our simulation, n=50 works beautifully for the symmetric Bernoulli(0.5) case. If we used a Bernoulli(0.01) instead, which is highly skewed, we would need a much larger n to see the bell shape emerge.

## Why This Matters

The Central Limit Theorem is the reason the normal distribution deserves its special status. It is not just one distribution among many. It is the **attractor** distribution for sums of independent random variables. Wherever independent effects accumulate, the normal distribution appears. This makes it the most important distribution in probability and statistics, and the CLT tells us why.
