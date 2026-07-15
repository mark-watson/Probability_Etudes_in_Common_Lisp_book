# Continuous Distributions {#continuous_distributions}

When a random variable can take any value in an interval, we call it a **continuous** random variable. The height of a randomly chosen person, the time until a radioactive decay, and the noise in an electrical signal are all continuous. Because a continuous variable can take uncountably many values, the probability of any single exact value is zero. We need a new tool: the **probability density function**.

The example program for this chapter is in the file **05_continuous_distributions.lisp**.

## Why Continuous Random Variables?

The step from discrete to continuous random variables is more than a technical convenience. Many quantities we care about are, at least in the model we use for them, genuinely continuous. Time, position, temperature, angle, mass, and voltage are all naturally described by real numbers, and no discrete PMF can capture their full range in a clean way.

There is a philosophical subtlety here. In the physical world, everything we measure is finite-resolution, so strictly speaking we could always model measurements as taking finitely many values. But the mathematics of continuous distributions is often much simpler than the mathematics of a very fine-grained discrete grid. Calculus tools such as differentiation, integration, and change of variable apply naturally to continuous distributions. So we treat continuous models as an idealization that trades a small amount of realism for a large gain in mathematical tractability.

## From Probability Mass to Probability Density

For a discrete random variable, the PMF gives P(X = x) directly. For a continuous random variable, P(X = x) = 0 for every single point x. Instead, we use the **probability density function** (PDF) f(x). Probability is now the **area under the density curve**:

    P(a <= X <= b) = integral from a to b of f(x) dx

A valid PDF satisfies two conditions: f(x) >= 0 everywhere, and the total integral of f(x) over the entire real line equals 1. The PDF itself is not a probability and can exceed 1. Only the integral over an interval is a probability.

This last point is worth emphasizing. If I tell you that the PDF of a certain distribution at x = 3 is f(3) = 2.5, that number does not mean the probability of X = 3. It means that the density of probability is 2.5 per unit length near x = 3. To get a probability, you integrate over an interval; to compare two densities you compare the ratio of densities in the same neighborhood.

The cumulative distribution function works the same way as in the discrete case: F(x) = P(X <= x) = integral from negative infinity to x of f(t) dt. The relationship between the PDF and CDF is that the PDF is the derivative of the CDF.

Formally, the CDF F is what really characterizes the distribution. The PDF is defined (when it exists) as the derivative of the CDF. Not every random variable has a PDF, but every real-valued random variable has a CDF.

### The Quantile Function

Sometimes we want the inverse question: given a probability, what value of X does it correspond to? The **quantile function** Q is the inverse of the CDF:

    Q(p) = inf { x : F(x) >= p }

Q(0.5) is the **median**, the value below which half the probability lies. Q(0.9) is the ninetieth percentile. Quantile functions are central to describing distributions in statistics; boxplots, confidence intervals, and value-at-risk calculations in finance all rely on quantiles. Note that the median and the mean are different in general. For symmetric distributions like the normal they coincide, but for skewed distributions like the exponential they can be very different.

## Expectation and Variance for Continuous Variables

The definitions of expectation and variance carry over to the continuous case, with integrals replacing sums:

    E[X] = integral of x * f(x) dx
    E[g(X)] = integral of g(x) * f(x) dx
    Var(X) = E[X^2] - (E[X])^2

All the properties from the discrete case still hold: linearity of expectation, additivity of variance for independent variables, and the concentration inequalities of Markov and Chebyshev.

## Numerical Integration

Since we cannot always compute integrals in closed form, the program uses **numerical integration** with the midpoint rule. The idea is to approximate the area under a curve by summing the areas of thin rectangles:

{lang="lisp",linenos=off}
~~~~~~~~
(defun integrate-rectangle (f a b &optional (n 100000))
  "Approximate the integral of f from a to b via the midpoint rule
   with N subintervals."
  (let ((h (/ (- b a) n)))
    (* h (loop for i from 0 below n
               sum (funcall f (+ a (* h (+ i 1/2))))))))
~~~~~~~~

With 100,000 subintervals, this method is accurate enough for our purposes. The error of the midpoint rule decreases as 1/N^2 for smooth functions, so doubling the number of intervals reduces the error by a factor of 4.

More sophisticated methods (Simpson's rule, Gaussian quadrature, adaptive quadrature) achieve much higher accuracy for the same number of function evaluations, but the midpoint rule is easy to code, easy to reason about, and good enough for pedagogy. Numerical integration is the standard way to answer questions like "what is P(a <= X <= b)?" for any distribution whose CDF is not available in closed form.

## The Uniform Distribution

The **uniform distribution** on the interval [a, b] has a constant density: f(x) = 1/(b-a) for a <= x <= b, and 0 elsewhere. Every value in the interval is equally likely.

    E[X] = (a + b) / 2
    Var(X) = (b - a)^2 / 12

The mean is the midpoint of the interval, which makes sense by symmetry. The variance formula involves the factor 1/12, which comes from integrating x^2 over [a, b] against the flat density.

The uniform distribution plays a foundational role in probability: uniform(0, 1) random variables are the raw material from which nearly every random-number generator constructs other distributions. If U is uniform on [0, 1] and F is the CDF of a distribution we want, then F^(-1)(U) has that distribution. This is called the **inverse-CDF method** or **inverse-transform sampling** and it is behind much of the machinery of Monte Carlo simulation.

The program computes the mean both by formula and by numerical integration, and they match exactly:

```
=== Uniform(0, 2) ===
  E[X] (formula) =  1.0
  E[X] (numeric) =  1.0
  Var(X) (formula) = .333
  P(0.5 <= X <= 1.5) =  0.5 (exact 0.5)
```

## The Exponential Distribution

The **exponential distribution** with rate parameter lambda models waiting times in a Poisson process. Its PDF is:

    f(x) = lambda * e^(-lambda * x)  for x >= 0

The CDF has a simple closed form: F(x) = 1 - e^(-lambda * x). The mean and variance are:

    E[X] = 1 / lambda
    Var(X) = 1 / lambda^2

A higher rate means a shorter expected wait. If the rate is 2, the expected wait is 1/2. The density decays exponentially, so most of the probability mass is concentrated near zero.

The exponential distribution shares the **memoryless property** with the geometric distribution. If you have been waiting for x minutes, the distribution of the remaining wait time is still exponential with the same rate. The past does not affect the future.

    P(X > s + t | X > s) = P(X > t)

The proof is a short computation using the closed form F(x) = 1 - e^(-lambda x). The exponential distribution is the **only** continuous distribution on [0, infinity) with this property. Its natural discrete counterpart is the geometric distribution.

The exponential distribution is closely related to the **Poisson distribution** we met briefly in the previous chapter. If events happen at random times such that the count of events in any interval of length t is Poisson with mean lambda * t, and the counts in disjoint intervals are independent, then the time between successive events is exponential with rate lambda. This unified picture is called the **Poisson process** and it describes radioactive decay, arrivals at a queue, and many other phenomena.

```
=== Exponential(lambda=2) ===
  E[X] (formula) =  0.5
  E[X] (numeric) =  0.5
  Var(X) (formula) = 0.25
  P(X <= 1) = .865 (CDF) vs .865 (numeric)
  Total probability (should be 1.0):  1.0
```

## The Normal Distribution

The **normal distribution** (also called the Gaussian distribution) is the most important distribution in all of probability theory. Its PDF is the famous bell curve:

    f(x) = (1 / (sigma * sqrt(2*pi))) * exp(-(x - mu)^2 / (2 * sigma^2))

The parameter mu is the mean (the center of the bell) and sigma^2 is the variance (how wide the bell is). The **standard normal** has mu = 0 and sigma = 1.

The normal distribution is universal: it appears as the limiting distribution of sums of independent random variables (the Central Limit Theorem in Chapter 7), as the maximum-entropy distribution given fixed mean and variance, and as the equilibrium distribution of many diffusive physical processes. Its density is smooth, symmetric, and unimodal, and it has the pleasant property that a sum of independent normals is again normal.

{lang="lisp",linenos=off}
~~~~~~~~
(defun normal-pdf (mu sigma x)
  "PDF of Normal(mu, sigma^2): the bell curve."
  (let ((z (/ (- x mu) sigma)))
    (/ (exp (- (/ (* z z) 2.0)))
       (* sigma (sqrt (* 2.0 pi))))))
~~~~~~~~

The normal CDF does not have a closed form in terms of elementary functions. It involves the **error function** (erf), which is defined as an integral. The program uses a numerical approximation from Abramowitz and Stegun that is accurate to about 7 decimal places:

{lang="lisp",linenos=off}
~~~~~~~~
(defun standard-normal-cdf (x)
  "CDF of the standard normal via a numerical approximation
   (Abramowitz and Stegun 26.2.17)."
  (let* ((sign (if (>= x 0) 1 -1))
         (ax (abs x))
         (t-val (/ 1.0 (+ 1.0 (* 0.3275911d0 ax))))
         (y (* t-val (+ 0.254829592d0
                       (* t-val (+ -0.284496736d0
                                   (* t-val (+ 1.421413741d0 ...)))))))
         (erf (- 1.0 (* y (exp (- (* ax ax)))))))
    (* 0.5 (+ 1.0 (* sign erf)))))
~~~~~~~~

### The 68-95-99.7 Rule

For any normal distribution, approximately:

- 68% of the probability mass lies within 1 standard deviation of the mean
- 95% lies within 2 standard deviations
- 99.7% lies within 3 standard deviations

The program verifies this rule by numerical integration:

```
=== Standard Normal (mu=0, sigma=1) ===
  E[X] (numeric) =    0.0 (should be 0)
  E[X^2] (numeric) =    1.0 (should be 1 = Var)
  Total probability (should be 1.0):  1.0
  68% rule P(-1 <= Z <= 1) = .683 (should be ~0.6827)
  95% rule P(-2 <= Z <= 2) = .954 (should be ~0.9545)
  99.7% rule P(-3 <= Z <= 3) = .997 (should be ~0.9973)
  CDF(0) =  0.5 (should be 0.5 by symmetry)
```

The numerical results match the theoretical values closely. The CDF at 0 is exactly 0.5 because the standard normal is symmetric about its mean of 0.

## Other Named Continuous Distributions

Beyond the uniform, exponential, and normal, several other distributions appear frequently in probability and statistics. This list is not exhaustive but gives a sense of the broader landscape.

The **Gamma distribution** generalizes the exponential: it models the total waiting time until the k-th event in a Poisson process. When k = 1 it is the exponential. When k is a positive integer it is sometimes called the **Erlang distribution**.

The **Beta distribution** on the interval [0, 1] has two shape parameters a and b. It can be flat (a = b = 1, the uniform), symmetric (a = b), or heavily skewed (very different a and b). The Beta is the conjugate prior for the Bernoulli likelihood in Bayesian inference, as we will see in Chapter 9.

The **Chi-squared distribution** with k degrees of freedom is the distribution of a sum of squares of k independent standard normals. It appears in hypothesis testing (chi-squared tests) and in confidence intervals for the variance of a normal.

The **Student's t distribution** with k degrees of freedom appears when estimating the mean of a normal population from a small sample. For large k it is nearly the standard normal; for small k it has heavier tails.

The **log-normal distribution** is the distribution of exp(X) where X is normal. It has a long right tail and is often used to model quantities that are positive and heavy-tailed, such as file sizes or income distributions.

These distributions form an interconnected web: many can be derived from one another by transformations, sums, or limits. Learning them one by one is less useful than understanding the general framework of PDFs, CDFs, expectations, and variance; the specific formulas become recognizable applications of a few underlying ideas.

## Transformations of Random Variables

If X is continuous with PDF f(x) and Y = g(X) for a strictly increasing function g, then the CDF of Y is F_Y(y) = F_X(g^{-1}(y)), and the PDF of Y is:

    f_Y(y) = f_X(g^{-1}(y)) * |d g^{-1}(y) / dy|

The absolute-value derivative is called the **Jacobian** of the transformation. This change-of-variables formula is essential for turning problems about one random variable into problems about another. For instance, the CDF of a squared standard normal can be derived from this formula, leading directly to the chi-squared distribution with 1 degree of freedom.

## Why This Matters

The uniform, exponential, and normal distributions are the three most commonly encountered continuous distributions. The uniform distribution is the simplest: flat and unconstrained. The exponential distribution models waiting times and decay processes. The normal distribution appears everywhere, thanks to the Central Limit Theorem that we will study in a later chapter. Understanding their PDFs, CDFs, means, and variances gives you the toolkit for working with continuous probability models.

## Problem Set

**Problem 5.1.** Verify by numerical integration that the PDF of Uniform(0, 5) integrates to 1 over its support. Then compute P(1 <= X <= 3) both directly and using the uniform CDF, and confirm the two answers match.

**Problem 5.2.** For the exponential distribution with rate lambda = 1/2, use the CDF to compute P(X <= 2), P(X <= 5), and P(2 <= X <= 5). Confirm your answers by numerical integration of the PDF.

**Problem 5.3 (Memoryless property).** For an exponential distribution with rate lambda = 1, verify by numerical calculation that P(X > 3 | X > 1) equals P(X > 2). Then prove the general identity P(X > s + t | X > s) = P(X > t) using the CDF.

**Problem 5.4 (Standard normal probabilities).** Using the program's standard normal CDF function, compute the following:
- P(Z <= 1.96) (this is the critical value for a 95% two-sided confidence interval)
- P(-1.645 <= Z <= 1.645) (should be about 0.90)
- P(Z >= 3) (a "three sigma" event)
- P(|Z| >= 3) (the two-sided version)

**Problem 5.5 (Converting normal to standard normal).** Let X be Normal(mu = 100, sigma = 15) (a common model for IQ scores). Using standardization Z = (X - mu) / sigma, compute P(X > 130), P(85 <= X <= 115), and the value x such that P(X > x) = 0.01.

**Problem 5.6 (Median vs. mean).** For the exponential distribution with rate lambda, compute the median (the value m such that F(m) = 1/2). Compare it to the mean 1/lambda. Which is larger, and why does the exponential have this asymmetry?

**Problem 5.7 (Numerical integration accuracy).** Rerun the numerical integration of the standard normal PDF over [-6, 6] using 100, 1000, 10000, and 100000 subintervals. Record the total integral in each case and note how quickly it converges to 1. Empirically, how does the error scale with the number of subintervals?

**Problem 5.8 (Inverse-CDF sampling).** Add a function `sample-exponential` to the example program that takes a rate lambda and returns a sample from the exponential distribution using the inverse-CDF method: draw U from Uniform(0, 1) and return -ln(1 - U) / lambda. Sample 10,000 values, compute their empirical mean and variance, and compare against 1/lambda and 1/lambda^2.

**Problem 5.9 (A transformation).** Let Z be a standard normal. Let Y = exp(Z). What are the possible values of Y? Using the change-of-variables formula, derive the PDF of Y. This is the standard log-normal distribution. Draw a rough sketch of its shape.
