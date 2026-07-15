# Continuous Distributions {#continuous_distributions}

When a random variable can take any value in an interval, we call it a **continuous** random variable. The height of a randomly chosen person, the time until a radioactive decay, and the noise in an electrical signal are all continuous. Because a continuous variable can take uncountably many values, the probability of any single exact value is zero. We need a new tool: the **probability density function**.

The example program for this chapter is in the file **05_continuous_distributions.lisp**.

## From Probability Mass to Probability Density

For a discrete random variable, the PMF gives P(X = x) directly. For a continuous random variable, P(X = x) = 0 for every single point x. Instead, we use the **probability density function** (PDF) f(x). Probability is now the **area under the density curve**:

    P(a <= X <= b) = integral from a to b of f(x) dx

A valid PDF satisfies two conditions: f(x) >= 0 everywhere, and the total integral of f(x) over the entire real line equals 1. The PDF itself is not a probability and can exceed 1. Only the integral over an interval is a probability.

The cumulative distribution function works the same way as in the discrete case: F(x) = P(X <= x) = integral from negative infinity to x of f(t) dt. The relationship between the PDF and CDF is that the PDF is the derivative of the CDF.

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

## The Uniform Distribution

The **uniform distribution** on the interval [a, b] has a constant density: f(x) = 1/(b-a) for a <= x <= b, and 0 elsewhere. Every value in the interval is equally likely.

    E[X] = (a + b) / 2
    Var(X) = (b - a)^2 / 12

The mean is the midpoint of the interval, which makes sense by symmetry. The variance formula involves the factor 1/12, which comes from integrating x^2 over [a, b] against the flat density.

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

## Why This Matters

The uniform, exponential, and normal distributions are the three most commonly encountered continuous distributions. The uniform distribution is the simplest: flat and unconstrained. The exponential distribution models waiting times and decay processes. The normal distribution appears everywhere, thanks to the Central Limit Theorem that we will study in a later chapter. Understanding their PDFs, CDFs, means, and variances gives you the toolkit for working with continuous probability models.
