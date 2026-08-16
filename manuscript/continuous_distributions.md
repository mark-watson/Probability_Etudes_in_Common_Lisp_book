# Continuous Distributions {#continuous_distributions}

When a random variable can take any value in an interval, we call it a **continuous** random variable. The height of a randomly chosen person, the time until a radioactive decay, and the noise in an electrical signal are all continuous. Because a continuous variable can take uncountably many values, the probability of any single exact value is zero. We need a new tool: the **probability density function**.

The example program for this chapter is in the file **05_continuous_distributions.lisp**.

## Why Continuous Random Variables?

The step from discrete to continuous random variables is more than a technical convenience. Many quantities we care about are, at least in the model we use for them, genuinely continuous. Time, position, temperature, angle, mass, and voltage are all naturally described by real numbers, and no discrete PMF can capture their full range in a clean way.

There is a philosophical subtlety here. In the physical world, everything we measure is finite-resolution, so strictly speaking we could always model measurements as taking finitely many values. But the mathematics of continuous distributions is often much simpler than the mathematics of a very fine-grained discrete grid. Calculus tools such as differentiation, integration, and change of variable apply naturally to continuous distributions. So we treat continuous models as an idealization that trades a small amount of realism for a large gain in mathematical tractability.

## From Probability Mass to Probability Density

For a discrete random variable, the PMF gives `P(X = x)`$ directly. For a continuous random variable, `P(X = x) = 0`$ for every single point `x`$. Instead, we use the **probability density function** (PDF) `f(x)`$. Probability is now the **area under the density curve**:

```$
P(a \leq X \leq b) = \int_{a}^{b} f(x) \, dx.
```

A valid PDF satisfies two conditions: `f(x) \geq 0`$ everywhere, and `\int_{-\infty}^{\infty} f(x) \, dx = 1`$. The PDF itself is not a probability and can exceed `1`$. Only the integral over an interval is a probability.

This last point is worth emphasizing. If I tell you that the PDF of a certain distribution at `x = 3`$ is `f(3) = 2.5`$, that number does not mean the probability of `X = 3`$. It means that the density of probability is `2.5`$ per unit length near `x = 3`$. To get a probability, you integrate over an interval; to compare two densities you compare the ratio of densities in the same neighborhood.

The cumulative distribution function works the same way as in the discrete case: `F(x) = P(X \leq x) = \int_{-\infty}^{x} f(t) \, dt`$. The relationship between the PDF and CDF is that the PDF is the derivative of the CDF.

Formally, the CDF `F`$ is what really characterizes the distribution. The PDF is defined (when it exists) as the derivative of the CDF. Not every random variable has a PDF, but every real-valued random variable has a CDF.

A random variable that has a density is called **absolutely continuous**, and the density is the derivative of the CDF wherever that derivative exists. Not every continuous CDF comes from a density: there are exotic **singular** distributions whose CDF rises continuously yet has zero derivative almost everywhere. There are also **mixed** distributions that place point masses at some values and spread density over others, such as a waiting time that equals exactly `0`$ with positive probability and is otherwise continuous. The CDF handles all of these uniformly, which is the deeper reason it, not the PDF, is taken as the defining object of a distribution. Every continuous example in this book is absolutely continuous, so a density always exists.

### The Quantile Function

Sometimes we want the inverse question: given a probability, what value of `X`$ does it correspond to? The **quantile function** `Q`$ is the inverse of the CDF:

```$
Q(p) = \inf \{ x : F(x) \geq p \}.
```

`Q(0.5)`$ is the **median**, the value below which half the probability lies. `Q(0.9)`$ is the ninetieth percentile. Quantile functions are central to describing distributions in statistics; boxplots, confidence intervals, and value-at-risk calculations in finance all rely on quantiles. Note that the median and the mean are different in general. For symmetric distributions like the normal they coincide, but for skewed distributions like the exponential they can be very different.

## Expectation and Variance for Continuous Variables

The definitions of expectation and variance carry over to the continuous case, with integrals replacing sums:

```$
\begin{aligned}
E[X]         &= \int_{-\infty}^{\infty} x \, f(x) \, dx, \\
E[g(X)]      &= \int_{-\infty}^{\infty} g(x) \, f(x) \, dx, \\
\mathrm{Var}(X) &= E[X^2] - (E[X])^2.
\end{aligned}
```

All the properties from the discrete case still hold: linearity of expectation, additivity of variance for independent variables, and the concentration inequalities of Markov and Chebyshev.

## Numerical Integration

Since we cannot always compute integrals in closed form, the program uses **numerical integration** with the midpoint rule. The idea is to approximate the area under a curve by summing the areas of thin rectangles:

```lisp
(defun integrate-rectangle (f a b &optional (n 100000))
  "Approximate the integral of f from a to b via the midpoint rule
   with N subintervals."
  (let ((h (/ (- b a) n)))
    (* h (loop for i from 0 below n
               sum (funcall f (+ a (* h (+ i 1/2))))))))
```

With `100{,}000`$ subintervals, this method is accurate enough for our purposes. The error of the midpoint rule decreases as `1/N^2`$ for smooth functions, so doubling the number of intervals reduces the error by a factor of `4`$.

More sophisticated methods (Simpson's rule, Gaussian quadrature, adaptive quadrature) achieve much higher accuracy for the same number of function evaluations, but the midpoint rule is easy to code, easy to reason about, and good enough for pedagogy. Numerical integration is the standard way to answer questions like "what is P(a <= X <= b)?" for any distribution whose CDF is not available in closed form.

To make the accuracy difference concrete, the program also includes **Simpson's rule**, which fits a parabola to each pair of subintervals and has error `O(1/N^4)`$ instead of `O(1/N^2)`$:

```lisp
(defun integrate-simpson (f a b &optional (n 100000))
  "Composite Simpson's rule: error O(1/N^4)."
  (let* ((n (if (evenp n) n (1+ n)))
         (h (/ (- b a) n))
         (s (+ (funcall f a) (funcall f b))))
    (loop for i from 1 below n
          do (incf s (* (if (oddp i) 4.0d0 2.0d0) (funcall f (+ a (* i h))))))
    (* (/ h 3.0d0) s)))
```

Integrating `e^x`$ over `[0, 1]`$ (exact value `e - 1`$) exposes the two convergence rates. Each time `N`$ doubles, the midpoint error falls by about `4`$ and the Simpson error by about `16`$:

```
=== Numerical Integration: Midpoint vs Simpson ===
  Integrating e^x over [0, 1] (exact = e - 1 = 1.718282):
    n     midpoint error   Simpson error
    4           4.467d-3        3.701d-5
    8           1.118d-3        2.326d-6
    16          2.796d-4        1.456d-7
    32          6.992d-5        9.103d-9
  Midpoint error falls ~4x per doubling (O(1/N^2)); Simpson ~16x (O(1/N^4)).
```

(The standard normal PDF is a poor test case for this comparison. Because that density and all its derivatives are essentially zero at the ends of a wide integration interval, the plain midpoint rule already converges extremely fast, and Simpson shows no clear advantage. The comparison needs an integrand with nonzero endpoint behaviour, such as `e^x`$.)

## The Uniform Distribution

The **uniform distribution** on the interval `[a, b]`$ has a constant density: `f(x) = \frac{1}{b - a}`$ for `a \leq x \leq b`$, and `0`$ elsewhere. Every value in the interval is equally likely.

```$
E[X] = \frac{a + b}{2}, \qquad \mathrm{Var}(X) = \frac{(b - a)^2}{12}.
```

The mean is the midpoint of the interval, which makes sense by symmetry. The variance formula involves the factor `1/12`$, which comes from integrating `x^2`$ over `[a, b]`$ against the flat density.

The uniform distribution plays a foundational role in probability: `\text{Uniform}(0, 1)`$ random variables are the raw material from which nearly every random-number generator constructs other distributions. If `U`$ is uniform on `[0, 1]`$ and `F`$ is the CDF of a distribution we want, then `F^{-1}(U)`$ has that distribution. This is called the **inverse-CDF method** or **inverse-transform sampling** and it is behind much of the machinery of Monte Carlo simulation.

The program computes the mean both by formula and by numerical integration, and they match exactly:

```
=== Uniform(0, 2) ===
  E[X] (formula) =  1.0
  E[X] (numeric) =  1.0
  Var(X) (formula) = .333
  P(0.5 <= X <= 1.5) =  0.5 (exact 0.5)
```

## The Exponential Distribution

The **exponential distribution** with rate parameter `\lambda`$ models waiting times in a Poisson process. Its PDF is:

```$
f(x) = \lambda e^{-\lambda x} \quad \text{for } x \geq 0.
```

The CDF has a simple closed form: `F(x) = 1 - e^{-\lambda x}`$. The mean and variance are:

```$
E[X] = \frac{1}{\lambda}, \qquad \mathrm{Var}(X) = \frac{1}{\lambda^2}.
```

A higher rate means a shorter expected wait. If the rate is `2`$, the expected wait is `1/2`$. The density decays exponentially, so most of the probability mass is concentrated near zero.

The exponential distribution shares the **memoryless property** with the geometric distribution. If you have been waiting for `x`$ minutes, the distribution of the remaining wait time is still exponential with the same rate. The past does not affect the future.

```$
P(X > s + t \mid X > s) = P(X > t).
```

The proof is a short computation. Since `P(X > x) = e^{-\lambda x}`$,

```$
P(X > s + t \mid X > s) = \frac{P(X > s + t)}{P(X > s)} = \frac{e^{-\lambda(s + t)}}{e^{-\lambda s}} = e^{-\lambda t} = P(X > t).
```

The exponential distribution is the **only** continuous distribution on `[0, \infty)`$ with this property, and its natural discrete counterpart is the geometric distribution. The link is exact in a limit: chop time into slices of width `\Delta t`$ and let each slice be an independent Bernoulli trial with success probability `\lambda\,\Delta t`$. The number of slices until the first success is geometric, and as `\Delta t \to 0`$ the waiting time converges in distribution to `\text{Exponential}(\lambda)`$. The memorylessness of the geometric passes to the exponential in the limit.

The exponential distribution is closely related to the **Poisson distribution** we met briefly in the previous chapter. If events happen at random times such that the count of events in any interval of length `t`$ is Poisson with mean `\lambda t`$, and the counts in disjoint intervals are independent, then the time between successive events is exponential with rate `\lambda`$. This unified picture is called the **Poisson process** and it describes radioactive decay, arrivals at a queue, and many other phenomena.

```
=== Exponential(lambda=2) ===
  E[X] (formula) =  0.5
  E[X] (numeric) =  0.5
  Var(X) (formula) = 0.25
  P(X <= 1) = .865 (CDF) vs .865 (numeric)
  Total probability (should be 1.0): .9999
```

The total probability comes out as `.9999`$ rather than a clean `1.0`$ because we integrate only over `[0, 20]`$ and the midpoint rule slightly underestimates the integral of a convex decreasing density. The missing mass is `e^{-40} \approx 4 \times 10^{-18}`$, so the shortfall we see is numerical, not a gap in the tail.

### Sampling by Inverse Transform

The uniform distribution is also the raw material for generating samples from other distributions. The **inverse-transform method** takes a `U \sim \text{Uniform}(0, 1)`$ draw and returns `F^{-1}(U)`$, which then has CDF `F`$. For the exponential, `F(x) = 1 - e^{-\lambda x}`$ inverts to `x = -\ln(1 - U)/\lambda`$:

```lisp
(defun sample-exponential (lambda-rate)
  "Draw Exponential(lambda) by inverse transform: -ln(1-U)/lambda."
  (/ (- (log (- 1.0d0 (random 1.0d0 *rng-state*)))) lambda-rate))
```

Drawing `100{,}000`$ samples this way and computing their mean and variance recovers `1/\lambda`$ and `1/\lambda^2`$ (Problem 5.8):

```
=== Inverse-Transform Sampling: Exponential(lambda=2) ===
  Drew 100000 samples via -ln(1-U)/lambda.
  empirical mean = .50123  (1/lambda   =    0.5)
  empirical var  = .25369  (1/lambda^2 =   0.25)
```

The exact figures shift from run to run because the samples are random, but they cluster around `0.5`$ and `0.25`$.

## The Normal Distribution

The **normal distribution** (also called the Gaussian distribution) is the most important distribution in all of probability theory. Its PDF is the famous bell curve:

```$
f(x) = \frac{1}{\sigma \sqrt{2\pi}} \exp\!\left( -\frac{(x - \mu)^2}{2 \sigma^2} \right).
```

The parameter `\mu`$ is the mean (the center of the bell) and `\sigma^2`$ is the variance (how wide the bell is). The **standard normal** has `\mu = 0`$ and `\sigma = 1`$.

The normal distribution is universal: it appears as the limiting distribution of sums of independent random variables (the Central Limit Theorem in Chapter 7), as the maximum-entropy distribution given fixed mean and variance, and as the equilibrium distribution of many diffusive physical processes. Its density is smooth, symmetric, and unimodal, and it has the pleasant property that a sum of independent normals is again normal.

```lisp
(defun normal-pdf (mu sigma x)
  "PDF of Normal(mu, sigma^2): the bell curve."
  (let ((z (/ (- x mu) sigma)))
    (/ (exp (- (/ (* z z) 2.0)))
       (* sigma (sqrt (* 2.0 pi))))))
```

The normal CDF does not have a closed form in terms of elementary functions. It is written through the **error function** (erf) by the identity `\Phi(x) = \tfrac{1}{2}\big(1 + \mathrm{erf}(x/\sqrt{2})\big)`$. The program approximates erf with a rational-times-Gaussian formula from Abramowitz and Stegun (7.1.26), whose absolute error is below `1.5 \times 10^{-7}`$:

```lisp
(defun standard-normal-cdf (x)
  "CDF of the standard normal, Phi(x) = 0.5 (1 + erf(x / sqrt 2)), using the
   Abramowitz & Stegun 7.1.26 approximation to erf."
  (let* ((sign (if (>= x 0) 1 -1))
         (z (/ (abs x) (sqrt 2.0d0)))       ; erf argument: |x| / sqrt 2
         (t-val (/ 1.0 (+ 1.0 (* 0.3275911d0 z))))
         (y (* t-val (+ 0.254829592d0
                       (* t-val (+ -0.284496736d0
                                   (* t-val (+ 1.421413741d0 ...)))))))
         (erf (- 1.0 (* y (exp (- (* z z)))))))
    (* 0.5 (+ 1.0 (* sign erf)))))
```

The `x/\sqrt{2}`$ scaling is what turns the error function into the standard normal CDF; without it the code would return the CDF of a normal with variance `1/2`$. The modern successor to Abramowitz and Stegun is the NIST Digital Library of Mathematical Functions at [dlmf.nist.gov](https://dlmf.nist.gov).

### The Moment Generating Function of the Normal

The normal distribution has moment generating function

```$
M_X(t) = \exp\!\left( \mu t + \tfrac{1}{2}\sigma^2 t^2 \right).
```

Differentiating at `t = 0`$ returns the mean `\mu`$, and the second derivative, after subtracting `\mu^2`$, returns the variance `\sigma^2`$. The exponential-of-a-quadratic shape explains two facts stated above. First, if `X \sim \text{Normal}(\mu_1, \sigma_1^2)`$ and `Y \sim \text{Normal}(\mu_2, \sigma_2^2)`$ are independent, then multiplying their MGFs adds the exponents,

```$
M_{X+Y}(t) = \exp\!\left( (\mu_1 + \mu_2)\,t + \tfrac{1}{2}(\sigma_1^2 + \sigma_2^2)\,t^2 \right),
```

which is again a normal MGF, so `X + Y \sim \text{Normal}(\mu_1 + \mu_2,\ \sigma_1^2 + \sigma_2^2)`$. The family is closed under adding independent members: means add and variances add. Second, an affine map `aX + b`$ is normal with mean `a\mu + b`$ and variance `a^2 \sigma^2`$; the special case `Z = (X - \mu)/\sigma`$ is the standardization that turns any normal into the standard normal. This closure under sums and affine maps is exactly what makes the normal the natural limit in the Central Limit Theorem of the next chapter.

### The 68-95-99.7 Rule

For any normal distribution, approximately:

- 68% of the probability mass lies within `1`$ standard deviation of the mean
- 95% lies within `2`$ standard deviations
- 99.7% lies within `3`$ standard deviations

The program verifies this rule two independent ways: by integrating the PDF over each interval, and by differencing the `standard-normal-cdf` at the interval endpoints. The two columns agree:

```
=== Standard Normal (mu=0, sigma=1) ===
  E[X] (numeric) =    0.0 (should be 0)
  E[X^2] (numeric) =    1.0 (should be 1 = Var)
  Total probability (should be 1.0):  1.0
  68% rule   P(-1<=Z<=1) = .683 (integ) .683 (CDF)  (~0.6827)
  95% rule   P(-2<=Z<=2) = .954 (integ) .954 (CDF)  (~0.9545)
  99.7% rule P(-3<=Z<=3) = .997 (integ) .997 (CDF)  (~0.9973)
  CDF(0) =  0.5 (should be 0.5 by symmetry)
```

The two methods matching is a useful check on both: the integration and the closed-form CDF approximation are computed by completely different code, yet they land on the same 68-95-99.7 figures. The CDF at `0`$ is exactly `0.5`$ because the standard normal is symmetric about its mean of `0`$.

## Other Named Continuous Distributions

Beyond the uniform, exponential, and normal, several other distributions appear frequently in probability and statistics. This list is not exhaustive but gives a sense of the broader landscape.

The **Gamma distribution** generalizes the exponential: it models the total waiting time until the `k`$-th event in a Poisson process. When `k = 1`$ it is the exponential. When `k`$ is a positive integer it is sometimes called the **Erlang distribution**.

The **Beta distribution** on the interval `[0, 1]`$ has two shape parameters `a`$ and `b`$. It can be flat (`a = b = 1`$, the uniform), symmetric (`a = b`$), or heavily skewed (very different `a`$ and `b`$). The Beta is the conjugate prior for the Bernoulli likelihood in Bayesian inference, as we will see in Chapter 9.

The **Chi-squared distribution** with `k`$ degrees of freedom is the distribution of a sum of squares of `k`$ independent standard normals. It appears in hypothesis testing (chi-squared tests) and in confidence intervals for the variance of a normal.

The **Student's t distribution** with `k`$ degrees of freedom appears when estimating the mean of a normal population from a small sample. For large `k`$ it is nearly the standard normal; for small `k`$ it has heavier tails.

The **log-normal distribution** is the distribution of `e^X`$ where `X`$ is normal. It has a long right tail and is often used to model quantities that are positive and heavy-tailed, such as file sizes or income distributions.

These distributions form an interconnected web: many can be derived from one another by transformations, sums, or limits. Learning them one by one is less useful than understanding the general framework of PDFs, CDFs, expectations, and variance; the specific formulas become recognizable applications of a few underlying ideas.

## The Maximum Entropy Viewpoint

A single principle picks out all three of our main distributions at once. The **differential entropy** of a continuous distribution with density `f`$ is

```$
h(X) = -\int_{-\infty}^{\infty} f(x) \ln f(x)\, dx,
```

a measure of how spread out, or how uncommitted, the distribution is. Among all distributions consistent with a given set of constraints, the one that maximizes `h`$ is the least presumptuous choice: it adds no structure beyond what the constraints force. The three workhorse distributions are exactly these maximum-entropy answers:

- Constrained only to live on a bounded interval `[a, b]`$, the maximum-entropy distribution is `\text{Uniform}(a, b)`$.
- Constrained to `[0, \infty)`$ with a fixed mean, it is the `\text{Exponential}`$ distribution.
- Constrained to the whole real line with a fixed mean and variance, it is the `\text{Normal}`$ distribution.

This is why the three appear so often. Each is the most honest distribution to assume when all you know is a support, a mean, or a mean and a variance. The maximum-entropy principle recurs throughout statistical physics, information theory, and Bayesian modeling as a systematic way to turn partial knowledge into a full distribution.

## Transformations of Random Variables

If `X`$ is continuous with PDF `f_X(x)`$ and `Y = g(X)`$ for a strictly increasing function `g`$, then the CDF of `Y`$ is `F_Y(y) = F_X(g^{-1}(y))`$, and the PDF of `Y`$ is:

```$
f_Y(y) = f_X(g^{-1}(y)) \left| \frac{d\, g^{-1}(y)}{dy} \right|.
```

The absolute-value derivative is called the **Jacobian** of the transformation. This change-of-variables formula is essential for turning problems about one random variable into problems about another. For instance, the CDF of a squared standard normal can be derived from this formula, leading directly to the chi-squared distribution with `1`$ degree of freedom.

## Why This Matters

The uniform, exponential, and normal distributions are the three most commonly encountered continuous distributions. The uniform distribution is the simplest: flat and unconstrained. The exponential distribution models waiting times and decay processes. The normal distribution appears everywhere, thanks to the Central Limit Theorem that we will study in a later chapter. Understanding their PDFs, CDFs, means, and variances gives you the toolkit for working with continuous probability models.

## Problem Set

**Problem 5.1.** Verify by numerical integration that the PDF of `\text{Uniform}(0, 5)`$ integrates to `1`$ over its support. Then compute `P(1 \leq X \leq 3)`$ both directly and using the uniform CDF, and confirm the two answers match.

**Problem 5.2.** For the exponential distribution with rate `\lambda = 1/2`$, use the CDF to compute `P(X \leq 2)`$, `P(X \leq 5)`$, and `P(2 \leq X \leq 5)`$. Confirm your answers by numerical integration of the PDF.

**Problem 5.3 (Memoryless property).** For an exponential distribution with rate `\lambda = 1`$, verify by numerical calculation that `P(X > 3 \mid X > 1)`$ equals `P(X > 2)`$. Then prove the general identity `P(X > s + t \mid X > s) = P(X > t)`$ using the CDF.

**Problem 5.4 (Standard normal probabilities).** Using the program's standard normal CDF function, compute the following:
- `P(Z \leq 1.96)`$ (this is the critical value for a 95% two-sided confidence interval)
- `P(-1.645 \leq Z \leq 1.645)`$ (should be about `0.90`$)
- `P(Z \geq 3)`$ (a "three sigma" event)
- `P(|Z| \geq 3)`$ (the two-sided version)

**Problem 5.5 (Converting normal to standard normal).** Let `X`$ be `\text{Normal}(\mu = 100,\, \sigma = 15)`$ (a common model for IQ scores). Using standardization `Z = (X - \mu)/\sigma`$, compute `P(X > 130)`$, `P(85 \leq X \leq 115)`$, and the value `x`$ such that `P(X > x) = 0.01`$.

**Problem 5.6 (Median vs. mean).** For the exponential distribution with rate `\lambda`$, compute the median (the value `m`$ such that `F(m) = 1/2`$). Compare it to the mean `1/\lambda`$. Which is larger, and why does the exponential have this asymmetry?

**Problem 5.7 (Numerical integration accuracy).** Rerun the numerical integration of the standard normal PDF over `[-6, 6]`$ using `100, 1000, 10000`$, and `100000`$ subintervals. Record the total integral in each case and note how quickly it converges to `1`$. Empirically, how does the error scale with the number of subintervals?

**Problem 5.8 (Inverse-CDF sampling).** Add a function `sample-exponential` to the example program that takes a rate `\lambda`$ and returns a sample from the exponential distribution using the inverse-CDF method: draw `U`$ from `\text{Uniform}(0, 1)`$ and return `-\ln(1 - U)/\lambda`$. Sample `10{,}000`$ values, compute their empirical mean and variance, and compare against `1/\lambda`$ and `1/\lambda^2`$.

**Problem 5.9 (A transformation).** Let `Z`$ be a standard normal. Let `Y = e^Z`$. What are the possible values of `Y`$? Using the change-of-variables formula, derive the PDF of `Y`$. This is the standard log-normal distribution. Draw a rough sketch of its shape.
