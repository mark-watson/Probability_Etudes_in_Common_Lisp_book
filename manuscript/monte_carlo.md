# Monte Carlo Methods {#monte_carlo}

Some quantities are easy to describe but hard to compute exactly. The **Monte Carlo** approach is to estimate them by repeated random sampling. The theoretical justification is the Law of Large Numbers: if we can express a quantity as an expectation, then the sample average of repeated random draws converges to that expectation.

The example program for this chapter is in the file **08_monte_carlo.lisp**.

## A Short History

The idea of using random sampling to estimate mathematical quantities has surprisingly deep roots. In 1777, the French naturalist Georges-Louis Leclerc, Comte de Buffon, posed a problem about dropping a needle onto a floor of parallel lines and computed the probability that the needle crosses a line. His answer involved pi, giving an early example of a randomized geometric estimator for a mathematical constant.

The modern practice of Monte Carlo methods was developed in the 1940s at Los Alamos National Laboratory by Stanislaw Ulam, John von Neumann, and Nicholas Metropolis, in the context of the Manhattan Project. They needed to solve integrals arising in neutron transport that had no analytic solution, and Ulam suggested simulating the neutron trajectories randomly on the ENIAC computer. The name "Monte Carlo" was chosen by Metropolis as a nod to Ulam's uncle, who often gambled at the Monte Carlo casino in Monaco.

Since then, Monte Carlo methods have become one of the most important computational techniques in all of science. They power simulations in physics, chemistry, biology, finance, engineering, and machine learning. The Markov chain Monte Carlo methods that underlie modern Bayesian statistics are direct descendants of the original Los Alamos work.

## The Core Idea

A Monte Carlo estimator works as follows. We want to estimate some quantity `Q`$. We find a random variable `X`$ whose expected value `E[X]`$ equals `Q`$. Then we draw many independent samples of `X`$ and average them:

```$
\hat{Q} = \frac{X_1 + X_2 + \cdots + X_n}{n}.
```

By the Law of Large Numbers, `\hat{Q}`$ converges to `Q`$ as `n`$ grows. The variance of the estimator is `\mathrm{Var}(X)/n`$, so the standard error decreases as `1/\sqrt{n}`$. By the Central Limit Theorem, `\hat{Q}`$ is approximately normal for large `n`$, so we can construct confidence intervals for `Q`$ using the normal-based methods from Chapter 7.

The design freedom in Monte Carlo lies in choosing the random variable `X`$. Any random variable whose expectation equals `Q`$ is a valid estimator, but different choices have different variances and therefore different accuracies for the same sample size. Much of the theory of Monte Carlo is really about finding good `X`$.

Two properties make this estimator trustworthy. It is **unbiased**: `E[\hat{Q}] = \frac{1}{n}\sum_i E[X_i] = Q`$ for every `n`$, so it is centered on the right answer rather than merely approaching it. And because it is unbiased, its **mean squared error** equals its variance,

```$
E\big[(\hat{Q} - Q)^2\big] = \mathrm{Var}(\hat{Q}) = \frac{\mathrm{Var}(X)}{n},
```

so the root-mean-square error is exactly the standard error `\sqrt{\mathrm{Var}(X)/n}`$ that the program reports. Reducing error therefore means reducing `\mathrm{Var}(X)`$, which is the whole point of the variance-reduction techniques below. Not every Monte Carlo estimator is unbiased, though: estimators built as ratios or other nonlinear functions of averages carry a bias of order `1/n`$, which the `1/\sqrt{n}`$ standard error dominates for large `n`$ but which matters for small samples.

## Monte Carlo Integration

The most common use of the core idea is computing integrals, because any integral is an expectation in disguise. To evaluate

```$
I = \int_{a}^{b} h(x)\, dx,
```

write it as `I = (b - a)\, E[h(U)]`$ with `U \sim \text{Uniform}(a, b)`$, since the uniform density on `[a, b]`$ is the constant `1/(b - a)`$. The estimator draws `U_1, \ldots, U_n`$ uniformly and averages:

```$
\hat{I} = \frac{b - a}{n} \sum_{i=1}^{n} h(U_i).
```

More generally, if `f`$ is any density we can sample from, then `\int g(x) f(x)\, dx = E_f[g(X)]`$, estimated by the plain average of `g(X_i)`$ over draws `X_i \sim f`$. The choice of sampling density is where the design freedom lives, because the same integral can be rewritten against *any* density `f`$ that is positive wherever the integrand is nonzero:

```$
\int g(x)\, dx = \int \frac{g(x)}{f(x)}\, f(x)\, dx = E_f\!\left[\frac{g(X)}{f(X)}\right].
```

Sampling from an `f`$ concentrated where `g`$ is large is exactly the **importance sampling** we return to below; picking `f`$ is the same as picking the estimator's variance.

## Estimating Pi

Our example estimates the value of `\pi`$. Consider a unit square `[0, 1] \times [0, 1]`$ (area `1`$) and the quarter disk of radius `1`$ centered at the origin (area `\pi/4`$). If we throw uniformly random points into the square, the probability that a point lands inside the quarter disk equals the ratio of the areas:

```$
P(\text{point in disk}) = \frac{\pi/4}{1} = \frac{\pi}{4}.
```

Let `X_i = 1`$ if the `i`$-th point lands in the disk, `0`$ otherwise. Then the `X_i`$ are i.i.d. `\text{Bernoulli}(\pi/4)`$, and by the LLN:

```$
\frac{1}{n} \sum_{i=1}^{n} X_i \to \frac{\pi}{4}.
```

So our estimator is:

```$
\hat{\pi} = 4 \cdot (\text{fraction of points inside the disk}).
```

```lisp
(defun random-point-in-unit-square ()
  "Return (x, y) with x,y ~ Uniform(0,1) independently."
  (values (random 1.0d0 *rng-state*)
          (random 1.0d0 *rng-state*)))

(defun in-quarter-disk-p (x y)
  "Is (x,y) inside the quarter disk of radius 1?
   The condition is x^2 + y^2 <= 1."
  (<= (+ (* x x) (* y y)) 1.0d0))

(defun estimate-pi (n)
  "Monte Carlo estimate of pi using n random points.
   Estimator = 4 * (count inside disk) / n."
  (let ((inside 0))
    (dotimes (i n)
      (multiple-value-bind (x y) (random-point-in-unit-square)
        (when (in-quarter-disk-p x y)
          (incf inside))))
    (* 4.0d0 (/ inside n 1.0d0))))
```

A point is inside the quarter disk if `x^2 + y^2 \leq 1`$, which is the Pythagorean distance from the origin. Uniform sampling over the square gives each region a probability equal to its area, which is the key property that makes the ratio of counts equal the ratio of areas.

## Standard Error

How accurate is our estimate? For a `\text{Bernoulli}(\pi/4)`$ variable, the standard error of the sample mean `M_n`$ is `\sqrt{p(1-p)/n}`$. The estimator is `4 M_n`$, so its standard error is:

```$
\mathrm{SE} = 4 \sqrt{\frac{\hat{p}(1 - \hat{p})}{n}},
```

where `\hat{p}`$ is the sample proportion. The program reports both the estimate and the standard error:

```lisp
(defun estimate-pi-with-se (n)
  "Estimate pi and also report a standard error for the estimate."
  (let* ((inside 0)
         (p-hat 0.0d0)
         (estimate 4.0d0)
         (se 0.0d0))
    (dotimes (i n)
      (multiple-value-bind (x y) (random-point-in-unit-square)
        (when (in-quarter-disk-p x y)
          (incf inside))))
    (setf p-hat (/ inside n 1.0d0))
    (setf estimate (* 4.0d0 p-hat))
    (setf se (* 4.0d0 (sqrt (/ (* p-hat (- 1.0d0 p-hat)) n 1.0d0))))
    (values estimate se)))
```

By the Central Limit Theorem, an approximate 95% confidence interval for `\pi`$ is `\hat{\pi} \pm 1.96 \cdot \mathrm{SE}`$. This lets us report the estimate along with a principled measure of its uncertainty.

## Running the Example

```
=== Monte Carlo Estimation of pi ===
True pi = 3.141593

       n     estimate     error    std-error
  1000      3.032000  0.109593  0.054175
  10000     3.144000  0.002407  0.016405
  100000    3.137640  0.003953  0.005202
  1000000   3.138252  0.003341  0.001645

Note: the error scales roughly as 1/sqrt(n). Doubling n
shrinks the standard error by ~0.707x; halving the error needs 4x n.
```

With `1{,}000`$ points, the estimate is off by about `0.1`$. With `1{,}000{,}000`$ points, the error drops to about `0.003`$. The standard error column tracks the actual error well: the error is typically within `1`$ or `2`$ standard errors of the true value. The random state is reseeded on each run, so your exact estimates will differ, while the `1/\sqrt{n}`$ shrinkage of the standard error stays the same.

## The Cost of Monte Carlo

The standard error decreases as `1/\sqrt{n}`$. This is a slow convergence rate. To halve the error, you need `4`$ times as many samples. To reduce the error by a factor of `10`$, you need `100`$ times as many samples.

This is the fundamental tradeoff of Monte Carlo methods. They are incredibly general: you can estimate almost anything by random sampling. But they are slow to converge. For problems where exact computation is feasible, you should prefer exact methods. Monte Carlo shines when exact computation is intractable, such as high-dimensional integrals or complex probabilistic models.

An important observation: the `1/\sqrt{n}`$ rate does **not** depend on the dimension of the problem. Classical numerical integration (Simpson's rule, Gaussian quadrature) has convergence rates that get worse in higher dimensions. In `d`$ dimensions, deterministic quadrature typically has an error that scales like `N^{-r/d}`$ for some fixed `r`$ depending on smoothness. When `d`$ is large (say, `d = 100`$), this deterministic error decreases painfully slowly with the number of function evaluations. Monte Carlo, in contrast, still has the same `1/\sqrt{n}`$ rate regardless of dimension. This is why Monte Carlo dominates in high-dimensional problems like statistical physics, Bayesian inference in complex models, and reinforcement learning.

## Variance Reduction

A large amount of Monte Carlo research is devoted to reducing variance without increasing the number of samples. Four common techniques:

**Antithetic variates**: pair each sample with its "opposite" so that variability partly cancels. For example, when sampling `U`$ from `\text{Uniform}(0, 1)`$, also use `1 - U`$. If the estimator is monotonic in `U`$, the two are negatively correlated and their average has lower variance.

**Control variates**: subtract a related quantity whose expectation is known. If we want `E[f(X)]`$ and know `E[g(X)]`$ exactly, we can estimate `E[f(X) - c(g(X) - E[g(X)])]`$ for a well-chosen `c`$. If `f`$ and `g`$ are correlated, the modified estimator has lower variance.

**Importance sampling**: sample from a different distribution and reweight. Especially useful when the region of interest (say, a tail event) has low probability under the natural distribution.

**Stratified sampling**: partition the sample space into strata and sample proportionally from each.

These techniques can dramatically improve Monte Carlo performance in practice, often by a factor of `10`$ or more. In our `\pi`$ example, we could use antithetic variates by pairing `(U_1, U_2)`$ with `(1 - U_1,\, 1 - U_2)`$; the correlation between the two indicator variables would reduce the variance of the estimator.

## Quasi-Random Sequences

An alternative to true Monte Carlo is **quasi-Monte Carlo**, which uses low-discrepancy sequences like Sobol or Halton sequences instead of pseudo-random numbers. These sequences are more evenly spread than random points and can give convergence rates closer to `1/n`$ rather than `1/\sqrt{n}`$ for smooth integrands. Quasi-Monte Carlo is widely used in high-dimensional finance and computer graphics, but its analysis is more delicate than classical Monte Carlo, and independence-based tools like the CLT do not apply directly.

## Buffon's Needle: An Older Pi Estimator

Buffon's original problem was to estimate the probability that a needle of length `L`$ dropped onto a floor with parallel lines spaced distance `d \geq L`$ apart will cross a line. Buffon showed that this probability is `\frac{2L}{\pi d}`$. Rearranging, a Monte Carlo estimate of `\pi`$ from `n`$ dropped needles that produce `k`$ crossings is:

```$
\hat{\pi} = \frac{2 L n}{d k}.
```

This is an even older and just as valid Monte Carlo estimator for `\pi`$ as the disk method. It has different variance properties and is a nice teaching example.

The program implements both the antithetic-variates trick and Buffon's needle, then compares three estimators by repeating each one `200`$ times and reporting the empirical standard deviation of the estimates:

```
=== Variance Reduction and Buffon's Needle ===
  Comparing estimators over 200 repetitions of 5000 points each:
    plain disk       mean= 3.1390  empirical SD= 0.0231
    antithetic disk  mean= 3.1400  empirical SD= 0.0202
    Buffon's needle  mean= 3.1373  empirical SD= 0.0334
```

The antithetic estimator has a smaller standard deviation than the plain one at the same cost, so pairing each `(u, v)`$ with its reflection `(1 - u, 1 - v)`$ genuinely reduces variance (Problem 8.6). Buffon's needle, using the same number of random draws, has the largest standard deviation of the three: a valid estimator but a less efficient one (Problem 8.4). The exact numbers vary from run to run.

## Applications Beyond Pi

Monte Carlo methods are used in physics (particle transport simulations), finance (option pricing), engineering (reliability analysis), and machine learning (Bayesian inference, reinforcement learning). The basic idea is always the same: express the quantity of interest as an expectation, then estimate it by averaging random samples.

Some concrete examples:

- **Physics**: simulate the trajectories of neutrons in a nuclear reactor to estimate the fraction that escape without absorption.
- **Finance**: simulate future price paths of an underlying asset to estimate the price of a complex option, especially path-dependent options.
- **Machine learning**: use Monte Carlo to approximate Bayesian posterior distributions in models where the posterior has no closed form (Markov chain Monte Carlo, particle filters).
- **Statistics**: bootstrap to estimate the variability of a sample statistic without needing a parametric model.
- **Computer graphics**: path tracing algorithms estimate the integral of light over all paths from a scene to the camera.

## Why This Matters

The `\pi`$ estimation example is pedagogically simple, but the same principle scales to problems of enormous complexity. The only requirement is that you can simulate the random variable whose expectation you want. Once you can do that, the Law of Large Numbers does the rest, and the Central Limit Theorem tells you how uncertain your estimate is.

## Problem Set

**Problem 8.1.** Run the example program with `n = 100, 400, 1600, 6400, 25600`$ sample points. Record the standard error at each step and verify that it shrinks by roughly the factor of `2`$ you would expect (since sample size quadruples each time and SE scales as `1/\sqrt{n}`$).

**Problem 8.2 (Reproducibility).** Modify the example to accept an optional seed for the random state and use it to reproduce the same estimate multiple times. Why is reproducibility particularly important in scientific Monte Carlo work?

**Problem 8.3 (Estimating an integral).** Use Monte Carlo to estimate `\int_{0}^{1} e^{-x^2/2} \, dx`$. Compare against the exact value, which is `\sqrt{2\pi}\,(\Phi(1) - \Phi(0))`$ where `\Phi`$ is the standard normal CDF. What sample size do you need to get within `0.001`$ of the true answer with 95% confidence?

**Problem 8.4 (Buffon's needle).** Implement Buffon's needle in Common Lisp. Drop a needle of length `L = 1`$ onto a floor with lines spaced `d = 1`$ apart. Estimate `\pi`$ from `n = 100{,}000`$ drops and report the standard error. Compare the accuracy to the disk method with the same `n`$.

**Problem 8.5 (Volume of a d-dimensional ball).** Use Monte Carlo to estimate the volume of the unit ball in `d = 5`$ dimensions. The exact answer is `\pi^{5/2}/\Gamma(7/2)`$, roughly `5.264`$. Now try `d = 10`$ and `d = 20`$. As `d`$ grows, most of the volume of the enclosing hypercube lies outside the ball, so the naive Monte Carlo estimator becomes very inefficient. Comment on your observations.

**Problem 8.6 (Antithetic variates).** Modify the `\pi`$ estimator to use antithetic variates: for each random `(U_1, U_2)`$, also use `(1 - U_1,\, 1 - U_2)`$. Report the estimate and standard error. Compare against the plain estimator with the same total sample size.

**Problem 8.7 (Confidence interval).** For `n = 10{,}000`$ sample points in the `\pi`$ estimation, report a 95% confidence interval for `\pi`$ using the CLT-based formula `\hat{\pi} \pm 1.96 \cdot \mathrm{SE}`$. Does the true `\pi = 3.14159\ldots`$ lie inside your interval?

**Problem 8.8 (High-dimensional integration).** Use Monte Carlo to estimate the average value of the function `f(x_1, x_2, \ldots, x_{10}) = \prod_{i=1}^{10} \sin(\pi x_i)`$ over the unit hypercube `[0, 1]^{10}`$. The exact answer is `(2/\pi)^{10}`$. Compare the Monte Carlo estimate against a direct grid-based numerical integration attempt.

**Problem 8.9 (Coding exercise).** Add a general function `monte-carlo-mean` that takes a thunk (a function of zero arguments returning a sample) and a sample size `n`$, and returns the estimated mean along with a 95% confidence interval using the CLT. Test it on the disk indicator for `\pi`$ and on the exponential from previous chapters.
