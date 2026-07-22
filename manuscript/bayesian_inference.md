# Bayesian Inference {#bayesian_inference}

In the chapter on conditional probability we encountered Bayes' theorem as a way to reverse the direction of conditioning. **Bayesian inference** extends this idea to a full framework for learning from data. We start with a belief about an unknown parameter, then update that belief as evidence arrives.

The example program for this chapter is in the file **09_bayesian_inference.lisp**.

## Two Paradigms of Statistics

Statistical inference broadly divides into two paradigms: **frequentist** and **Bayesian**. Understanding the difference is important because they answer subtly different questions.

The **frequentist** view treats a parameter as an unknown but fixed constant. The data are random, and any statement about the parameter must be phrased in terms of the sampling distribution of an estimator. A frequentist 95% confidence interval means: if the experiment were repeated many times, 95% of the intervals we construct this way would contain the true parameter. It says nothing directly about the probability that the specific interval we computed contains the true value.

The **Bayesian** view treats the parameter itself as a random variable with a probability distribution that expresses our belief. Data are observed, but the parameter distribution is what we update. A Bayesian 95% credible interval means: given the data, there is a 95% probability that the parameter lies in the interval, under our prior beliefs.

Neither paradigm is universally "correct"; they are two different frameworks for reasoning under uncertainty. Bayesian methods are natural when we have meaningful prior information or when we want to make direct probability statements about parameters. Frequentist methods can be more objective when there is disagreement about the prior. In modern practice, most working statisticians and machine-learning engineers use both, choosing the paradigm that best fits the problem at hand.

## The Bayesian Framework

In Bayesian inference, an unknown parameter `\theta`$ is treated as a random variable with its own probability distribution. The core equation is:

```$
p(\theta \mid \text{data}) \propto p(\text{data} \mid \theta) \, p(\theta).
```

Reading from right to left:

- `p(\theta)`$ is the **prior**: what we believe about `\theta`$ before seeing any data.
- `p(\text{data} \mid \theta)`$ is the **likelihood**: how likely the observed data is if `\theta`$ has a particular value.
- `p(\theta \mid \text{data})`$ is the **posterior**: what we believe about `\theta`$ after seeing the data.

The proportionality constant is `p(\text{data})`$, which normalizes the posterior so it sums (or integrates) to `1`$. In practice, we often work with the unnormalized posterior and normalize at the end.

The posterior distribution is the complete answer that a Bayesian inference provides. From it we can compute point estimates (like the posterior mean or median), credible intervals, and predictions about future data. The posterior is a much richer object than a single point estimate; it captures how much we know and how much we still do not know.

## The Choice of Prior

The prior expresses beliefs about theta before observing data. In some problems, we have genuine prior information (from previous studies, physical constraints, or expert opinion), and the prior encodes that. In other problems, we want the data to speak for itself and use an **uninformative** prior.

Several classical choices of uninformative prior appear in the literature:

- **Uniform prior**: assign equal density to all values of `\theta`$. This is what we use in the example (`\text{Beta}(1, 1)`$) for a coin bias. Uniform priors seem to convey no information, but they are not invariant under reparameterization: a uniform prior on the probability `p`$ is not uniform on the odds `p/(1 - p)`$.
- **Jeffreys prior**: constructed to be invariant under reparameterization. For a Bernoulli likelihood, the Jeffreys prior is `\text{Beta}(1/2,\, 1/2)`$, which is peaked near `0`$ and `1`$.
- **Reference priors**: derived to maximize the information gain from the data, again invariant under transformations.

The choice of uninformative prior can matter for small sample sizes but usually washes out as more data arrive. This convergence of the posterior toward the true value regardless of the prior is one of the most attractive features of Bayesian inference.

The prior can also express strong beliefs. In a clinical trial with a new drug, a prior might reflect the fact that most new drugs do not work well; this "skeptical" prior would then require substantial evidence in the data before concluding that the drug is effective. This is not cheating; it is a principled way to incorporate what we already know.

## Conjugate Priors

For certain combinations of likelihood and prior, the posterior is in the same family as the prior. This is called **conjugacy**, and it makes Bayesian updating remarkably simple.

For a Bernoulli or binomial likelihood with unknown success probability `\theta`$, the conjugate prior is the **Beta distribution**. If the prior is `\text{Beta}(a, b)`$ and we observe `s`$ successes and `f`$ failures, the posterior is:

```$
\text{Beta}(a + s,\, b + f).
```

The `\text{Beta}(a, b)`$ distribution has a density on the interval `[0, 1]`$:

```$
f(\theta) = \frac{\theta^{a-1} (1 - \theta)^{b-1}}{B(a, b)},
```

where `B(a, b)`$ is the Beta function (the normalizing constant). The parameters `a`$ and `b`$ act as **pseudo-counts**: `a - 1`$ prior successes and `b - 1`$ prior failures.

The conjugacy is a one-line calculation. The likelihood of `s`$ successes and `f`$ failures is `p(\text{data} \mid \theta) = \theta^{s}(1 - \theta)^{f}`$, and the prior density is proportional to `\theta^{a-1}(1 - \theta)^{b-1}`$. Multiplying them,

```$
p(\theta \mid \text{data}) \propto \theta^{s}(1 - \theta)^{f} \cdot \theta^{a-1}(1 - \theta)^{b-1} = \theta^{(a+s)-1}(1 - \theta)^{(b+f)-1},
```

which is the unnormalized `\text{Beta}(a + s,\, b + f)`$ density. The two factors have the same functional form in `\theta`$, so their product stays in the Beta family; only the exponents move. This is exactly why the update is just addition of counts, and why we never need to touch the awkward normalizing constant `B(a, b)`$ while updating.

The mean and variance are:

```$
E[\theta] = \frac{a}{a + b}, \qquad \mathrm{Var}(\theta) = \frac{a b}{(a + b)^2 (a + b + 1)}.
```

The program implements these formulas:

{lang="lisp",linenos=off}
~~~~~~~~
(defun beta-mean (a b)
  "Posterior/prior mean of Beta(a,b): E[theta] = a / (a + b)."
  (/ a (+ a b)))

(defun beta-variance (a b)
  "Var(theta) for Beta(a,b): a b / ((a+b)^2 (a+b+1))."
  (/ (* a b) (* (expt (+ a b) 2) (+ a b 1))))
~~~~~~~~

### Other Conjugate Pairs

Beta and Bernoulli are not the only conjugate pair. Several others show up frequently:

- **Normal-Normal**: if the data are normal with unknown mean but known variance, and the prior on the mean is normal, then the posterior on the mean is also normal.
- **Gamma-Poisson**: if the data are Poisson counts with unknown rate, and the prior on the rate is Gamma, then the posterior is also Gamma.
- **Dirichlet-Multinomial**: the multivariate generalization of Beta-Bernoulli, used for categorical distributions with more than two outcomes.

These pairs are not a coincidence. A conjugate prior exists whenever the likelihood belongs to an **exponential family**, the class of distributions whose density can be written as `p(x \mid \theta) = h(x)\exp\!\big(\eta(\theta)\cdot T(x) - A(\theta)\big)`$ for a natural parameter `\eta`$, a sufficient statistic `T`$, and a log-partition function `A`$. The Bernoulli, Poisson, normal, and multinomial are all exponential families, which is why each has a tidy conjugate partner. The conjugate prior is built to share the algebraic form of the likelihood in `\theta`$, so that multiplying prior by likelihood updates the parameters and leaves the shape untouched. The sufficient statistic `T(x)`$ is what the pseudo-counts accumulate: for the Bernoulli it is the success count, which is why our update simply adds `s`$ and `f`$.

Conjugacy is a mathematical convenience: it lets us do inference in closed form. In modern practice, however, we rarely have neat conjugate models. Bayesian inference in complex models usually uses Monte Carlo methods (Markov chain Monte Carlo, variational inference) to approximate the posterior. But conjugate models remain valuable for building intuition and as building blocks in larger hierarchical models.

## The Update Rule

The beauty of conjugacy is that the update rule is just addition. Each observed success increments a, and each observed failure increments b:

{lang="lisp",linenos=off}
~~~~~~~~
(defun bayesian-update (prior-a prior-b successes failures)
  "Conjugate update for a Beta prior with Bernoulli/binomial data.
   Prior  Beta(a, b)  ->  Posterior Beta(a + s, b + f)."
  (values (+ prior-a successes) (+ prior-b failures)))
~~~~~~~~

This is the same result whether we update all at once (batch) or one observation at a time (sequential). Conjugacy guarantees that the final posterior is the same either way. This equivalence of batch and sequential updating is a very useful property; it means we can process data as it arrives without waiting for it all to be collected.

### The Posterior Mean as a Weighted Average

The conjugate update has an interpretation that explains the prior-sensitivity behaviour we return to later. Write the prior mean as `\mu_0 = a/(a + b)`$ and the sample proportion, which is the maximum-likelihood estimate, as `\hat{\theta} = s/n`$ with `n = s + f`$. The posterior mean rearranges into a convex combination:

```$
E[\theta \mid \text{data}] = \frac{a + s}{a + b + n} = \underbrace{\frac{a + b}{a + b + n}}_{w}\,\mu_0 + \underbrace{\frac{n}{a + b + n}}_{1 - w}\,\hat{\theta}.
```

The posterior mean sits between the prior mean and the data's own estimate, with weights set by the prior strength `a + b`$ and the sample size `n`$. The Bayesian estimate is **shrunk** from the raw proportion toward the prior mean. When data are scarce (`n \ll a + b`$) the prior dominates; when data are plentiful (`n \gg a + b`$) the weight on the prior fades like `(a + b)/n`$ and the estimate reduces to the sample proportion. This one formula is the precise sense in which "the prior washes out," and it justifies reading `a + b`$ as a count of prior observations that compete on equal footing with the `n`$ real ones.

## Point Estimates and Credible Intervals

The posterior distribution is the complete answer, but sometimes we want a single number or a range as a summary.

The **posterior mean** `E[\theta \mid \text{data}]`$ is the standard point estimate. It minimizes squared-error loss.

The **maximum a posteriori (MAP) estimate** is the value of `\theta`$ that maximizes the posterior density. It corresponds to the "mode" of the posterior. When the prior is flat, the MAP estimate coincides with the maximum-likelihood estimate.

A **credible interval** is a range that contains the parameter with a specified posterior probability. A 95% credible interval `[L, U]`$ satisfies `P(L \leq \theta \leq U \mid \text{data}) = 0.95`$. Credible intervals answer the question that most practitioners really want: given what I have observed, what is a plausible range for the parameter?

## Prediction: The Posterior Predictive

Bayesian inference gives more than parameter estimates; it also gives principled predictions about future data. The **posterior predictive distribution** for a future observation `X_{\text{new}}`$ given past data `D`$ is:

```$
p(X_{\text{new}} \mid D) = \int p(X_{\text{new}} \mid \theta) \, p(\theta \mid D) \, d\theta.
```

This averages the likelihood over the posterior on `\theta`$. Unlike a frequentist point-prediction, the posterior predictive automatically accounts for uncertainty in `\theta`$.

For the Beta-Bernoulli model, if the posterior after observations is `\text{Beta}(a, b)`$, then the probability of a success on the next trial is `E[\theta \mid \text{data}] = a / (a + b)`$. This is exactly the posterior mean, which is why the posterior mean is such a natural point estimate.

## The Example: Estimating a Coin's Bias

We start with a `\text{Beta}(1, 1)`$ prior, which is the uniform distribution on `[0, 1]`$. This represents maximum ignorance: before seeing any data, every value of `\theta`$ is equally likely.

Then we simulate `100`$ flips of a coin with true bias `\theta = 0.7`$ and update our posterior. The program shows the posterior at several intermediate stages:

{lang="lisp",linenos=off}
~~~~~~~~
(defun main ()
  (let ((prior-a 1) (prior-b 1))      ; Beta(1,1) = Uniform(0,1)
    (print-beta prior-a prior-b "Prior         ")
    (let* ((rs (make-random-state t))
           (data (loop for i below 100
                       collect (if (< (random 1.0d0 rs) 0.7d0) 1 0))))
      (let ((successes (count 1 data))
            (failures (count 0 data)))
        (multiple-value-bind (pa pb) (bayesian-update prior-a prior-b
                                                       successes failures)
          (print-beta pa pb "Posterior     "))
        ;; Show intermediate stages
        ...))))
~~~~~~~~

## Running the Example

```
=== Bayesian Inference for a Coin Bias ===
Likelihood: Bernoulli(theta). Prior: Beta(1,1) = Uniform on [0,1].
True theta (used only to generate data) = 0.7.

  Prior         : Beta(1, 1)  mean= 0.5  var=0.0833

Observed 74 successes and 26 failures in 100 flips.
  Posterior     : Beta(75, 27)  mean=.735  var=0.0019

Sequential updates (belief concentrates as data arrives):
  After 10   flips: Beta(8, 4)  mean=.667  var=0.0171
  After 50   flips: Beta(37, 15)  mean=.712  var=0.0039
  After 100  flips: Beta(75, 27)  mean=.735  var=0.0019

True theta = 0.7. As n grows the posterior mean converges to the
true value and the variance shrinks to 0 (Bernoulli's theorem).
```

Watch what happens as data accumulates:

- **Prior**: `\text{Beta}(1, 1)`$ with mean `0.5`$ and variance `0.083`$. We know nothing yet.
- **After `10`$ flips**: `\text{Beta}(8, 4)`$ with mean `0.667`$ and variance `0.017`$. The mean has shifted toward `0.7`$, and the variance has shrunk by a factor of `5`$.
- **After `50`$ flips**: `\text{Beta}(37, 15)`$ with mean `0.712`$ and variance `0.004`$. We are getting closer, and the variance is still shrinking.
- **After `100`$ flips**: `\text{Beta}(75, 27)`$ with mean `0.735`$ and variance `0.002`$. The posterior mean is close to the true value of `0.7`$, and the variance is tiny.

The posterior mean (`0.735`$) is not exactly `0.7`$ because we only have `100`$ data points. With `1000`$ flips, it would be even closer. The variance continues to shrink toward zero as more data arrives, reflecting increasing confidence in our estimate.

## The Laplace Rule of Succession

With a uniform prior `\text{Beta}(1, 1)`$, the posterior mean after `s`$ successes and `f`$ failures is:

```$
E[\theta \mid \text{data}] = \frac{s + 1}{s + f + 2}.
```

This is the famous **Laplace rule of succession**. If you have seen `s`$ successes in `s + f`$ trials, your best estimate of the success probability is `(s + 1)/(s + f + 2)`$, not `s/(s + f)`$. The `+1`$ and `+2`$ come from the prior pseudo-counts. This rule prevents overconfidence from small samples: if you flip a coin once and get heads, the rule estimates the bias as `2/3`$ rather than `1`$.

Laplace himself used this rule to estimate the probability that the sun would rise tomorrow given that it has risen every day so far. The result is silly if taken too literally, but the underlying mathematical idea, that small samples should not be treated as certain, is enduring and important.

## Prior Sensitivity

An important practical question is how much the choice of prior affects the posterior. The general answer is: for small data, the prior matters a lot; for large data, it washes out.

To see this concretely, compare the posterior mean after `s`$ successes in `n = s + f`$ trials for two different priors:

- `\text{Beta}(1, 1)`$ prior: posterior mean `= (s + 1)/(n + 2)`$
- `\text{Beta}(10, 10)`$ prior: posterior mean `= (s + 10)/(n + 20)`$

For `n = 10`$ and `s = 7`$, the first gives `8/12 = 0.667`$, and the second gives `17/30 = 0.567`$. These are quite different. But for `n = 10000`$ and `s = 7000`$, the first gives `7001/10002 = 0.700`$, and the second gives `7010/10020 = 0.700`$. They are essentially the same. This washout effect is a hallmark of Bayesian inference with lots of data.

## Why This Matters

Bayesian inference is the foundation of modern machine learning in many domains. Spam filters use it to classify emails. Medical tests use it to interpret results. A/B testing platforms use it to decide which variant is better. Recommendation systems use it to model user preferences.

The key insight of Bayesian inference is that learning is the process of updating beliefs. You start with what you know (the prior), you observe data (the likelihood), and you update your knowledge (the posterior). This cycle can be repeated indefinitely: each posterior becomes the prior for the next round of data. As data accumulates, the posterior concentrates around the true value, and the prior's influence fades away. This convergence is guaranteed by Bernoulli's theorem, a special case of the Law of Large Numbers applied to the posterior distribution.

## Problem Set

**Problem 9.1.** Starting from a `\text{Beta}(1, 1)`$ prior and observing `3`$ successes in `10`$ flips, compute the posterior parameters, the posterior mean, and the posterior variance. Compare with the raw sample proportion `3/10`$.

**Problem 9.2 (Sequential vs. batch).** Observe the sequence of coin flips `1, 0, 1, 1, 0`$. Starting from `\text{Beta}(1, 1)`$, update the posterior after each flip and record the `(a, b)`$ pair at each step. Then compute the batch update after all `5`$ flips. Verify that the two answers agree.

**Problem 9.3 (Prior sensitivity).** Suppose the data are `7`$ successes in `10`$ trials. Compute the posterior mean for three priors: `\text{Beta}(1, 1), \text{Beta}(10, 10)`$, and `\text{Beta}(1, 100)`$. Interpret each of these priors in words and explain how they influence the posterior estimate.

**Problem 9.4 (Credible interval).** For a `\text{Beta}(75, 27)`$ posterior, compute an approximate 95% credible interval `[L, U]`$. Since the Beta CDF is not elementary, use numerical integration of the density (as the example program does for normalization) and find the values `L`$ and `U`$ such that `P(\theta < L) = 0.025`$ and `P(\theta > U) = 0.025`$.

**Problem 9.5 (MAP estimate).** For a `\text{Beta}(a, b)`$ distribution with `a, b > 1`$, the mode is at `(a - 1)/(a + b - 2)`$. Verify this by differentiating the log-density and setting the derivative to zero. Compute the MAP estimate for a `\text{Beta}(75, 27)`$ posterior and compare with the posterior mean.

**Problem 9.6 (Posterior predictive).** For a `\text{Beta}(75, 27)`$ posterior, what is the probability that the next flip is heads? Now compute the probability that the next two flips are both heads. Hint: the answer is not `(75/102)^2`$. Use the fact that the two flips are not independent given the posterior on `\theta`$; average the joint likelihood over the posterior.

**Problem 9.7 (Jeffreys prior).** Rerun the example with a `\text{Beta}(1/2,\, 1/2)`$ prior instead of `\text{Beta}(1, 1)`$. Compare the posterior after `10`$ flips (`7`$ successes, `3`$ failures) under both priors. When do the two priors give noticeably different answers?

**Problem 9.8 (A biased coin).** Suppose you have strong prior belief that a coin is fair, expressed as a `\text{Beta}(50, 50)`$ prior (mean `0.5`$ and quite concentrated). You then observe `30`$ heads in `40`$ flips. What is your posterior mean? Is the data enough to override your prior?

**Problem 9.9 (Coding exercise).** Extend the example program to compute the posterior predictive probability of the next observation being a success. Also extend it to plot (or print a text histogram of) the posterior density at several intermediate stages so you can visually watch the belief concentrate.

**Problem 9.10 (A different conjugate pair).** Suppose you are counting the number of shooting stars per hour and model it as `\text{Poisson}(\lambda)`$. A conjugate prior for the rate `\lambda`$ is `\text{Gamma}(\alpha, \beta)`$, and the posterior after observing counts `x_1, \ldots, x_n`$ is `\text{Gamma}(\alpha + \sum x_i,\, \beta + n)`$. Starting from a `\text{Gamma}(1, 1)`$ prior and observing counts `5, 7, 4, 6`$, compute the posterior parameters and the posterior mean. Compare with the sample mean of the counts.
