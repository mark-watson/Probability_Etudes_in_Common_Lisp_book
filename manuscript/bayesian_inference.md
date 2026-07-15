# Bayesian Inference {#bayesian_inference}

In the chapter on conditional probability we encountered Bayes' theorem as a way to reverse the direction of conditioning. **Bayesian inference** extends this idea to a full framework for learning from data. We start with a belief about an unknown parameter, then update that belief as evidence arrives.

The example program for this chapter is in the file **09_bayesian_inference.lisp**.

## The Bayesian Framework

In Bayesian inference, an unknown parameter theta is treated as a random variable with its own probability distribution. The core equation is:

    p(theta | data) is proportional to p(data | theta) * p(theta)

Reading from right to left:

- **p(theta)** is the **prior**: what we believe about theta before seeing any data.
- **p(data | theta)** is the **likelihood**: how likely the observed data is if theta has a particular value.
- **p(theta | data)** is the **posterior**: what we believe about theta after seeing the data.

The proportionality constant is p(data), which normalizes the posterior so it sums (or integrates) to 1. In practice, we often work with the unnormalized posterior and normalize at the end.

## Conjugate Priors

For certain combinations of likelihood and prior, the posterior is in the same family as the prior. This is called **conjugacy**, and it makes Bayesian updating remarkably simple.

For a Bernoulli or binomial likelihood with unknown success probability theta, the conjugate prior is the **Beta distribution**. If the prior is Beta(a, b) and we observe s successes and f failures, the posterior is:

    Beta(a + s, b + f)

The Beta(a, b) distribution has a density on the interval [0, 1]:

    f(theta) = theta^(a-1) * (1-theta)^(b-1) / B(a, b)

where B(a, b) is the Beta function (the normalizing constant). The parameters a and b act as **pseudo-counts**: a - 1 prior successes and b - 1 prior failures. The mean and variance are:

    E[theta] = a / (a + b)
    Var(theta) = a * b / ((a+b)^2 * (a+b+1))

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

## The Update Rule

The beauty of conjugacy is that the update rule is just addition. Each observed success increments a, and each observed failure increments b:

{lang="lisp",linenos=off}
~~~~~~~~
(defun bayesian-update (prior-a prior-b successes failures)
  "Conjugate update for a Beta prior with Bernoulli/binomial data.
   Prior  Beta(a, b)  ->  Posterior Beta(a + s, b + f)."
  (values (+ prior-a successes) (+ prior-b failures)))
~~~~~~~~

This is the same result whether we update all at once (batch) or one observation at a time (sequential). Conjugacy guarantees that the final posterior is the same either way.

## The Example: Estimating a Coin's Bias

We start with a Beta(1, 1) prior, which is the uniform distribution on [0, 1]. This represents maximum ignorance: before seeing any data, every value of theta is equally likely.

Then we simulate 100 flips of a coin with true bias theta = 0.7 and update our posterior. The program shows the posterior at several intermediate stages:

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

- **Prior**: Beta(1, 1) with mean 0.5 and variance 0.083. We know nothing yet.
- **After 10 flips**: Beta(8, 4) with mean 0.667 and variance 0.017. The mean has shifted toward 0.7, and the variance has shrunk by a factor of 5.
- **After 50 flips**: Beta(37, 15) with mean 0.712 and variance 0.004. We are getting closer, and the variance is still shrinking.
- **After 100 flips**: Beta(75, 27) with mean 0.735 and variance 0.002. The posterior mean is close to the true value of 0.7, and the variance is tiny.

The posterior mean (0.735) is not exactly 0.7 because we only have 100 data points. With 1,000 flips, it would be even closer. The variance continues to shrink toward zero as more data arrives, reflecting increasing confidence in our estimate.

## The Laplace Rule of Succession

With a uniform prior Beta(1, 1), the posterior mean after s successes and f failures is:

    E[theta | data] = (s + 1) / (s + f + 2)

This is the famous **Laplace rule of succession**. If you have seen s successes in s + f trials, your best estimate of the success probability is (s+1)/(s+f+2), not s/(s+f). The "+1" and "+2" come from the prior pseudo-counts. This rule prevents overconfidence from small samples: if you flip a coin once and get heads, the rule estimates the bias as 2/3 rather than 1.

## Why This Matters

Bayesian inference is the foundation of modern machine learning in many domains. Spam filters use it to classify emails. Medical tests use it to interpret results. A/B testing platforms use it to decide which variant is better. Recommendation systems use it to model user preferences.

The key insight of Bayesian inference is that learning is the process of updating beliefs. You start with what you know (the prior), you observe data (the likelihood), and you update your knowledge (the posterior). This cycle can be repeated indefinitely: each posterior becomes the prior for the next round of data. As data accumulates, the posterior concentrates around the true value, and the prior's influence fades away. This convergence is guaranteed by Bernoulli's theorem, a special case of the Law of Large Numbers applied to the posterior distribution.
