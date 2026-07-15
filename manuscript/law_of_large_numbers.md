# Law of Large Numbers {#lln}

If you flip a fair coin many times, the fraction of heads gets closer and closer to 1/2. If you roll a fair die many times, the average roll gets closer and closer to 3.5. This is not just an empirical observation. It is a mathematical theorem, and it is one of the pillars of probability theory.

The example program for this chapter is in the file **06_law_of_large_numbers.lisp**.

## The Sample Mean

Let X1, X2, X3, ... be independent and identically distributed (i.i.d.) random variables with finite mean mu = E[Xi]. The **sample mean** after n observations is:

    M_n = (X1 + X2 + ... + Xn) / n

By linearity of expectation, E[M_n] = mu for every n. The sample mean is always centered on the true mean. But what about its variability?

## The Weak Law of Large Numbers

The **Weak Law of Large Numbers** states that for any positive number epsilon, no matter how small:

    P(|M_n - mu| > epsilon) approaches 0 as n approaches infinity

In plain language: as the sample size grows, the probability that the sample mean deviates from the true mean by more than any fixed amount shrinks to zero. The sample mean **converges in probability** to the true mean.

## The Strong Law of Large Numbers

The **Strong Law of Large Numbers** goes further. It states that:

    M_n approaches mu almost surely (with probability 1)

This means that for almost every possible sequence of outcomes, the sample mean eventually settles down to mu and stays there. The strong law is a stronger statement than the weak law because almost-sure convergence implies convergence in probability.

## Why This Matters

The Law of Large Numbers is the theoretical justification for using sample averages as estimates of true means. When a pollster surveys 1,000 people and reports that 52% support a candidate, they are relying on the LLN: with a large enough sample, the sample proportion will be close to the true population proportion. When a scientist repeats an experiment many times and averages the results, the LLN guarantees that the average converges to the expected value.

The LLN also explains why casinos always make money in the long run. Each individual bet is random, but over thousands of bets, the average outcome converges to the house edge. The randomness averages out.

## The Simulation

The program demonstrates the LLN by simulating a fair die and a Bernoulli process, tracking the sample mean as n grows from 10 to 1,000,000:

{lang="lisp",linenos=off}
~~~~~~~~
(defun roll-die ()
  "Simulate one roll of a fair six-sided die: uniform on {1,...,6}."
  (+ 1 (random 6 *rng-state*)))

(defun sample-mean-of-rolls (n)
  "Compute M_n = average of n die rolls. By the LLN this approaches 3.5."
  (/ (loop for i below n sum (roll-die)) n))
~~~~~~~~

For the fair die, the true mean is mu = (1+2+3+4+5+6)/6 = 3.5. For the Bernoulli process with p = 0.3, the true mean is mu = 0.3.

## Running the Example

```
=== Law of Large Numbers: Fair Die (mu = 3.5) ===
       n     M_n     |M_n - mu|
  10        3.8000    0.3000
  100       3.3300    0.1700
  1000      3.4510    0.0490
  10000     3.5415    0.0415
  100000    3.5015    0.0015
  1000000   3.4960    0.0040

=== Law of Large Numbers: Bernoulli(p=0.3) ===
       n     M_n     |M_n - p|
  10        0.2000    0.1000
  100       0.3500    0.0500
  1000      0.3160    0.0160
  10000     0.3003    0.0003
  100000    0.2977    0.0023
  1000000   0.3000    0.0000
```

Watch the column |M_n - mu| shrink as n grows. With only 10 rolls, the sample mean can be off by 0.3 or more. By 100,000 rolls, the deviation is typically under 0.002. By 1,000,000, it is essentially zero. This is the Law of Large Numbers in action.

The convergence is not perfectly monotonic. You can see that the die simulation at n=10000 has a slightly larger deviation than at n=100000. This is expected: the LLN guarantees convergence in the long run, but individual samples can fluctuate. The trend is what matters.

## A Practical Observation

Notice that the deviation shrinks roughly as 1/sqrt(n), not as 1/n. Going from 10 to 100 samples reduces the error by about a factor of 3 (roughly sqrt(10)), not a factor of 10. This slow convergence rate is a fundamental limitation of averaging. To halve the error, you need about 4 times as much data. We will see this 1/sqrt(n) rate appear again in the Monte Carlo chapter.
