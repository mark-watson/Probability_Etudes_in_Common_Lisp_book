# Binomial and Geometric Distributions {#binomial_geometric}

Two of the most important discrete distributions arise from repeating the same simple experiment over and over. If we flip a biased coin repeatedly, we can ask two natural questions: how many successes happen in a fixed number of flips, and how long do we wait for the first success? These questions lead to the **binomial** and **geometric** distributions.

The example program for this chapter is in the file **04_binomial_geometric.lisp**.

## Bernoulli Trials

A **Bernoulli trial** is the simplest random experiment with two outcomes: "success" (with probability `p`$) and "failure" (with probability `1 - p`$). If we let `X = 1`$ on success and `X = 0`$ on failure, then `X`$ has the `\text{Bernoulli}(p)`$ distribution with:

```$
E[X] = p, \qquad \mathrm{Var}(X) = p(1 - p).
```

The Bernoulli trial is the atom from which we build more complex distributions. The Bernoulli distribution itself is named after Jakob Bernoulli, the seventeenth-century Swiss mathematician whose posthumous *Ars Conjectandi* (1713) contains the first published proof of what we would now call the Weak Law of Large Numbers for Bernoulli trials.

A key modeling question is when a Bernoulli trial is a reasonable idealization of reality. It works well for a coin flip or a fair random selection, but many real situations only approximate the assumptions of a Bernoulli trial. Two trials are only truly "the same" if the mechanism that produces them is unchanged and the outcome of one does not influence the other. Whenever we build a binomial or geometric model, we are implicitly making these assumptions.

## The Binomial Distribution

If we run `n`$ **independent** `\text{Bernoulli}(p)`$ trials and count the number of successes, the count `X`$ has the **binomial distribution**, written `\text{Binomial}(n, p)`$. The probability of getting exactly `k`$ successes is:

```$
P(X = k) = \binom{n}{k} p^k (1 - p)^{n - k},
```

where `\binom{n}{k} = \frac{n!}{k!\,(n-k)!}`$ is the **binomial coefficient**, the number of ways to choose which `k`$ of the `n`$ trials are the successes. The factor `p^k (1-p)^{n-k}`$ is the probability of any one specific sequence with `k`$ successes and `n - k`$ failures, and `\binom{n}{k}`$ counts how many such sequences there are.

Two ingredients are worth pausing over. First, the factor `p^k (1-p)^{n-k}`$ arises because the trials are independent: the joint probability of a specific sequence of outcomes is the product of the individual probabilities. If the trials were correlated, this factorization would fail. Second, the binomial coefficient `\binom{n}{k}`$ is exactly the number of arrangements of the `k`$ successes among the `n`$ positions. Multiplying gives us the total probability of the event "exactly `k`$ successes."

The program computes binomial coefficients iteratively to avoid huge intermediate factorials:

{lang="lisp",linenos=off}
~~~~~~~~
(defun binomial-coefficient (n k)
  "C(n,k) = n!/(k!(n-k)!). Computed iteratively to avoid huge
   intermediate factorials."
  (if (or (< k 0) (> k n))
      0
      (loop with result = 1
            for i from 1 to k
            do (setf result (* result (/ (+ n (- k) i) i)))
            finally (return (round result)))))
~~~~~~~~

The mean and variance of the binomial distribution are:

```$
E[X] = np, \qquad \mathrm{Var}(X) = np(1 - p).
```

The mean formula has an intuitive explanation via linearity of expectation. Each trial contributes `p`$ to the expected count on average, and there are `n`$ independent trials, so the expected total is `np`$. Similarly, since the trials are independent, the variance of the count is the sum of the variances of the individual Bernoulli indicators, which is `np(1 - p)`$.

### The Binomial as a Sum of Bernoullis

Formally, if `X_1, X_2, \ldots, X_n`$ are i.i.d. `\text{Bernoulli}(p)`$ random variables and `S = X_1 + X_2 + \cdots + X_n`$, then `S`$ has the `\text{Binomial}(n, p)`$ distribution. This decomposition is often the fastest way to prove properties of the binomial: linearity of expectation gives `E[S] = np`$ in one line, and independence gives `\mathrm{Var}(S) = np(1-p)`$ in another. Any theorem you know about sums of independent random variables applies immediately to binomial random variables.

The same decomposition hands us the moment generating function for free. A single `\text{Bernoulli}(p)`$ variable has MGF `E[e^{tX_i}] = (1 - p) + p e^{t}`$, and because the trials are independent the MGF of their sum is the product of the individual MGFs:

```$
M_S(t) = \left(1 - p + p e^{t}\right)^{n}.
```

Reading the mean and variance off the first two derivatives at `t = 0`$ recovers `np`$ and `np(1 - p)`$ yet again. The product form also proves a **reproductive property**: if `X \sim \text{Binomial}(m, p)`$ and `Y \sim \text{Binomial}(n, p)`$ are independent with the *same* `p`$, then `X + Y \sim \text{Binomial}(m + n, p)`$, since their MGFs multiply to `(1 - p + p e^{t})^{m + n}`$. The trial picture makes this obvious too: pooling `m`$ trials with `n`$ more of the same kind simply gives `m + n`$ trials.

### The Normal Approximation

For large `n`$, the binomial distribution is closely approximated by a normal distribution with mean `np`$ and variance `np(1 - p)`$. This follows from the Central Limit Theorem (Chapter 7). A useful rule of thumb is that the normal approximation is reasonable when both `np`$ and `n(1 - p)`$ are at least `5`$. For a single coin flip (`n = 1`$), the binomial is nowhere near normal; for `1000`$ flips, the normal is nearly indistinguishable from the exact binomial.

### The Poisson Limit

There is a beautiful limit that connects the binomial distribution to another famous distribution. Suppose `n`$ grows large and `p`$ shrinks such that the product `\lambda = np`$ is held constant. Then the binomial distribution converges to the **Poisson distribution** with parameter `\lambda`$:

```$
P(X = k) = \frac{\lambda^k e^{-\lambda}}{k!}.
```

This is why Poisson distributions describe rare events counted over long stretches: the number of decays of a radioactive sample in one second, the number of typos on a page, the number of calls arriving at a call center in a minute. Each event is one of a huge number of possibilities, each with tiny individual probability; the total count is nearly Poisson.

The program makes this concrete. It tabulates `\text{Binomial}(50, 0.06)`$, `\text{Binomial}(500, 0.006)`$, and `\text{Poisson}(3)`$ side by side (all with mean `\lambda = 3`$), and the binomial columns march toward the Poisson one as `n`$ grows. It also reports the **mode** of `\text{Binomial}(10, 0.3)`$, the most likely count, from the closed form `\lfloor (n+1)p \rfloor`$ (Problem 4.9).

## The Geometric Distribution

The **geometric distribution** answers a different question: if we run `\text{Bernoulli}(p)`$ trials until the **first success** occurs, how many trials do we need? If `Y`$ is the number of trials (including the successful one), then:

```$
P(Y = k) = (1 - p)^{k - 1} p, \qquad k = 1, 2, 3, \ldots
```

This formula reads as "`k - 1`$ failures followed by one success." The geometric distribution is a valid PMF because it sums to a geometric series: `\sum_{k=1}^{\infty} (1 - p)^{k-1} p = p \cdot \frac{1}{1 - (1-p)} = 1`$.

A note on conventions: some books define `Y`$ as the number of failures **before** the first success, so `Y`$ takes values `0, 1, 2, \ldots`$ rather than `1, 2, 3, \ldots`$. The formulas are almost identical but the mean shifts by `1`$. Always check which convention a book or library uses.

The mean and variance are:

```$
E[Y] = \frac{1}{p}, \qquad \mathrm{Var}(Y) = \frac{1 - p}{p^2}.
```

The mean `1/p`$ makes intuitive sense. If the success probability is `1/5`$, you wait about `5`$ trials on average. The smaller the success probability, the longer you expect to wait. The variance grows as `p`$ shrinks, so waiting times for rare events are both long and highly variable.

### Why the Mean Is 1/p

The formula `E[Y] = 1/p`$ has a derivation that uses no series at all, only the structure of the experiment. Condition on the first trial. With probability `p`$ it succeeds and `Y = 1`$. With probability `1 - p`$ it fails, one trial is spent, and the remaining wait is a fresh, statistically identical copy of `Y`$. Hence

```$
E[Y] = p \cdot 1 + (1 - p)\big(1 + E[Y]\big).
```

Solving for `E[Y]`$ gives `E[Y] = 1/p`$. The same conditioning trick applied to `E[Y^2]`$ produces the variance `(1 - p)/p^2`$ without summing a single geometric series. This self-consistency argument, in which a quantity is expressed in terms of itself one step later, is the discrete seed of the first-step analysis we will use for Markov chains in the final chapter.

The tail probability `P(Y > k)`$ has an especially clean form:

```$
P(Y > k) = (1 - p)^k.
```

This says: the probability that we still have not had a success after `k`$ trials is the probability that all `k`$ trials failed, which is `(1 - p)^k`$ by independence. We will use this formula in the proof of the memoryless property below.

## The Memoryless Property

The geometric distribution has a remarkable property called **memorylessness**. It says that the past does not affect the future. Specifically:

```$
P(Y > m + n \mid Y > m) = P(Y > n).
```

If you have already waited `m`$ trials without a success, the probability of waiting at least `n`$ more trials is the same as if you had just started. The distribution "forgets" how long you have been waiting.

The program demonstrates this property directly:

{lang="lisp",linenos=off}
~~~~~~~~
(defun demonstrate-memoryless-property (p m n)
  "Show P(Y > m+n | Y > m) = P(Y > n) for a geometric random variable Y."
  (let ((conditional (/ (geometric-tail p (+ m n)) (geometric-tail p m))))
    (format t "  Memoryless property check (p=~a, m=~a, n=~a):~%" p m n)
    (format t "    P(Y > m+n | Y > m) = ~a = ~4f~%" conditional (float conditional))
    (format t "    P(Y > n)          = ~a = ~4f~%" (geometric-tail p n)
            (float (geometric-tail p n)))))
~~~~~~~~

The proof is a one-liner. Since `P(Y > k) = (1 - p)^k`$, we have:

```$
P(Y > m + n \mid Y > m) = \frac{P(Y > m + n)}{P(Y > m)} = \frac{(1 - p)^{m + n}}{(1 - p)^m} = (1 - p)^n = P(Y > n).
```

The geometric is the only discrete distribution on `\{1, 2, 3, \ldots\}`$ with this property. In the continuous world, its analogue is the exponential distribution, which we will meet in the next chapter.

## Other Distributions from Bernoulli Trials

The binomial and geometric are two members of a small family of distributions built from Bernoulli trials. Two others are worth naming, even briefly.

The **negative binomial distribution** generalizes the geometric: instead of waiting for the first success, we wait for the `r`$-th success. The number of failures before the `r`$-th success has a negative binomial distribution. The geometric distribution is the special case `r = 1`$.

The **hypergeometric distribution** replaces sampling with replacement (where each trial is truly independent) with sampling without replacement from a finite population. If we draw n cards from a shuffled deck without replacement and count how many are aces, the count is hypergeometric, not binomial, because each draw changes the composition of the remaining deck. When the population is much larger than the sample, the hypergeometric distribution is well approximated by the binomial.

## Real-World Applications

The binomial distribution is the workhorse of any situation involving repeated trials with a fixed success probability. It underlies:

- **Opinion polling**: the number of respondents in a random sample who support a candidate.
- **Quality control**: the number of defective items in a batch.
- **A/B testing**: the number of users in a test group who click a button.
- **Genetics**: the number of offspring with a particular trait, under simple Mendelian assumptions.

The geometric distribution is the workhorse of any situation involving waiting for a first success:

- **Sales**: the number of sales calls until the first sale.
- **Networking**: the number of packet transmissions until a successful acknowledgment.
- **Reliability**: the number of trials until a machine failure, if the failure probability is constant per trial.
- **Games**: the number of rolls until you get a specific number on a die.

The memoryless property is a strong assumption in these applications. It is only true if the underlying success probability really is constant across trials. In many real settings the true success probability changes over time, and the geometric distribution is only a first approximation.

## Running the Example

```
=== Binomial Distribution: n=10, p=0.3 ===
  E[X] = n p = 3,  Var(X) = n p (1-p) = 21/10
  P(X=0) = 0.0282   CDF F(0) = 0.0282
  P(X=1) = 0.1211   CDF F(1) = 0.1493
  P(X=2) = 0.2335   CDF F(2) = 0.3828
  P(X=3) = 0.2668   CDF F(3) = 0.6496
  ...
  P(X=10) = 0.0000   CDF F(10) = 1.0000

=== Geometric Distribution: p=0.2 ===
  E[Y] = 1/p = 5,  Var(Y) = (1-p)/p^2 = 20
  P(Y=1) = 1/5 =  0.2   P(Y>1) = 4/5 =  0.8
  P(Y=2) = 4/25 = 0.16   P(Y>2) = 16/25 = 0.64
  ...

=== Memoryless Property ===
  Memoryless property check (p=1/5, m=3, n=2):
    P(Y > m+n | Y > m) = 16/25 = 0.64
    P(Y > n)          = 16/25 = 0.64

=== Poisson Limit of the Binomial (lambda = n p = 3) ===
  As n grows with n p = 3 fixed, Binomial(n, 3/n) -> Poisson(3).
   k   Binom(50,0.06)   Binom(500,0.006)   Poisson(3)
  0         0.04533          0.04934         0.04979
  1         0.14467          0.14891         0.14936
  2         0.22624          0.22427         0.22404
  3         0.23106          0.22472         0.22404
  4         0.17329          0.16854         0.16803
  5         0.10176          0.10092         0.10082
  6         0.04872          0.05026         0.05041

=== Binomial Mode ===
  Most likely k for Binomial(10, 0.3) = 3 (the peak of the PMF above)
```

The binomial distribution with `n = 10`$ and `p = 0.3`$ is peaked around `k = 3`$ (the mean), with the probability declining on both sides. The CDF rises from near `0`$ to exactly `1`$, confirming that the PMF sums to `1`$.

The geometric distribution with `p = 0.2`$ is monotonically decreasing: the most likely outcome is `Y = 1`$ (a success on the very first trial) with probability `0.2`$, and longer waits become progressively less likely.

The memoryless property check confirms that `P(Y > 5 \mid Y > 3)`$ equals `P(Y > 2)`$, both being `16/25 = 0.64`$. The past waiting time of `3`$ trials has no effect on the future.

## Why This Matters

The binomial distribution appears whenever we count successes in repeated independent trials: opinion polling, quality control, A/B testing, and many other applications. The geometric distribution models waiting times: how many calls until a sale, how many attempts until a success, how many packets until a collision. Together, these two distributions cover a huge range of practical discrete probability problems.

## Problem Set

**Problem 4.1.** A fair coin is flipped `20`$ times. Using the binomial PMF from the example program, compute the probability of getting exactly `10`$ heads, at least `15`$ heads, and at most `5`$ heads. Then verify that the probabilities of "exactly `k`$ heads" for `k = 0, 1, \ldots, 20`$ sum to `1`$.

**Problem 4.2.** A political poll interviews `1000`$ randomly chosen voters. Suppose the true fraction supporting a candidate is `0.5`$. Using the binomial distribution, what is the probability that the observed proportion in the poll lies within one percentage point of `0.5`$ (that is, between `490`$ and `510`$ supporters)? You will probably want to use the normal approximation for this one; check your answer with the exact binomial as well.

**Problem 4.3 (Comparing means and variances).** Compute the mean and variance of `\text{Binomial}(100, 0.1)`$ and of `\text{Binomial}(100, 0.5)`$. Which has the larger variance? Explain in words why symmetry (`p = 1/2`$) maximizes the variance among all `\text{Bernoulli}(p)`$ distributions.

**Problem 4.4 (Geometric waiting).** A die is rolled repeatedly. Let `Y`$ be the number of rolls until the first `6`$ appears. What is `E[Y]`$? What is `P(Y \leq 5)`$? What is `P(Y > 10)`$? Now compute `P(Y > 10 \mid Y > 5)`$ and confirm it equals `P(Y > 5)`$, the memoryless property.

**Problem 4.5 (Memorylessness and gambler's fallacy).** A slot machine hits a jackpot with probability `p`$ on each pull, independently of past pulls. A gambler has just experienced `20`$ losing pulls in a row and reasons that a jackpot must be "due" soon. Use the memoryless property to explain why this reasoning is wrong. What does the memoryless property say about the distribution of remaining pulls until the next jackpot?

**Problem 4.6 (Simulation vs. theory).** Extend the example program to simulate `100{,}000`$ trials of `\text{Binomial}(20, 0.4)`$. Compare the empirical PMF against the theoretical PMF. Also compute the empirical mean and variance and compare against `np`$ and `np(1 - p)`$.

**Problem 4.7 (Poisson approximation).** Compute the exact `\text{Binomial}(1000, 0.003)`$ probability that `X = 3`$. Then compute the Poisson approximation with `\lambda = np = 3`$, using `P(X = 3) = 3^3 e^{-3} / 3!`$. Compare the two values. The Poisson formula uses no factorial-of-`1000`$ and is much cheaper to compute.

**Problem 4.8 (Waiting for multiple successes).** In a game where each trial succeeds with probability `p = 0.2`$, what is the expected number of trials needed to accumulate `5`$ successes? Hint: by linearity of expectation, this is the sum of `5`$ independent geometric waiting times. Use this to compute the answer directly.

**Problem 4.9 (Coding exercise).** Write a function `binomial-mode` that, given `n`$ and `p`$, returns the value of `k`$ that maximizes the binomial PMF. There is a closed-form expression involving `\lfloor (n + 1) p \rfloor`$. Verify your function against a direct search for a few values of `n`$ and `p`$.
