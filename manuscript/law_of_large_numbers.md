# Law of Large Numbers {#lln}

If you flip a fair coin many times, the fraction of heads gets closer and closer to 1/2. If you roll a fair die many times, the average roll gets closer and closer to 3.5. This is not just an empirical observation. It is a mathematical theorem, and it is one of the pillars of probability theory.

The example program for this chapter is in the file **06_law_of_large_numbers.lisp**.

## Historical Roots

Jakob Bernoulli, in his posthumous *Ars Conjectandi* (1713), proved the first version of what we now call the Law of Large Numbers for Bernoulli trials. He called his result the "Golden Theorem" and considered it his most important contribution. In modern language, Bernoulli's theorem says that the sample proportion of successes converges to the true probability `p`$ as the number of trials grows. Bernoulli was proud enough of his result to spend twenty years polishing it before publication.

Later refinements by Chebyshev in the nineteenth century, and eventually by Khinchin and Kolmogorov in the twentieth, extended the theorem to arbitrary i.i.d. random variables with finite mean. Chebyshev's contribution was a simple but powerful inequality (Chapter 3) that gives a proof of the Weak Law of Large Numbers in a few lines. Kolmogorov's Strong Law is more delicate but rests on the same intuition: averages of many i.i.d. quantities converge to the underlying mean.

## The Sample Mean

Let `X_1, X_2, X_3, \ldots`$ be independent and identically distributed (i.i.d.) random variables with finite mean `\mu = E[X_i]`$. The **sample mean** after `n`$ observations is:

```$
M_n = \frac{X_1 + X_2 + \cdots + X_n}{n}.
```

By linearity of expectation, `E[M_n] = \mu`$ for every `n`$. The sample mean is always centered on the true mean. But what about its variability?

If the individual variance is `\sigma^2`$, then by additivity of variance for independent variables:

```$
\mathrm{Var}(M_n) = \frac{\sigma^2}{n}.
```

The variance of the sample mean shrinks as `n`$ grows. Its standard deviation, which is the natural scale for measuring "typical" deviations, shrinks as `\sigma / \sqrt{n}`$. This `1/\sqrt{n}`$ rate is fundamental. It says the sample mean concentrates around the true mean, but only slowly: reducing the uncertainty by a factor of `10`$ requires `100`$ times more data.

## The Weak Law of Large Numbers

The **Weak Law of Large Numbers** states that for any positive number `\epsilon`$, no matter how small:

```$
P(|M_n - \mu| > \epsilon) \to 0 \quad \text{as } n \to \infty.
```

In plain language: as the sample size grows, the probability that the sample mean deviates from the true mean by more than any fixed amount shrinks to zero. The sample mean **converges in probability** to the true mean.

### A Proof Sketch Using Chebyshev

The Weak Law has a beautifully short proof using Chebyshev's inequality (Chapter 3). Applied to `M_n`$:

```$
P(|M_n - \mu| \geq \epsilon) \leq \frac{\mathrm{Var}(M_n)}{\epsilon^2} = \frac{\sigma^2}{n \, \epsilon^2}.
```

As `n`$ grows, the right-hand side goes to zero, and so the left-hand side must too. This gives convergence in probability directly. The proof relies only on Chebyshev's inequality and additivity of variance, both of which are elementary. Bernoulli's original proof used a much more intricate combinatorial argument; Chebyshev's cleaner path is the one usually taught today.

## The Strong Law of Large Numbers

The **Strong Law of Large Numbers** goes further. It states that:

```$
M_n \to \mu \quad \text{almost surely (with probability 1).}
```

This means that for almost every possible sequence of outcomes, the sample mean eventually settles down to `\mu`$ and stays there. The strong law is a stronger statement than the weak law because almost-sure convergence implies convergence in probability but not vice versa.

The difference between the two laws is subtle. The Weak Law says that at each large `n`$, the sample mean is *probably* close to `\mu`$. The Strong Law says that if we watch the whole sequence `M_1, M_2, M_3, \ldots`$, the sequence itself converges to `\mu`$ with probability `1`$. For most practical purposes the two laws give the same guarantee, but the strong law is what the theoretical statistician wants.

## Modes of Convergence

The distinction between weak and strong laws is a special case of a broader theme in probability: different notions of "convergence of random variables." A brief tour:

- **Convergence in probability**: for every `\epsilon > 0`$, `P(|X_n - X| > \epsilon) \to 0`$. This is what the Weak Law provides.
- **Almost sure convergence**: the sequence `X_n`$ converges to `X`$ pointwise, on all sample paths outside a set of probability zero. This is what the Strong Law provides.
- **Convergence in distribution**: the CDFs `F_n`$ converge to `F`$ at every continuity point of `F`$. This is the weakest notion; it is what the Central Limit Theorem in the next chapter provides.
- **Convergence in mean square**: `E[|X_n - X|^2] \to 0`$.

These modes are related by implications: almost sure implies in probability; in probability implies in distribution; mean square implies in probability. Whenever a probability textbook talks about convergence of random variables, it is worth pausing to check which mode is being used.

## When the Law of Large Numbers Fails

The Law of Large Numbers requires finite mean. When the mean does not exist, the sample mean does **not** converge to any deterministic value. The classic pathological example is the **Cauchy distribution**, which has PDF:

```$
f(x) = \frac{1}{\pi (1 + x^2)}.
```

The Cauchy distribution has infinitely heavy tails and no finite mean. Sample means of Cauchy-distributed variables do not settle down; instead, they themselves are Cauchy-distributed. Simulating a million Cauchy draws and averaging them gives an answer that is just as noisy as a single draw. This is a striking counterexample and a reminder that "average enough data" is not a universal solution.

For distributions with finite mean but infinite variance, the Weak Law still holds (in a slightly more delicate form) but the variance-based Chebyshev proof does not apply, and the rate of convergence can be much slower than `1/\sqrt{n}`$. Heavy-tailed distributions appear in finance, network traffic, and other real settings where the LLN's guarantees are still true but weaker than the standard textbook picture suggests.

## Why This Matters

The Law of Large Numbers is the theoretical justification for using sample averages as estimates of true means. When a pollster surveys 1,000 people and reports that 52% support a candidate, they are relying on the LLN: with a large enough sample, the sample proportion will be close to the true population proportion. When a scientist repeats an experiment many times and averages the results, the LLN guarantees that the average converges to the expected value.

The LLN also explains why casinos always make money in the long run. Each individual bet is random, but over thousands of bets, the average outcome converges to the house edge. The randomness averages out.

The LLN is also the theoretical foundation of Monte Carlo methods (Chapter 8), where we estimate expectations by simulating many samples and averaging them. Every Monte Carlo estimator is a direct application of the LLN.

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

For the fair die, the true mean is `\mu = (1 + 2 + 3 + 4 + 5 + 6)/6 = 3.5`$. For the Bernoulli process with `p = 0.3`$, the true mean is `\mu = 0.3`$.

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

Watch the column `|M_n - \mu|`$ shrink as `n`$ grows. With only `10`$ rolls, the sample mean can be off by `0.3`$ or more. By `100{,}000`$ rolls, the deviation is typically under `0.002`$. By `1{,}000{,}000`$, it is essentially zero. This is the Law of Large Numbers in action.

The convergence is not perfectly monotonic. You can see that the die simulation at `n = 10000`$ has a slightly larger deviation than at `n = 100000`$. This is expected: the LLN guarantees convergence in the long run, but individual samples can fluctuate. The trend is what matters.

## A Practical Observation

Notice that the deviation shrinks roughly as `1/\sqrt{n}`$, not as `1/n`$. Going from `10`$ to `100`$ samples reduces the error by about a factor of `3`$ (roughly `\sqrt{10}`$), not a factor of `10`$. This slow convergence rate is a fundamental limitation of averaging. To halve the error, you need about `4`$ times as much data. We will see this `1/\sqrt{n}`$ rate appear again in the Monte Carlo chapter.

The `1/\sqrt{n}`$ rate is not a defect of a particular simulation; it is a consequence of the fact that `\mathrm{Var}(M_n) = \sigma^2 / n`$. Standard deviations scale as the square root of variance, so typical deviations of `M_n`$ from `\mu`$ scale as `\sigma / \sqrt{n}`$. This is a universal law for i.i.d. averages with finite variance.

## Problem Set

**Problem 6.1.** For the fair-die simulation, plot (on paper if necessary) the sample mean `M_n`$ as a function of `n`$. Add a horizontal line at `\mu = 3.5`$ and dashed lines at `\mu \pm \sigma/\sqrt{n}`$ using `\sigma^2 = 35/12`$. Describe how the sample mean fluctuates relative to those confidence bands.

**Problem 6.2 (Chebyshev bound).** For a fair die (`\mu = 3.5`$, `\sigma^2 = 35/12`$), use Chebyshev's inequality to bound `P(|M_{100} - 3.5| \geq 0.5)`$. Compare against the empirical frequency of this event by running the simulation `1000`$ times.

**Problem 6.3 (Different distributions).** Modify the example program to simulate the sample mean of an exponential distribution with rate `\lambda = 1`$. True mean is `1`$. Watch the sample mean converge to `1`$ as `n`$ grows. Because the exponential is more heavily right-skewed than the fair die, does convergence appear faster or slower?

**Problem 6.4 (The Cauchy failure).** Add a function that samples from a Cauchy distribution using inverse-CDF: draw `U`$ from `\text{Uniform}(0, 1)`$ and return `\tan(\pi(U - 0.5))`$. Compute the sample mean for `n = 10, 100, 1000, \ldots, 1{,}000{,}000`$ Cauchy samples. Does the sample mean converge, or does it keep jumping around? Explain what you observe in terms of the LLN's assumption of finite mean.

**Problem 6.5 (Empirical variance of M_n).** Run `1000`$ independent simulations of a sample mean of `100`$ die rolls. Record the `1000`$ sample means and compute their variance. Compare against the theoretical prediction `\mathrm{Var}(M_n) = \sigma^2/n = 35/(12 \cdot 100)`$.

**Problem 6.6 (Rate of convergence).** For the `\text{Bernoulli}(0.5)`$ coin, empirically estimate the rate at which `|M_n - 0.5|`$ shrinks. Fit a curve of the form `c/n^r`$ to the deviations at `n = 10, 100, 1000, \ldots`$, and estimate `r`$. Does your estimate come out close to `1/2`$ as the theory predicts?

**Problem 6.7 (Coding exercise).** Extend the example program to a function `lln-demo` that takes a sampling procedure (a thunk returning one sample) and a true mean, and prints a table of the same form as the example's die output. Test your function on the die, the Bernoulli, and the exponential.

**Problem 6.8 (Almost sure vs. in probability).** Describe informally what it would look like if the sample mean of a die converged in probability but not almost surely. Would the running sample mean have to keep spiking away from `3.5`$ forever, or only rarely? Even without a rigorous treatment, this thought experiment sharpens the distinction between the two modes of convergence.
