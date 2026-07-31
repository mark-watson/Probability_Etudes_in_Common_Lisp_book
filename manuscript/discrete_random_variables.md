# Discrete Random Variables {#discrete_rv}

So far we have been working with events as sets of outcomes. A **random variable** gives us a more convenient way to talk about randomness by assigning a number to each outcome. This lets us work with averages, spreads, and distributions instead of individual events.

The example program for this chapter is in the file **03_discrete_random_variables.lisp**.

## From Events to Random Variables

Working with events directly gets cumbersome once questions become quantitative. Consider rolling two dice and asking about the sum. We could describe the events "sum equals 2," "sum equals 3," ..., "sum equals 12" one at a time and compute their probabilities separately. But it is much cleaner to introduce a single object, the sum-of-dice random variable, and treat all these events as questions about it.

Random variables give us the language for summary statistics like mean and variance, for comparing distributions, for talking about correlations between quantities, and for stating limit theorems. They are the mathematical objects that a probability calculation is usually really about.

## What Is a Random Variable?

A **random variable** `X`$ is a function from the sample space to the real numbers. It assigns a numerical value to each outcome. For example, if we roll two dice, we can define `X`$ to be the sum of the two dice. Then `X`$ maps the outcome `(3, 4)`$ to the value `7`$, the outcome `(1, 1)`$ to the value `2`$, and so on.

The name "random variable" is slightly misleading. A random variable is not itself random; it is a deterministic function. The randomness comes from the underlying experiment. Once the outcome `\omega \in \Omega`$ is fixed, `X(\omega)`$ is a definite number. This distinction matters when we discuss expectation and variance: those are properties of the function `X`$ together with the probability measure `P`$, not of any particular value.

The **support** of a random variable is the set of values it can actually take. For the sum of two dice, the support is `\{2, 3, 4, \ldots, 12\}`$. For a coin flip encoded as `0`$ or `1`$, the support is `\{0, 1\}`$. A **discrete** random variable takes values in a countable set, often a finite set of integers.

A technical requirement hides behind the word "function." For `P(X = x)`$ and `P(X \leq x)`$ to make sense, the sets `\{\omega : X(\omega) = x\}`$ and `\{\omega : X(\omega) \leq x\}`$ must be events, that is, members of the `\sigma`$-algebra `\mathcal{F}`$ from the previous chapter. A function `X : \Omega \to \mathbb{R}`$ with this property is called **measurable**, and only measurable functions qualify as random variables. In the discrete case, where `\mathcal{F}`$ is the full power set, every function is measurable and the requirement is automatic, so we will not raise it again until continuous distributions make it matter.

Any function of a random variable is again a random variable. If `X`$ is a random variable and `g`$ is a function, then `g(X)`$ assigns `g(X(\omega))`$ to each outcome `\omega`$. So `X^2, |X|, e^{X}`$, and `X \bmod 2`$ are all random variables built from `X`$.

## Indicator Random Variables

The simplest and most useful random variable is the **indicator** of an event. Given an event `A`$, its indicator `I_A`$ is defined as:

```$
I_A(\omega) = \begin{cases} 1 & \text{if } \omega \in A, \\ 0 & \text{otherwise.} \end{cases}
```

Indicator random variables are the bridge between events and expectations. Note that `E[I_A] = 1 \cdot P(A) + 0 \cdot P(A^c) = P(A)`$. Every probability can be written as an expectation of an indicator. This turns out to be extraordinarily useful when combined with linearity of expectation, as we will see below.

## The Probability Mass Function

The **probability mass function** (PMF) of a discrete random variable `X`$ gives the probability that `X`$ equals a specific value `x`$:

```$
p(x) = P(X = x) = P(\{ \omega : X(\omega) = x \}).
```

The PMF satisfies two properties: every value `p(x)`$ is non-negative, and the sum of all `p(x)`$ values equals `1`$. These are the discrete analogues of the Kolmogorov axioms.

In the example program, we represent a PMF as an alist mapping values to probabilities:

{lang="lisp",linenos=off}
~~~~~~~~
(defstruct pmf
  "A discrete probability mass function represented as an alist
   mapping value -> probability."
  table)
~~~~~~~~

The function **pmf-total-probability** verifies the normalization axiom by summing all the probabilities. For a valid PMF, this sum must equal 1.

## The Cumulative Distribution Function

The **cumulative distribution function** (CDF) accumulates probability up to a given value:

```$
F(x) = P(X \leq x) = \sum_{t \leq x} p(t).
```

The CDF is a non-decreasing function: as `x`$ increases, `F(x)`$ can only go up or stay the same. It starts at `0`$ (for values below the smallest possible outcome) and ends at `1`$ (for values above the largest possible outcome). For discrete random variables, the CDF is a step function that jumps at each value in the support of `X`$ by an amount equal to `p(x)`$.

The CDF is more universal than the PMF: it is defined for any real-valued random variable, discrete or continuous. Two random variables that share the same CDF share the same distribution. When we talk about "the distribution of `X`$," the CDF is the most general representation.

{lang="lisp",linenos=off}
~~~~~~~~
(defun cdf (p x)
  "Cumulative distribution function F_X(x) = P(X <= x)."
  (reduce #'+ (mapcar (lambda (entry)
                        (if (<= (car entry) x) (cdr entry) 0))
                      (pmf-table p))))
~~~~~~~~

## Expected Value

The **expected value** (or **mean**) of a random variable is the probability-weighted average of all its possible values:

```$
E[X] = \sum_{x} x \, p(x).
```

The expected value represents the long-run average if we repeated the experiment many times. If we rolled two dice a million times and averaged all the sums, we would get a number very close to `E[X] = 7`$.

{lang="lisp",linenos=off}
~~~~~~~~
(defun expectation (p)
  "E[X] = sum of x * p_X(x), the mean of the distribution."
  (reduce #'+ (mapcar (lambda (entry) (* (car entry) (cdr entry)))
                      (pmf-table p))))
~~~~~~~~

More generally, for any function `g`$ of a random variable `X`$, the **law of the unconscious statistician** says:

```$
E[g(X)] = \sum_{x} g(x) \, p(x).
```

This lets us compute the expected value of any function of `X`$ without needing to derive the distribution of `g(X)`$ first.

### Linearity of Expectation

A key property of expectation is **linearity**. For any constants `a`$ and `b`$:

```$
E[aX + b] = a\, E[X] + b.
```

And for any two random variables `X`$ and `Y`$, even if they are dependent:

```$
E[X + Y] = E[X] + E[Y].
```

Linearity of expectation is one of the most useful tools in probability. It lets us compute the expected value of a sum without knowing anything about the relationship between the variables. This holds even if `X`$ and `Y`$ are strongly correlated, which is often surprising to newcomers.

One striking application uses indicator random variables. Suppose we want the expected number of pairs among `n`$ coin flips that both show heads. Let `I_{ij}`$ be the indicator that flips `i`$ and `j`$ both show heads. The count of matching heads pairs is `\sum_{i < j} I_{ij}`$. By linearity, the expected count is `\binom{n}{2} \cdot 1/4`$, regardless of whether we treat the flips as independent or not. Indicator variables plus linearity turn many counting-in-expectation problems into one-line calculations.

## Variance and Standard Deviation

The **variance** measures how spread out a distribution is around its mean:

```$
\mathrm{Var}(X) = E[(X - E[X])^2].
```

There is a computational shortcut that is often easier to use:

```$
\mathrm{Var}(X) = E[X^2] - (E[X])^2.
```

This identity comes from expanding `(X - \mu)^2`$ and using linearity of expectation. The program uses this shortcut:

{lang="lisp",linenos=off}
~~~~~~~~
(defun variance (p)
  "Var(X) = E[X^2] - (E[X])^2."
  (let ((ex (expectation p)))
    (- (expectation-of-square p) (* ex ex))))
~~~~~~~~

The **standard deviation** `\sigma`$ is the square root of the variance. It has the same units as the random variable itself, which makes it more interpretable. A standard deviation of `2.4`$ for the sum of two dice means that a typical roll deviates from the mean of `7`$ by about `2.4`$ points.

Unlike expectation, variance is **not** linear in general. For a constant `c`$:

```$
\mathrm{Var}(cX) = c^2 \, \mathrm{Var}(X), \qquad \mathrm{Var}(X + c) = \mathrm{Var}(X).
```

For a sum of two random variables:

```$
\mathrm{Var}(X + Y) = \mathrm{Var}(X) + \mathrm{Var}(Y) + 2 \, \mathrm{Cov}(X, Y),
```

where `\mathrm{Cov}(X, Y) = E[(X - E[X])(Y - E[Y])]`$ is the **covariance**. When `X`$ and `Y`$ are independent, the covariance is zero and `\mathrm{Var}(X + Y) = \mathrm{Var}(X) + \mathrm{Var}(Y)`$. This additivity of variance for independent variables is one of the reasons independence is so useful.

The covariance carries the product of the units of `X`$ and `Y`$, which makes its numerical value hard to interpret on its own. Dividing by the two standard deviations strips the units and yields the **correlation coefficient**:

```$
\rho(X, Y) = \frac{\mathrm{Cov}(X, Y)}{\sigma_X \, \sigma_Y}.
```

The Cauchy-Schwarz inequality forces `-1 \leq \rho \leq 1`$, with `\rho = \pm 1`$ exactly when `Y`$ is an affine function of `X`$. A correlation of `0`$ means `X`$ and `Y`$ are **uncorrelated**. Independence implies zero correlation, but the converse fails: uncorrelated variables can still be dependent, because `\rho`$ detects only *linear* association.

For a sum of many variables the two-term rule generalizes to

```$
\mathrm{Var}\!\left(\sum_{i=1}^{n} X_i\right) = \sum_{i=1}^{n} \mathrm{Var}(X_i) + 2 \sum_{i < j} \mathrm{Cov}(X_i, X_j).
```

When the variables are pairwise uncorrelated, every cross term drops out and the variance of the sum is just the sum of the variances. This is the fact that makes the variance of a sample mean shrink like `1/n`$, the engine of the Law of Large Numbers two chapters from now.

## Standardization

Given a random variable `X`$ with mean `\mu`$ and standard deviation `\sigma`$, its **standardized** form is:

```$
Z = \frac{X - \mu}{\sigma}.
```

By construction, `Z`$ has mean `0`$ and variance `1`$. Standardization lets us compare distributions on a common scale and is essential for stating limit theorems like the Central Limit Theorem in a clean form.

## Moments and the Moment Generating Function

The mean and variance are the first two members of a larger family. The `k`$-th **moment** of `X`$ is `E[X^k]`$, and the `k`$-th **central moment** is `E[(X - \mu)^k]`$. The mean is the first moment; the variance is the second central moment. The third and fourth central moments, once standardized, describe the shape of a distribution beyond its location and spread:

```$
\text{skewness} = \frac{E[(X-\mu)^3]}{\sigma^3}, \qquad \text{kurtosis} = \frac{E[(X-\mu)^4]}{\sigma^4}.
```

Skewness measures asymmetry: a positive value signals a long right tail. Kurtosis measures how heavy the tails are relative to a normal distribution, whose kurtosis is `3`$.

All of these moments are packaged together by the **moment generating function** (MGF):

```$
M_X(t) = E[e^{tX}] = \sum_{x} e^{tx} \, p(x),
```

defined for those `t`$ where the sum is finite. The name comes from a small miracle: differentiating `M_X`$ at `t = 0`$ reproduces the moments one at a time,

```$
M_X^{(k)}(0) = E[X^k].
```

This follows from expanding `e^{tX} = \sum_k (tX)^k / k!`$ and taking the expectation term by term. The MGF has two properties that make it a workhorse of the later chapters. First, it turns sums of *independent* variables into products: if `X`$ and `Y`$ are independent then `M_{X+Y}(t) = M_X(t)\,M_Y(t)`$, because the expectation of a product of independent quantities factors. Second, it determines the distribution uniquely: two random variables with the same MGF on an interval around `0`$ have the same distribution. Together these facts let us identify the distribution of a sum by multiplying MGFs, and they underlie one standard proof of the Central Limit Theorem. A close relative, the **characteristic function** `\varphi_X(t) = E[e^{itX}]`$, exists for *every* distribution even when the MGF does not, and is the tool of choice in the rigorous theory.

## Concentration Inequalities

Given only the mean and variance of a random variable, we can already bound how much probability can lie far from the mean. Two classic bounds appear again and again.

**Markov's inequality**: for a non-negative random variable `X`$ and any `a > 0`$:

```$
P(X \geq a) \leq \frac{E[X]}{a}.
```

**Chebyshev's inequality**: for any random variable `X`$ with mean `\mu`$ and finite variance `\sigma^2`$ and any `k > 0`$:

```$
P(|X - \mu| \geq k \sigma) \leq \frac{1}{k^2}.
```

Chebyshev's inequality follows from Markov's inequality applied to the non-negative random variable `(X - \mu)^2`$. It says that the probability of being more than `k`$ standard deviations from the mean is at most `1/k^2`$, no matter what the distribution is. This is the key ingredient in the proof of the Weak Law of Large Numbers, which we will meet in a later chapter.

These inequalities are usually quite loose: for a specific distribution you can often do much better. But their power is that they work for **every** distribution with finite mean or variance.

## Three Example Distributions

The program demonstrates three distributions. The first is a **loaded die** where the probability of rolling a `6`$ is `1/2`$ and the remaining probabilities are `1/10`$ each. The mean is `4.5`$ (higher than the fair die's `3.5`$) because the loaded die favors high values.

The second is a **Bernoulli distribution** with `p = 1/2`$. A Bernoulli random variable takes only two values: `1`$ (success) with probability `p`$, and `0`$ (failure) with probability `1 - p`$. Its mean is `p`$ and its variance is `p(1-p)`$. The Bernoulli distribution is the simplest discrete distribution and the building block for the binomial distribution we will study in the next chapter.

The third is the **sum of two fair dice**. This is the classic triangular distribution: the probability of rolling a `7`$ is highest (`6/36 = 1/6`$) because there are more ways to make `7`$ than any other sum, while the probability of rolling a `2`$ or `12`$ is lowest (`1/36`$ each) because there is only one way to make each.

## Additional Computations in the Program

The program does more than tabulate PMFs. It also demonstrates four ideas from earlier in this chapter, each in a few lines.

**Convolution.** The PMF of a sum of independent variables is the convolution of their PMFs, `p_{X+Y}(s) = \sum_{x} p_X(x)\, p_Y(s - x)`$. The function `convolve-pmfs` rebuilds the two-dice distribution from two single dice, and its mean (`7`$) and variance (`35/6`$) match the direct computation. This is the coding exercise from Problem 3.8.

**Standardization.** `standardize` maps each value `x`$ to `(x - \mu)/\sigma`$, producing the PMF of `Z`$. The program confirms that `Z`$ has mean `0`$ and variance `1`$, as Problem 3.7 asks.

**Chebyshev's inequality.** For the dice sum the program prints the actual `P(|X - \mu| \geq k\sigma)`$ next to the distribution-free bound `1/k^2`$, showing the bound holds, and is loose, at every `k`$.

**Covariance and correlation.** Using a joint PMF over value pairs, the program computes `\mathrm{Cov}(X, Y)`$ and `\rho(X, Y)`$ for two independent dice (both `0`$) and for a single die paired with the two-dice sum, where the covariance is `\mathrm{Var}(X) = 35/12`$ and the correlation is `\rho = 1/\sqrt{2} \approx 0.707`$.

## Running the Example

```
=== Loaded Die ===
Distribution of LoadedDie:
  Normalization check Σ p(x) = 1
  x=1 : P(X=x)=1/10  P(X<=x)=1/10
  x=2 : P(X=x)=1/10  P(X<=x)=1/5
  ...
  x=6 : P(X=x)=1/2  P(X<=x)=1
  E[LoadedDie]   = 9/2 =  4.5
  Var(LoadedDie) = 13/4 = 3.25
  sigma   =  1.8

=== Sum of Two Fair Dice ===
  E[DiceSum]   = 7 =  7.0
  Var(DiceSum) = 35/6 = 5.83
  sigma   = 2.42

=== Convolution: single die + single die ===
  E[die+die]   = 7 =  7.0 (matches DiceSum mean 7)
  Var(die+die) = 35/6 = 5.83 (matches DiceSum var 35/6)

=== Standardization: Z = (X - mu)/sigma ===
  mean(Z) =  0.000 (should be 0)
  var(Z)  =  1.000 (should be 1)

=== Chebyshev's Inequality (dice sum) ===
  k     actual P(|X-mu|>=k sigma)   bound 1/k^2
  1.0               0.3333            1.0000
  1.5               0.1667            0.4444
  2.0               0.0556            0.2500
  2.5               0.0000            0.1600

=== Covariance and Correlation ===
  Two independent dice (X, Y):
    Cov =  0.000  rho =  0.000 (independent => 0)
  A die and the two-dice sum (X, X+Y):
    Cov =  2.917  rho =  0.707 (rho = 1/sqrt(2) ~ 0.707)
```

The normalization check confirms that the PMF sums to `1`$. The expected value of the dice sum is exactly `7`$, and the variance is `35/6`$, which is about `5.83`$. These exact rational results come from Lisp's built-in rational arithmetic.

Notice that the mean of the sum of two dice is `7`$, which equals `3.5 + 3.5`$, the sum of the means. This is linearity of expectation in action. The variance is `35/6`$, which equals `35/12 + 35/12`$, the sum of the variances of each fair die (since the two rolls are independent). This is the additivity of variance for independent random variables.

## Why This Matters

Expected value and variance are the two most important summary statistics for any random variable. Throughout the rest of this book, we will compute means and variances for every distribution we encounter. The formulas `E[X] = \sum_x x \, p(x)`$ and `\mathrm{Var}(X) = E[X^2] - (E[X])^2`$ are tools you will use again and again. Linearity of expectation and Chebyshev's inequality will show up in the proofs of the Law of Large Numbers and the Central Limit Theorem in later chapters.

## Problem Set

**Problem 3.1.** Compute `E[X]`$ and `\mathrm{Var}(X)`$ for a fair six-sided die directly from the definition. Confirm your answer against the program's output for the sum of two dice (using linearity of expectation and additivity of variance for independent variables).

**Problem 3.2.** For the loaded die in the example (`P(X = 6) = 1/2`$, all other faces `1/10`$), what is `P(X \geq 5)`$? What is `E[X^2]`$? Verify the variance formula `\mathrm{Var}(X) = E[X^2] - (E[X])^2`$ using your computed values.

**Problem 3.3 (Linearity of expectation).** Consider rolling three fair dice. Let `X, Y, Z`$ be the three roll values, and let `S = X + Y + Z`$. Use linearity of expectation to compute `E[S]`$. Now use additivity of variance for independent variables to compute `\mathrm{Var}(S)`$. Verify your answers by extending the example program to enumerate the 216-outcome sample space and directly compute the mean and variance of `S`$.

**Problem 3.4 (Indicator random variables).** A word is chosen uniformly at random from a 26-letter alphabet, one letter at a time, for a word of length `n`$. Let `X`$ be the number of vowels in the word (a, e, i, o, u count as vowels). Using indicator random variables and linearity of expectation, compute `E[X]`$ for a word of length `10`$.

**Problem 3.5 (Chebyshev's inequality).** For a distribution with mean `100`$ and standard deviation `10`$, use Chebyshev's inequality to bound the probability that `X`$ lies outside the interval `[70, 130]`$. Then use the same inequality to bound the probability that `X`$ lies outside `[80, 120]`$. Comment on how the bound tightens as the interval widens.

**Problem 3.6 (Comparing distributions).** A game offers two prizes: Game A pays `1`$ dollar with probability `1/2`$ and pays `0`$ dollars with probability `1/2`$. Game B pays `100`$ dollars with probability `1/200`$ and pays `0`$ dollars with probability `199/200`$. Compute the mean and variance for each game. Which game has the same expected payout? Which has more variance? Which would you rather play if you could play it just once? If you could play it a million times?

**Problem 3.7 (Coding exercise).** Extend the PMF library with a function `standardize` that takes a PMF for `X`$ and returns a PMF for `Z = (X - \mu)/\sigma`$. Verify numerically that the new PMF has mean `0`$ and variance `1`$.

**Problem 3.8 (Sum of two independent PMFs).** Add a function `convolve-pmfs` to the example program that takes two PMFs (representing independent random variables `X`$ and `Y`$) and returns the PMF of the sum `X + Y`$. Use it to construct the PMF of the sum of two fair dice starting from the PMF of a single die. Confirm that the resulting mean and variance match the direct computation.
