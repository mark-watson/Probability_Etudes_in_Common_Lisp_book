# Discrete Random Variables {#discrete_rv}

So far we have been working with events as sets of outcomes. A **random variable** gives us a more convenient way to talk about randomness by assigning a number to each outcome. This lets us work with averages, spreads, and distributions instead of individual events.

The example program for this chapter is in the file **03_discrete_random_variables.lisp**.

## What Is a Random Variable?

A **random variable** X is a function from the sample space to the real numbers. It assigns a numerical value to each outcome. For example, if we roll two dice, we can define X to be the sum of the two dice. Then X maps the outcome (3, 4) to the value 7, the outcome (1, 1) to the value 2, and so on.

A **discrete** random variable takes values in a countable set, often a finite set of integers. The sum of two dice is a discrete random variable that takes values from 2 to 12.

## The Probability Mass Function

The **probability mass function** (PMF) of a discrete random variable X gives the probability that X equals a specific value x:

    p(x) = P(X = x) = P({ outcomes where X equals x })

The PMF satisfies two properties: every value p(x) is non-negative, and the sum of all p(x) values equals 1. These are the discrete analogues of the Kolmogorov axioms.

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

    F(x) = P(X <= x) = sum of p(t) for all t <= x

The CDF is a non-decreasing function: as x increases, F(x) can only go up or stay the same. It starts at 0 (for values below the smallest possible outcome) and ends at 1 (for values above the largest possible outcome).

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

    E[X] = sum of x * p(x) over all x

The expected value represents the long-run average if we repeated the experiment many times. If we rolled two dice a million times and averaged all the sums, we would get a number very close to E[X] = 7.

{lang="lisp",linenos=off}
~~~~~~~~
(defun expectation (p)
  "E[X] = sum of x * p_X(x), the mean of the distribution."
  (reduce #'+ (mapcar (lambda (entry) (* (car entry) (cdr entry)))
                      (pmf-table p))))
~~~~~~~~

A key property of expectation is **linearity**. For any constants a and b:

    E[aX + b] = a * E[X] + b

And for any two random variables X and Y, even if they are dependent:

    E[X + Y] = E[X] + E[Y]

Linearity of expectation is one of the most useful tools in probability. It lets us compute the expected value of a sum without knowing anything about the relationship between the variables.

## Variance and Standard Deviation

The **variance** measures how spread out a distribution is around its mean:

    Var(X) = E[(X - mean)^2]

There is a computational shortcut that is often easier to use:

    Var(X) = E[X^2] - (E[X])^2

This identity comes from expanding (X - mean)^2 and using linearity of expectation. The program uses this shortcut:

{lang="lisp",linenos=off}
~~~~~~~~
(defun variance (p)
  "Var(X) = E[X^2] - (E[X])^2."
  (let ((ex (expectation p)))
    (- (expectation-of-square p) (* ex ex))))
~~~~~~~~

The **standard deviation** is the square root of the variance. It has the same units as the random variable itself, which makes it more interpretable. A standard deviation of 2.4 for the sum of two dice means that a typical roll deviates from the mean of 7 by about 2.4 points.

## Three Example Distributions

The program demonstrates three distributions. The first is a **loaded die** where the probability of rolling a 6 is 1/2 and the remaining probabilities are 1/10 each. The mean is 4.5 (higher than the fair die's 3.5) because the loaded die favors high values.

The second is a **Bernoulli distribution** with p = 1/2. A Bernoulli random variable takes only two values: 1 (success) with probability p, and 0 (failure) with probability 1-p. Its mean is p and its variance is p(1-p). The Bernoulli distribution is the simplest discrete distribution and the building block for the binomial distribution we will study in the next chapter.

The third is the **sum of two fair dice**. This is the classic triangular distribution: the probability of rolling a 7 is highest (6/36 = 1/6) because there are more ways to make 7 than any other sum, while the probability of rolling a 2 or 12 is lowest (1/36 each) because there is only one way to make each.

## Running the Example

```
=== Loaded Die ===
Distribution of LoadedDie:
  Normalization check: sum p(x) = 1
  x=1 : P(X=x)=1/10  P(X<=x)=1/10
  x=2 : P(X=x)=1/10  P(X<=x)=1/5
  ...
  x=6 : P(X=x)=1/2   P(X<=x)=1
  E[LoadedDie]   = 9/2 =  4.5
  Var(LoadedDie) = 13/4 = 3.25
  sigma   =  1.8

=== Sum of Two Fair Dice ===
  E[DiceSum]   = 7 =  7.0
  Var(DiceSum) = 35/6 = 5.83
  sigma   = 2.42
```

The normalization check confirms that the PMF sums to 1. The expected value of the dice sum is exactly 7, and the variance is 35/6, which is about 5.83. These exact rational results come from Lisp's built-in rational arithmetic.

## Why This Matters

Expected value and variance are the two most important summary statistics for any random variable. Throughout the rest of this book, we will compute means and variances for every distribution we encounter. The formulas E[X] = sum of x * p(x) and Var(X) = E[X^2] - (E[X])^2 are tools you will use again and again.
