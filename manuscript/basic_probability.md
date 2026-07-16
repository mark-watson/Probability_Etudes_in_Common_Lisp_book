# Basic Probability {#basic_probability}

Probability theory begins with a simple but powerful idea: we can assign numbers between 0 and 1 to events to measure how likely they are. This chapter introduces the fundamental vocabulary of probability theory and the three axioms that everything else is built upon.

The example program for this chapter is in the file **01_basic_probability.lisp**.

## A Brief History

The mathematical study of probability grew out of very practical questions. In the sixteenth century the Italian polymath Gerolamo Cardano wrote *Liber de Ludo Aleae* (The Book on Games of Chance), analyzing dice games and card games with an eye toward gambling strategy. In 1654 the French mathematicians Blaise Pascal and Pierre de Fermat exchanged a famous series of letters about the "problem of points": how should the stakes be divided if a game of chance is interrupted before it is finished? Their correspondence is often cited as the birth of probability theory as a mathematical discipline.

Over the next three centuries, Jakob Bernoulli, Abraham de Moivre, Pierre-Simon Laplace, Carl Friedrich Gauss, Pafnuty Chebyshev, and many others developed probability into a sophisticated body of theory. But for a long time the field lacked rigorous foundations. Different authors used different definitions, and paradoxes appeared whenever the intuitive notion of "probability" was pushed hard enough.

In 1933 the Russian mathematician Andrey Kolmogorov published a short monograph, *Grundbegriffe der Wahrscheinlichkeitsrechnung* (Foundations of the Theory of Probability), that finally put probability on firm mathematical footing. He grounded it in measure theory and reduced its foundations to just three axioms. Everything we do in this book, and everything in modern probability theory, ultimately traces back to Kolmogorov's axioms.

## Interpretations of Probability

Before we plunge into the mathematics, it is worth pausing to ask: what does a probability actually mean? There are several answers, and they are not always in conflict.

The **classical interpretation**, associated with Laplace, says that if there are `n`$ equally likely outcomes and `k`$ of them are favorable to some event, then the probability of the event is `k/n`$. This is the interpretation used when we say a fair die has probability 1/6 of showing a 4. It works well for games of chance where symmetry gives us equally likely outcomes for free.

The **frequentist interpretation** says that the probability of an event is the long-run fraction of times it occurs in repeated trials. If we flip a coin ten thousand times and see 4,972 heads, we would estimate the probability of heads as roughly 0.4972. This interpretation grounds probability in observable frequencies and is the philosophical basis for much of classical statistics.

The **subjective (or Bayesian) interpretation** treats probability as a degree of belief. Under this view, it is meaningful to talk about the probability that a particular historical event happened, or the probability that a specific hypothesis is true, even if there is no notion of "repeated trials." Bayesians update their beliefs using Bayes' theorem as new evidence arrives.

The remarkable thing about the Kolmogorov axioms is that they are compatible with all three interpretations. The mathematics works the same way regardless of which interpretation you prefer, and each interpretation illuminates different applications.

## Sample Spaces and Events

A **random experiment** is any procedure whose outcome we cannot predict with certainty. Flipping a coin, rolling a die, and drawing a card from a deck are all random experiments.

The **sample space**, written `\Omega`$, is the set of all possible outcomes of a random experiment. When we roll a single die, the sample space is `\{1, 2, 3, 4, 5, 6\}`$. When we roll two dice, the sample space consists of all `36`$ ordered pairs `(i, j)`$ where `i`$ and `j`$ range from `1`$ to `6`$.

An **event** is any subset of the sample space. When rolling two dice, the event "the sum is 7" is the set of all outcomes where the two dice add up to `7`$. There are `6`$ such outcomes: `(1,6), (2,5), (3,4), (4,3), (5,2)`$, and `(6,1)`$.

Note that events are just sets of outcomes. The distinction between an outcome (a single element of `\Omega`$) and an event (a subset of `\Omega`$) is important. The outcome `(3,4)`$ is a single point of the sample space; the event "the first die is 3" is the six-element set `\{(3,1), (3,2), (3,3), (3,4), (3,5), (3,6)\}`$. A single outcome can also be viewed as a singleton event, but not every event corresponds to a single outcome.

In the example program, we build the sample space for two dice as a list of all 36 ordered pairs:

{lang="lisp",linenos=off}
~~~~~~~~
(defun make-two-dice-sample-space ()
  "Build the sample space Omega for rolling two distinguishable dice.
   Omega = { (i,j) : 1 <= i <= 6, 1 <= j <= 6 } has |Omega| = 36 outcomes."
  (loop for i from 1 to 6
        append (loop for j from 1 to 6 collect (list i j))))
~~~~~~~~

The events are defined by filtering the sample space with a predicate. For example, the event "sum equals 7" keeps only the outcomes where the two dice add up to 7:

{lang="lisp",linenos=off}
~~~~~~~~
(defun event-sum-is (omega target)
  "Event: outcomes whose dice sum equals TARGET."
  (remove-if-not (lambda (outcome) (= (reduce #'+ outcome) target)) omega))
~~~~~~~~

## Set Operations on Events

Because events are sets, we can combine them using set operations, and each operation has a natural probabilistic meaning:

- **Union** `A \cup B`$: the event that `A`$ or `B`$ (or both) occurs. In Lisp we compute it with `union`.
- **Intersection** `A \cap B`$: the event that both `A`$ and `B`$ occur. In Lisp we compute it with `intersection`.
- **Complement** `A^c`$ (or "not `A`$"): the event that `A`$ does not occur. Its outcomes are those in `\Omega`$ but not in `A`$.
- **Difference** `A \setminus B`$: the event that `A`$ occurs but `B`$ does not. It equals `A \cap B^c`$.

Two events are called **disjoint** (or mutually exclusive) if their intersection is empty. Disjoint events cannot both happen at the same trial. For example, "sum equals 7" and "sum equals 11" are disjoint events on the two-dice sample space.

**De Morgan's laws** connect complementation to union and intersection:

```$
(A \cup B)^c = A^c \cap B^c, \qquad (A \cap B)^c = A^c \cup B^c.
```

In words, "not (`A`$ or `B`$)" is the same as "not `A`$ and not `B`$," and "not (`A`$ and `B`$)" is the same as "not `A`$ or not `B`$." These laws are extraordinarily useful when you find it easier to reason about complements.

## The Kolmogorov Axioms

In 1933, the Russian mathematician Andrey Kolmogorov formulated three axioms that serve as the foundation of modern probability theory. A **probability measure** `P`$ is a function that assigns a number to each event, satisfying these three rules:

1. **Non-negativity**: `P(A) \geq 0`$ for every event `A`$. Probabilities are never negative.
2. **Normalization**: `P(\Omega) = 1`$. The probability of the entire sample space is `1`$, meaning that some outcome must occur.
3. **Countable additivity**: If `A_1, A_2, A_3, \ldots`$ are pairwise disjoint events (no two of them share an outcome), then `P(A_1 \cup A_2 \cup \cdots) = \sum_i P(A_i)`$.

These three axioms are remarkably compact, yet they generate the entire edifice of probability theory. Every theorem we will encounter in this book can be traced back to these three statements.

### Consequences of the Axioms

Several important properties follow immediately from the three axioms. Each one is worth internalizing because we will use them constantly.

**Probability of the empty set is zero.** Since `\Omega`$ and the empty set are disjoint and their union is `\Omega`$, by additivity `P(\Omega) + P(\emptyset) = P(\Omega)`$, so `P(\emptyset) = 0`$.

**Probabilities are bounded above by 1.** Since `A`$ and its complement are disjoint and their union is `\Omega`$, we have `P(A) + P(A^c) = 1`$. Non-negativity of `P(A^c)`$ forces `P(A) \leq 1`$. Combining with axiom 1, every probability lies in the interval `[0, 1]`$.

**Complement rule.** As already noted, `P(A^c) = 1 - P(A)`$. This is often the fastest way to compute a probability when the direct calculation is awkward.

**Monotonicity.** If `A \subseteq B`$, then `P(A) \leq P(B)`$. This says that adding outcomes to an event can only make it more likely, which matches intuition.

**Inclusion-exclusion for two events.** For any events `A`$ and `B`$:

```$
P(A \cup B) = P(A) + P(B) - P(A \cap B).
```

The term `P(A \cap B)`$ corrects for double-counting the outcomes that lie in both events. For three or more events, a longer inclusion-exclusion formula alternates signs; we will not need it in this book.

**Subadditivity (Boole's inequality).** For any events `A`$ and `B`$, whether disjoint or not:

```$
P(A \cup B) \leq P(A) + P(B).
```

More generally, the probability of a union is at most the sum of the individual probabilities. This is sometimes called the **union bound** and is extremely useful in probabilistic analysis of algorithms.

## The Classical Definition of Probability

When all outcomes in the sample space are **equally likely**, the probability of an event `A`$ reduces to a simple counting problem:

```$
P(A) = \frac{|A|}{|\Omega|}.
```

That is, the probability of event `A`$ is the number of favorable outcomes divided by the total number of outcomes. This is called the **classical definition of probability**, and it is what most people learn first.

For two fair dice, each of the 36 outcomes is equally likely, so we can compute probabilities by counting:

{lang="lisp",linenos=off}
~~~~~~~~
(defun probability-classical (event omega)
  "Compute P(event) under the equally-likely-outcomes model:
       P(A) = |A| / |Omega|."
  (/ (length event) (length omega)))
~~~~~~~~

The function returns an exact rational number. Lisp's rational arithmetic is a nice fit for probability calculations because many probabilities are fractions like 1/6 or 11/36.

### A Word About Counting

The classical formula reduces every problem to counting. Two combinatorial tools appear again and again:

**Permutations** count ordered arrangements. The number of ways to arrange `n`$ distinct objects in a row is `n!`$ (`n`$ factorial). The number of ordered arrangements of `k`$ objects chosen from `n`$ is `n! / (n - k)!`$.

**Combinations** count unordered selections. The number of ways to choose `k`$ objects from `n`$ distinct objects without regard to order is the binomial coefficient `\binom{n}{k} = \frac{n!}{k!\,(n - k)!}`$. We will meet the binomial coefficient again in the chapter on the binomial distribution.

For the two-dice sample space, no combinatorial machinery is needed: the sample space has exactly `6 \cdot 6 = 36`$ elements, one for each ordered pair. But in a card problem, or in a problem involving many coin flips, you almost always start by counting the sample space with permutations or combinations.

## The Complement Rule

A direct consequence of the axioms is the **complement rule**:

```$
P(A^c) = 1 - P(A).
```

This follows because an event `A`$ and its complement (everything in the sample space that is not in `A`$) partition the sample space. By the additivity axiom, `P(A) + P(A^c) = P(\Omega) = 1`$, so `P(A^c) = 1 - P(A)`$.

The complement rule is surprisingly useful. Sometimes it is easier to compute the probability that something does NOT happen, and then subtract from 1. For example, the probability of rolling "at least one 6" in two dice is easier to compute via its complement "no 6 appears at all":

{lang="lisp",linenos=off}
~~~~~~~~
(defun demonstrate-complement-rule (event omega)
  "The COMPLEMENT RULE is a consequence of the axioms:
       P(not A) = 1 - P(A)."
  (- 1 (probability-classical event omega)))
~~~~~~~~

## Running the Example

When you load the program, it builds the two-dice sample space and computes several probabilities:

```
Sample space Omega: two dice, |Omega| = 36
Axiom check P(Omega) = 1

Event 'sum = 7': 6 favorable outcome(s)
  P(sum = 7) = 1/6 = .167

Event 'at least one 6': 11 favorable outcome(s)
  P(at least one 6) = 11/36 = .306
  Complement P(no 6) = 1 - P(at least one 6) = .694

Event 'sum >= 10': 6 favorable outcome(s)
  P(sum >= 10) = 1/6 = .167
```

Notice that `P(\Omega) = 1`$, confirming the normalization axiom. The probability of rolling a sum of `7`$ is `1/6`$, which makes sense because there are `6`$ favorable outcomes out of `36`$. The complement rule gives us `P(\text{no } 6) = 1 - 11/36 = 25/36`$, which is approximately `0.694`$.

## Why This Matters

The concepts in this chapter may seem simple, but they are the foundation for everything that follows. Conditional probability, random variables, distributions, and the great theorems of probability all build on the ideas of sample spaces, events, and the Kolmogorov axioms. When you find yourself confused by a more advanced concept, returning to these foundations often clears things up.

## Problem Set

The problems below give you a chance to apply the ideas of this chapter. All of them can be solved on paper, but I encourage you to also verify your answers by extending the Common Lisp program. Building the sample space and filtering by a predicate is a very general recipe that works for many discrete probability problems.

**Problem 1.1.** For the two-dice sample space, compute `P(\text{sum is even})`$ and `P(\text{sum is odd})`$. Verify that they add to `1`$, as the axioms require.

**Problem 1.2.** Let `A`$ be the event "the first die is 3" and `B`$ be the event "the sum is 5." Compute `P(A), P(B), P(A \cap B)`$, and `P(A \cup B)`$ directly. Then verify the inclusion-exclusion formula `P(A \cup B) = P(A) + P(B) - P(A \cap B)`$.

**Problem 1.3.** Compute the probability that at least one die shows an even number. Do this two ways: (a) directly by counting outcomes, and (b) via the complement rule using the event "both dice are odd." Confirm the answers agree.

**Problem 1.4.** Extend the program to three dice. The new sample space has `6^3 = 216`$ outcomes. Compute the probability that the sum is exactly `10`$. (For those who want to check by hand, the answer is `27/216 = 1/8`$.)

**Problem 1.5.** Draw a card from a standard 52-card deck. Let `A`$ be the event "the card is a heart" and `B`$ be the event "the card is a face card (jack, queen, or king)." Assuming all `52`$ cards are equally likely, compute `P(A), P(B), P(A \cap B)`$, and `P(A \cup B)`$. Are `A`$ and `B`$ disjoint? Explain.

**Problem 1.6.** A fair coin is flipped `5`$ times. Using the classical definition, what is the probability of getting exactly `3`$ heads? Hint: the sample space has `2^5 = 32`$ outcomes, and the number of ways to arrange `3`$ heads among `5`$ flips is the binomial coefficient `\binom{5}{3}`$. Now compute the probability of "at least 3 heads" and the probability of "at most 3 heads," and use these to check the complement rule.

**Problem 1.7 (De Morgan in action).** For any two events `A`$ and `B`$ on any sample space, prove that `P(A^c \cap B^c) = 1 - P(A \cup B)`$. Then use this identity to give an alternative derivation of the inclusion-exclusion formula.

**Problem 1.8 (Coding exercise).** Modify the example program to accept an arbitrary predicate and print (a) the number of favorable outcomes, (b) the probability as an exact rational, and (c) the probability as a decimal. Then use this generalized reporter to compute the probability that the two dice differ by exactly `2`$.
