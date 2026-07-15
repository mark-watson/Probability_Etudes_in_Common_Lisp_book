# Basic Probability {#basic_probability}

Probability theory begins with a simple but powerful idea: we can assign numbers between 0 and 1 to events to measure how likely they are. This chapter introduces the fundamental vocabulary of probability theory and the three axioms that everything else is built upon.

The example program for this chapter is in the file **01_basic_probability.lisp**.

## Sample Spaces and Events

A **random experiment** is any procedure whose outcome we cannot predict with certainty. Flipping a coin, rolling a die, and drawing a card from a deck are all random experiments.

The **sample space**, written as the Greek letter Omega, is the set of all possible outcomes of a random experiment. When we roll a single die, the sample space is {1, 2, 3, 4, 5, 6}. When we roll two dice, the sample space consists of all 36 ordered pairs (i, j) where i and j range from 1 to 6.

An **event** is any subset of the sample space. When rolling two dice, the event "the sum is 7" is the set of all outcomes where the two dice add up to 7. There are 6 such outcomes: (1,6), (2,5), (3,4), (4,3), (5,2), and (6,1).

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

## The Kolmogorov Axioms

In 1933, the Russian mathematician Andrey Kolmogorov formulated three axioms that serve as the foundation of modern probability theory. A **probability measure** P is a function that assigns a number to each event, satisfying these three rules:

1. **Non-negativity**: P(A) >= 0 for every event A. Probabilities are never negative.
2. **Normalization**: P(Omega) = 1. The probability of the entire sample space is 1, meaning that some outcome must occur.
3. **Countable additivity**: If A1, A2, A3, ... are pairwise disjoint events (no two of them share an outcome), then P(A1 union A2 union ...) = P(A1) + P(A2) + ...

These three axioms are remarkably compact, yet they generate the entire edifice of probability theory. Every theorem we will encounter in this book can be traced back to these three statements.

## The Classical Definition of Probability

When all outcomes in the sample space are **equally likely**, the probability of an event A reduces to a simple counting problem:

    P(A) = |A| / |Omega|

That is, the probability of event A is the number of favorable outcomes divided by the total number of outcomes. This is called the **classical definition of probability**, and it is what most people learn first.

For two fair dice, each of the 36 outcomes is equally likely, so we can compute probabilities by counting:

{lang="lisp",linenos=off}
~~~~~~~~
(defun probability-classical (event omega)
  "Compute P(event) under the equally-likely-outcomes model:
       P(A) = |A| / |Omega|."
  (/ (length event) (length omega)))
~~~~~~~~

The function returns an exact rational number. Lisp's rational arithmetic is a nice fit for probability calculations because many probabilities are fractions like 1/6 or 11/36.

## The Complement Rule

A direct consequence of the axioms is the **complement rule**:

    P(not A) = 1 - P(A)

This follows because an event A and its complement (everything in the sample space that is not in A) partition the sample space. By the additivity axiom, P(A) + P(not A) = P(Omega) = 1, so P(not A) = 1 - P(A).

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

Notice that P(Omega) = 1, confirming the normalization axiom. The probability of rolling a sum of 7 is 1/6, which makes sense because there are 6 favorable outcomes out of 36. The complement rule gives us P(no 6) = 1 - 11/36 = 25/36, which is approximately 0.694.

## Why This Matters

The concepts in this chapter may seem simple, but they are the foundation for everything that follows. Conditional probability, random variables, distributions, and the great theorems of probability all build on the ideas of sample spaces, events, and the Kolmogorov axioms. When you find yourself confused by a more advanced concept, returning to these foundations often clears things up.
