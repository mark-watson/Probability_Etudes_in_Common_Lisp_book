# Conditional Probability {#conditional_probability}

In the previous chapter we learned to compute probabilities by counting outcomes. But what happens when we gain partial information about an experiment? If we know that one event has occurred, does that change the probability of another event? The answer is yes, and the tool for reasoning about this is **conditional probability**.

The example program for this chapter is in the file **02_conditional_probability.lisp**.

## Why Conditional Probability?

Almost every real-world probability question is really a conditional question. A doctor does not ask "what is the probability this patient has diabetes?" in the abstract. She asks "what is the probability this patient has diabetes given the observed symptoms and lab results?" A jury does not evaluate whether a defendant is guilty in a vacuum. It evaluates guilt given the evidence presented at trial. A weather forecaster does not report the marginal probability of rain averaged over all possible tomorrows; she reports the probability of rain given the current radar image and pressure readings.

Unconditional (or **marginal**) probabilities are useful as inputs, but almost all interesting reasoning is conditional. Learning to work fluently with conditional probability is the single most important skill in this book.

## The Definition of Conditional Probability

The **conditional probability** of event A given event B, written P(A | B), measures how likely A is once we know that B has occurred. The formula is:

    P(A | B) = P(A and B) / P(B)

The intuition is that once we know B happened, B becomes our new sample space. We ask: what fraction of B also contains A? The numerator is the probability that both A and B occur, and the denominator is the probability that B occurs at all.

Note that P(A | B) is defined only when P(B) > 0. Conditioning on an event of probability zero is not meaningful in the discrete setting. There is a more elaborate theory of conditional expectation that handles the continuous case, but we will not need it.

In the example program, we compute conditional probability by taking the intersection of the two events and dividing by the size of the conditioning event:

{lang="lisp",linenos=off}
~~~~~~~~
(defun conditional-probability (a b omega)
  "P(A | B) = P(A intersection B) / P(B).
   The intersection A intersection B is the set of outcomes in BOTH events."
  (let ((intersection (intersection a b :test #'equal)))
    (/ (length intersection) (length b))))
~~~~~~~~

### Conditional Probability as a Probability Measure

An important sanity check: P(. | B) is itself a bona-fide probability measure. Fixing B and letting A vary, the function P(A | B) satisfies all three Kolmogorov axioms on the sample space B. Non-negativity is immediate; normalization holds because P(B | B) = P(B) / P(B) = 1; and additivity carries over from the additivity of P. This means every theorem we know about probability measures applies equally well to conditional probabilities. In particular, the complement rule P(not A | B) = 1 - P(A | B) still holds.

## The Multiplication Rule and the Chain Rule

Rearranging the definition of conditional probability gives the **multiplication rule**:

    P(A and B) = P(B) * P(A | B) = P(A) * P(B | A)

This is often the easiest way to compute the joint probability of two events. Rather than dealing with the intersection directly, you compute one probability and then a conditional probability.

Extending this to more events gives the **chain rule of probability**:

    P(A1 and A2 and ... and An) = P(A1) * P(A2 | A1) * P(A3 | A1 and A2) * ... * P(An | A1 and ... and An-1)

The chain rule is the backbone of most probabilistic model-building. Any complicated joint distribution can be factored into a product of simpler conditional distributions. Bayesian networks and hidden Markov models rely on this idea directly.

## Independence

Two events A and B are **independent** if knowing one tells you nothing about the other. Formally, A and B are independent if and only if:

    P(A and B) = P(A) * P(B)

Equivalently, P(A | B) = P(A) when A and B are independent. The condition does not change the probability.

Independence is one of the most important concepts in probability. Many powerful results, such as the Law of Large Numbers and the Central Limit Theorem, require independence as a hypothesis. Whenever you assume events are independent, you are making a real assumption about the world; whenever a proof assumes independence, that assumption is doing real work.

### Independence Versus Disjoint

A common source of confusion is the difference between independent events and disjoint (mutually exclusive) events. These concepts are almost opposites.

**Disjoint** events cannot happen together: P(A and B) = 0. If A occurs, then B cannot occur. So if we learn that A occurred, we can be certain that B did not, meaning that P(B | A) = 0, not P(B). Disjoint events are highly dependent.

**Independent** events can happen together, and their joint probability is exactly the product P(A) * P(B). Learning that A occurred does not change our estimate of B.

Disjoint events with positive probability are never independent, and independent events with positive probability are never disjoint. This is worth pausing to internalize; many mistakes in probability come from confusing these two concepts.

In our two-dice example, the event "first die is 1" and the event "second die is 1" are independent because the two dice do not affect each other. But the event "first die is 1" and the event "sum is 7" are also independent, which is less obvious. Let us see why: if the first die is 1, the second die must be 6 for the sum to be 7, and that happens with probability 1/6, which is the same as the unconditional probability P(sum = 7) = 6/36 = 1/6.

The program checks independence by comparing P(A and B) with P(A) times P(B):

{lang="lisp",linenos=off}
~~~~~~~~
(defun independent-p (a b omega)
  "Are A and B independent?  Check whether P(A intersection B) = P(A) * P(B)."
  (let ((intersection (intersection a b :test #'equal)))
    (= (probability-classical intersection omega)
       (* (probability-classical a omega)
          (probability-classical b omega)))))
~~~~~~~~

### Mutual Independence for Three or More Events

For three or more events the definition of independence becomes more subtle. Events A, B, and C are **mutually independent** if the joint probability factors for every subset:

    P(A and B) = P(A) * P(B)
    P(A and C) = P(A) * P(C)
    P(B and C) = P(B) * P(C)
    P(A and B and C) = P(A) * P(B) * P(C)

**Pairwise** independence (the first three conditions) does not imply mutual independence. There are classical examples of three events that are pairwise independent but for which the joint probability P(A and B and C) is not equal to the product P(A) * P(B) * P(C). Whenever a theorem or model calls for "independent" events, it almost always means mutually independent.

## The Law of Total Probability

The **law of total probability** is a tool for decomposing a hard probability into simpler pieces. If B1, B2, ..., Bn partition the sample space (they are disjoint and cover everything), then for any event A:

    P(A) = sum of P(A | Bi) * P(Bi) over all i

Each term P(A | Bi) * P(Bi) answers the question: what is the probability of A and Bi happening together? Summing over all the Bi accounts for all the ways A can happen.

The law of total probability is often applied when you know the conditional probabilities P(A | Bi) and the probabilities P(Bi) of the pieces of the partition, but you do not directly know P(A). The formula stitches these together into a single number.

We will use this law in the next section to compute the denominator in Bayes' theorem.

## Bayes' Theorem

**Bayes' theorem** is one of the most important results in probability theory. It tells us how to reverse the direction of conditioning. If we know P(A | B), we can compute P(B | A):

    P(B | A) = P(A | B) * P(B) / P(A)

Here P(B) is called the **prior**: what we believed about B before seeing evidence. P(A | B) is the **likelihood**: how likely the evidence A is if B is true. P(B | A) is the **posterior**: what we believe about B after seeing the evidence. The denominator P(A) is computed using the law of total probability.

A more general form of Bayes' theorem partitions the space with several hypotheses B1, B2, ..., Bn:

    P(Bi | A) = P(A | Bi) * P(Bi) / (sum over j of P(A | Bj) * P(Bj))

This is the form used in classification problems, where each Bi represents a possible class label and A represents the observed features.

### Bayes' Theorem in Odds Form

There is an alternative and often more convenient form of Bayes' theorem that works with **odds** instead of probabilities. The odds in favor of an event B are defined as P(B) / P(not B). In this form Bayes' theorem becomes:

    posterior odds = prior odds * likelihood ratio

where the **likelihood ratio** is P(A | B) / P(A | not B). This form is elegant because the normalizing constant P(A) drops out. Doctors, forensic scientists, and information retrieval systems often work directly in log-odds; each new piece of evidence adds a constant amount to a running log-odds score.

### The Medical Screening Example

The classic application of Bayes' theorem is medical testing. Suppose we have a disease that affects 1% of the population. We have a test with 99% sensitivity (it correctly identifies 99% of sick people) and 95% specificity (it correctly identifies 95% of healthy people). If a person tests positive, what is the probability they actually have the disease?

Most people intuitively guess around 99%, but the answer is dramatically lower. Let us work through it:

{lang="lisp",linenos=off}
~~~~~~~~
(defun bayes-medical-test ()
  (let* ((p-disease 1/100)            ; prior P(D)
         (p-healthy 99/100)           ; P(no D) = 1 - P(D)
         (p-pos-given-disease 99/100) ; likelihood P(+ | D)
         (p-pos-given-healthy 5/100)) ; likelihood P(+ | no D)
    ;; Total probability of a positive test:
    ;;   P(+) = P(+ | D) P(D) + P(+ | no D) P(no D)
    (let ((p-positive (+ (* p-pos-given-disease p-disease)
                         (* p-pos-given-healthy p-healthy))))
      ;; Bayes: P(D | +) = P(+ | D) P(D) / P(+)
      (let ((posterior (/ (* p-pos-given-disease p-disease) p-positive)))
        ...))))
~~~~~~~~

The prior probability of having the disease is only 1/100. The test has a 5% false-positive rate, and since 99% of the population is healthy, those false positives add up. The total probability of a positive test is:

    P(+) = (99/100)(1/100) + (5/100)(99/100) = 99/10000 + 495/10000 = 594/10000

The posterior probability is:

    P(D | +) = (99/100)(1/100) / (594/10000) = 99/10000 / 594/10000 = 1/6

So a positive test means only about a 16.7% chance of actually having the disease. This counterintuitive result arises because the disease is rare and false positives from the large healthy population dominate the true positives.

### The Base Rate Fallacy

The medical-testing surprise illustrates the **base rate fallacy**, a widespread cognitive bias. When people evaluate evidence, they tend to focus on the likelihood P(evidence | hypothesis) and neglect the prior P(hypothesis). The base rate (the prior) is easy to overlook, yet it dominates the calculation when the hypothesis is rare.

The same fallacy appears in criminal trials (the **prosecutor's fallacy**, in which the probability of the DNA match given innocence is confused with the probability of innocence given the DNA match), in security screening (a very accurate test can still generate mostly false alarms when hunting for extremely rare threats), and in many other high-stakes settings. Understanding Bayes' theorem is not just mathematically satisfying; it is a defense against widely-made real-world reasoning errors.

## Running the Example

When you load the program, you see:

```
=== Conditional Probability with Two Dice ===

B = 'first die is 1',  |B| = 6
A = 'sum is 7',        |A| = 6
P(A)         = 1/6 = .167
P(A | B)     = 1/6 = .167
Independent? yes

e1 = 'first die is 1', e2 = 'second die is 1'
P(e1) = 1/6, P(e2) = 1/6, P(e1 intersection e2) = 1/36
Independent? yes (rolling separate dice never affects each other)

=== Bayes' Theorem ===

Medical screening example (Bayes' theorem):
  Prior P(D)           = 1/100 = 0.01
  Sensitivity P(+|D)   = 99/100 = 0.99
  False-positive P(+|no D) = 1/20 = 0.05
  Total P(+)           = 297/5000 = .059
  Posterior P(D | +)   = 1/6 = .167
  => A positive test means only 16.7% chance of disease!
```

The first part shows that P(A | B) equals P(A) when A and B are independent, confirming the theory. The second part demonstrates the surprising power of Bayes' theorem to correct our intuitions about evidence and prior beliefs.

## Why This Matters

Conditional probability and Bayes' theorem are the engines that power statistical inference, machine learning classification, spam filtering, medical diagnosis, and many other applications. Whenever you need to update your beliefs in light of new evidence, you are applying Bayes' theorem, whether you realize it or not. In a later chapter we will see how Bayesian inference extends this idea to continuous updating of beliefs as data streams in.

## Problem Set

**Problem 2.1.** For the two-dice sample space, let A be the event "sum is 8" and B be the event "first die is 4." Compute P(A), P(B), P(A | B), and P(B | A). Are A and B independent? Give an intuitive explanation for your answer.

**Problem 2.2.** Two cards are drawn from a standard deck without replacement. Compute the probability that both are aces. Solve the problem two ways: (a) by counting the number of favorable ordered pairs and dividing by the total number of ordered pairs, and (b) by using the multiplication rule P(both aces) = P(first ace) * P(second ace | first ace).

**Problem 2.3 (Independence vs. disjointness).** Give an example of two events on the two-dice sample space that are disjoint but not independent. Then give an example of two events that are independent but not disjoint. Explain why disjoint events with positive probability can never be independent.

**Problem 2.4 (The Monty Hall problem).** You are on a game show. There are three doors. Behind one door is a car; behind the other two are goats. You pick door 1. The host, who knows what is behind each door, opens door 3 to reveal a goat. He then offers you the chance to switch to door 2. Should you switch? Use Bayes' theorem to compute the probability of winning if you switch and the probability of winning if you stay. Then modify the example program to simulate the game 100,000 times and empirically confirm your answer.

**Problem 2.5 (The two-child problem).** A family has two children. You are told that at least one of them is a boy. What is the probability that both are boys? (Assume boys and girls are equally likely and the sexes of the two children are independent.) Now suppose instead you are told that the older child is a boy. What is the probability that both are boys? Compare the two answers and explain why they differ.

**Problem 2.6 (A rarer disease).** Redo the medical-testing calculation with a disease prevalence of 1 in 10,000, still using a test with 99% sensitivity and 95% specificity. What is P(D | +) now? What lesson does this teach about screening for rare conditions?

**Problem 2.7 (Sequential testing).** Continuing from the original example (prevalence 1%, sensitivity 99%, specificity 95%), suppose that after the first positive test, the patient is given a second, independent test with the same characteristics, and it also comes back positive. What is the posterior probability of disease now? Hint: use the posterior from the first test as the prior for the second.

**Problem 2.8 (Three-hypothesis Bayes).** A factory has three machines producing widgets: machine A produces 50% of widgets with a 1% defect rate, machine B produces 30% of widgets with a 2% defect rate, and machine C produces 20% of widgets with a 5% defect rate. A widget is selected at random and found to be defective. What is the probability it came from each machine? Verify that your three posteriors sum to 1.

**Problem 2.9 (Coding exercise).** Extend the example program with a function `bayes-update` that takes a prior probability and a likelihood ratio and returns the posterior probability. Use it to reproduce the medical-screening result. Then apply your function to Problem 2.7 and confirm your answer.
