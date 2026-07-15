# Conditional Probability {#conditional_probability}

In the previous chapter we learned to compute probabilities by counting outcomes. But what happens when we gain partial information about an experiment? If we know that one event has occurred, does that change the probability of another event? The answer is yes, and the tool for reasoning about this is **conditional probability**.

The example program for this chapter is in the file **02_conditional_probability.lisp**.

## The Definition of Conditional Probability

The **conditional probability** of event A given event B, written P(A | B), measures how likely A is once we know that B has occurred. The formula is:

    P(A | B) = P(A and B) / P(B)

The intuition is that once we know B happened, B becomes our new sample space. We ask: what fraction of B also contains A? The numerator is the probability that both A and B occur, and the denominator is the probability that B occurs at all.

In the example program, we compute conditional probability by taking the intersection of the two events and dividing by the size of the conditioning event:

{lang="lisp",linenos=off}
~~~~~~~~
(defun conditional-probability (a b omega)
  "P(A | B) = P(A intersection B) / P(B).
   The intersection A intersection B is the set of outcomes in BOTH events."
  (let ((intersection (intersection a b :test #'equal)))
    (/ (length intersection) (length b))))
~~~~~~~~

## Independence

Two events A and B are **independent** if knowing one tells you nothing about the other. Formally, A and B are independent if and only if:

    P(A and B) = P(A) * P(B)

Equivalently, P(A | B) = P(A) when A and B are independent. The condition does not change the probability.

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

## The Law of Total Probability

The **law of total probability** is a tool for decomposing a hard probability into simpler pieces. If B1, B2, ..., Bn partition the sample space (they are disjoint and cover everything), then for any event A:

    P(A) = sum of P(A | Bi) * P(Bi) over all i

Each term P(A | Bi) * P(Bi) answers the question: what is the probability of A and Bi happening together? Summing over all the Bi accounts for all the ways A can happen.

We will use this law in the next section to compute the denominator in Bayes' theorem.

## Bayes' Theorem

**Bayes' theorem** is one of the most important results in probability theory. It tells us how to reverse the direction of conditioning. If we know P(A | B), we can compute P(B | A):

    P(B | A) = P(A | B) * P(B) / P(A)

Here P(B) is called the **prior**: what we believed about B before seeing evidence. P(A | B) is the **likelihood**: how likely the evidence A is if B is true. P(B | A) is the **posterior**: what we believe about B after seeing the evidence. The denominator P(A) is computed using the law of total probability.

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

So a positive test means only about a 16.7% chance of actually having the disease! This counterintuitive result arises because the disease is rare and false positives from the large healthy population dominate the true positives.

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
