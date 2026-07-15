# Binomial and Geometric Distributions {#binomial_geometric}

Two of the most important discrete distributions arise from repeating the same simple experiment over and over. If we flip a biased coin repeatedly, we can ask two natural questions: how many successes happen in a fixed number of flips, and how long do we wait for the first success? These questions lead to the **binomial** and **geometric** distributions.

The example program for this chapter is in the file **04_binomial_geometric.lisp**.

## Bernoulli Trials

A **Bernoulli trial** is the simplest random experiment with two outcomes: "success" (with probability p) and "failure" (with probability 1-p). If we let X = 1 on success and X = 0 on failure, then X has the Bernoulli(p) distribution with:

    E[X] = p
    Var(X) = p(1-p)

The Bernoulli trial is the atom from which we build more complex distributions.

## The Binomial Distribution

If we run n **independent** Bernoulli(p) trials and count the number of successes, the count X has the **binomial distribution**, written Binomial(n, p). The probability of getting exactly k successes is:

    P(X = k) = C(n,k) * p^k * (1-p)^(n-k)

where C(n,k) = n! / (k! * (n-k)!) is the **binomial coefficient**, the number of ways to choose which k of the n trials are the successes. The factor p^k * (1-p)^(n-k) is the probability of any one specific sequence with k successes and n-k failures, and C(n,k) counts how many such sequences there are.

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

    E[X] = n * p
    Var(X) = n * p * (1-p)

The mean formula has an intuitive explanation via linearity of expectation. Each trial contributes p to the expected count on average, and there are n independent trials, so the expected total is n * p.

## The Geometric Distribution

The **geometric distribution** answers a different question: if we run Bernoulli(p) trials until the **first success** occurs, how many trials do we need? If Y is the number of trials (including the successful one), then:

    P(Y = k) = (1-p)^(k-1) * p

This formula reads as "(k-1) failures followed by one success." The geometric distribution is a valid PMF because it sums to a geometric series: the sum of (1-p)^(k-1) * p for k from 1 to infinity equals p * 1/(1-(1-p)) = 1.

The mean and variance are:

    E[Y] = 1/p
    Var(Y) = (1-p) / p^2

The mean 1/p makes intuitive sense. If the success probability is 1/5, you wait about 5 trials on average. The smaller the success probability, the longer you expect to wait.

## The Memoryless Property

The geometric distribution has a remarkable property called **memorylessness**. It says that the past does not affect the future. Specifically:

    P(Y > m+n | Y > m) = P(Y > n)

If you have already waited m trials without a success, the probability of waiting at least n more trials is the same as if you had just started. The distribution "forgets" how long you have been waiting.

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

The proof is a one-liner. Since P(Y > k) = (1-p)^k, we have:

    P(Y > m+n | Y > m) = P(Y > m+n) / P(Y > m) = (1-p)^(m+n) / (1-p)^m = (1-p)^n = P(Y > n)

## Running the Example

```
=== Binomial Distribution: n=10, p=0.3 ===
  E[X] = n p = 3,  Var(X) = n p (1-p) = 21/10
  P(X=0) = .028   CDF F(0) = 0.028
  P(X=1) = .121   CDF F(1) = 0.149
  P(X=2) = .233   CDF F(2) = 0.383
  P(X=3) = .267   CDF F(3) = 0.650
  ...
  P(X=10) = 0.0   CDF F(10) = 1.0

=== Geometric Distribution: p=0.2 ===
  E[Y] = 1/p = 5,  Var(Y) = (1-p)/p^2 = 20
  P(Y=1) = 0.2   P(Y>1) = 0.8
  P(Y=2) = 0.16  P(Y>2) = 0.64
  ...

=== Memoryless Property ===
  Memoryless property check (p=1/5, m=3, n=2):
    P(Y > m+n | Y > m) = 16/25 = 0.64
    P(Y > n)          = 16/25 = 0.64
```

The binomial distribution with n=10 and p=0.3 is peaked around k=3 (the mean), with the probability declining on both sides. The CDF rises from near 0 to exactly 1, confirming that the PMF sums to 1.

The geometric distribution with p=0.2 is monotonically decreasing: the most likely outcome is Y=1 (a success on the very first trial) with probability 0.2, and longer waits become progressively less likely.

The memoryless property check confirms that P(Y > 5 | Y > 3) equals P(Y > 2), both being 16/25 = 0.64. The past waiting time of 3 trials has no effect on the future.

## Why This Matters

The binomial distribution appears whenever we count successes in repeated independent trials: opinion polling, quality control, A/B testing, and many other applications. The geometric distribution models waiting times: how many calls until a sale, how many attempts until a success, how many packets until a collision. Together, these two distributions cover a huge range of practical discrete probability problems.
