# Monte Carlo Methods {#monte_carlo}

Some quantities are easy to describe but hard to compute exactly. The **Monte Carlo** approach is to estimate them by repeated random sampling. The theoretical justification is the Law of Large Numbers: if we can express a quantity as an expectation, then the sample average of repeated random draws converges to that expectation.

The example program for this chapter is in the file **08_monte_carlo.lisp**.

## The Core Idea

A Monte Carlo estimator works as follows. We want to estimate some quantity Q. We find a random variable X whose expected value E[X] equals Q. Then we draw many independent samples of X and average them:

    Q_hat = (X1 + X2 + ... + Xn) / n

By the Law of Large Numbers, Q_hat converges to Q as n grows. The variance of the estimator is Var(X) / n, so the standard error decreases as 1/sqrt(n).

## Estimating Pi

Our example estimates the value of pi. Consider a unit square [0,1] x [0,1] (area 1) and the quarter disk of radius 1 centered at the origin (area pi/4). If we throw uniformly random points into the square, the probability that a point lands inside the quarter disk equals the ratio of the areas:

    P(point in disk) = (pi/4) / 1 = pi/4

Let X_i = 1 if the i-th point lands in the disk, 0 otherwise. Then the X_i are i.i.d. Bernoulli(pi/4), and by the LLN:

    (1/n) * sum of Xi  converges to  pi/4

So our estimator is:

    pi_hat = 4 * (fraction of points inside the disk)

{lang="lisp",linenos=off}
~~~~~~~~
(defun random-point-in-unit-square ()
  "Return (x, y) with x,y ~ Uniform(0,1) independently."
  (values (random 1.0d0 *rng-state*)
          (random 1.0d0 *rng-state*)))

(defun in-quarter-disk-p (x y)
  "Is (x,y) inside the quarter disk of radius 1?
   The condition is x^2 + y^2 <= 1."
  (<= (+ (* x x) (* y y)) 1.0d0))

(defun estimate-pi (n)
  "Monte Carlo estimate of pi using n random points.
   Estimator = 4 * (count inside disk) / n."
  (let ((inside 0))
    (dotimes (i n)
      (multiple-value-bind (x y) (random-point-in-unit-square)
        (when (in-quarter-disk-p x y)
          (incf inside))))
    (* 4.0d0 (/ inside n 1.0d0))))
~~~~~~~~

A point is inside the quarter disk if x^2 + y^2 <= 1, which is the Pythagorean distance from the origin. Uniform sampling over the square gives each region a probability equal to its area, which is the key property that makes the ratio of counts equal the ratio of areas.

## Standard Error

How accurate is our estimate? For a Bernoulli(pi/4) variable, the standard error of the sample mean M_n is sqrt(p(1-p)/n). The estimator is 4 * M_n, so its standard error is:

    SE = 4 * sqrt( p_hat * (1 - p_hat) / n )

where p_hat is the sample proportion. The program reports both the estimate and the standard error:

{lang="lisp",linenos=off}
~~~~~~~~
(defun estimate-pi-with-se (n)
  "Estimate pi and also report a standard error for the estimate."
  (let* ((inside 0)
         (p-hat 0.0d0)
         (estimate 4.0d0)
         (se 0.0d0))
    (dotimes (i n)
      (multiple-value-bind (x y) (random-point-in-unit-square)
        (when (in-quarter-disk-p x y)
          (incf inside))))
    (setf p-hat (/ inside n 1.0d0))
    (setf estimate (* 4.0d0 p-hat))
    (setf se (* 4.0d0 (sqrt (/ (* p-hat (- 1.0d0 p-hat)) n 1.0d0))))
    (values estimate se)))
~~~~~~~~

## Running the Example

```
=== Monte Carlo Estimation of pi ===
True pi = 3.141593

       n     estimate     error    std-error
  1000      3.032000  0.109593  0.054175
  10000     3.144000  0.002407  0.016405
  100000    3.137640  0.003953  0.005202
  1000000   3.138252  0.003341  0.001645

Note: the error scales roughly as 1/sqrt(n). Doubling n
shrinks the standard error by ~0.707x; halving the error needs 4x n.
```

With 1,000 points, the estimate is off by about 0.1. With 1,000,000 points, the error drops to about 0.003. The standard error column tracks the actual error well: the error is typically within 1 or 2 standard errors of the true value.

## The Cost of Monte Carlo

The standard error decreases as 1/sqrt(n). This is a slow convergence rate. To halve the error, you need 4 times as many samples. To reduce the error by a factor of 10, you need 100 times as many samples.

This is the fundamental tradeoff of Monte Carlo methods. They are incredibly general: you can estimate almost anything by random sampling. But they are slow to converge. For problems where exact computation is feasible, you should prefer exact methods. Monte Carlo shines when exact computation is intractable, such as high-dimensional integrals or complex probabilistic models.

## Why This Matters

Monte Carlo methods are used in physics (particle transport simulations), finance (option pricing), engineering (reliability analysis), and machine learning (Bayesian inference, reinforcement learning). The basic idea is always the same: express the quantity of interest as an expectation, then estimate it by averaging random samples.

The pi estimation example is pedagogically simple, but the same principle scales to problems of enormous complexity. The only requirement is that you can simulate the random variable whose expectation you want. Once you can do that, the Law of Large Numbers does the rest.
