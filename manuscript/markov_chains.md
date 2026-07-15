# Markov Chains {#markov_chains}

Many random processes unfold over time: the weather changes day by day, a stock price moves tick by tick, a customer navigates through a website page by page. A **Markov chain** is the simplest model for such sequential randomness. It captures the idea that the future depends on the present but not on the distant past.

The example program for this chapter is in the file **10_markov_chains.lisp**.

## The Markov Property

A Markov chain is a sequence of random variables X0, X1, X2, ... taking values in a **state space**. The defining property is the **Markov property**:

    P(X_{t+1} = j | X_t = i, X_{t-1}, ...) = P(X_{t+1} = j | X_t = i)

The next state depends only on the current state, not on the history of how we got there. This is sometimes called "memorylessness of the future given the present."

For a finite state space, the chain is described by a **transition matrix** P, where the entry P[i][j] is the probability of moving from state i to state j in one step. Each row of P is a probability distribution: the entries are non-negative and sum to 1.

## A Weather Model

Our example uses a simple weather model with two states: Sunny and Rainy. The transition matrix is:

- From Sunny: 80% chance of staying Sunny, 20% chance of becoming Rainy
- From Rainy: 40% chance of becoming Sunny, 60% chance of staying Rainy

{lang="lisp",linenos=off}
~~~~~~~~
(let* ((P (make-transition-matrix
            '((0.8 0.2)   ; from Sunny: 80% stay Sunny, 20% -> Rainy
              (0.4 0.6)))) ; from Rainy: 40% -> Sunny, 60% stay Rainy
       (start #(1.0d0 0.0d0)))  ; start certainly Sunny
  ...)
~~~~~~~~

If today is Sunny, there is an 80% chance tomorrow is Sunny too. But if today is Rainy, there is still a 40% chance of a Sunny tomorrow. The weather has persistence: Sunny days tend to follow Sunny days, and Rainy days tend to follow Rainy days.

## Evolving the Distribution

If the state distribution at time t is a row vector v_t, then one step of the chain is a matrix multiplication:

    v_{t+1} = v_t * P

After n steps, the distribution is v_n = v_0 * P^n. The program computes this by repeated multiplication:

{lang="lisp",linenos=off}
~~~~~~~~
(defun step-distribution (dist P)
  "Advance one step: v_{t+1} = v_t P. Returns a new distribution vector."
  (let* ((n (length dist))
         (result (make-array n :initial-element 0.0d0)))
    (dotimes (j n)
      (setf (aref result j)
            (reduce #'+ (loop for i below n
                              collect (* (aref dist i) (aref (aref P i) j))))))
    result))

(defun iterate-chain (dist P steps)
  "Compute v_steps = v_0 P^steps by repeated multiplication."
  (let ((d (copy-seq dist)))
    (dotimes (s steps)
      (setf d (step-distribution d P)))
    d))
~~~~~~~~

## The Stationary Distribution

A **stationary distribution** pi is a row vector that satisfies:

    pi = pi * P

If the chain starts distributed as pi, it stays distributed as pi forever. The stationary distribution is a fixed point of the chain's evolution.

For an **irreducible** chain (you can get from any state to any other state) that is **aperiodic** (not locked into a deterministic cycle), the stationary distribution exists, is unique, and the chain converges to it from any starting state:

    v_t approaches pi  as  t approaches infinity

This is the **ergodic theorem** for Markov chains. No matter where you start, the chain eventually settles into its equilibrium distribution.

## Finding the Stationary Distribution

The program finds the stationary distribution two ways. The first is by **long iteration**: just run the chain for 10,000 steps and read off the distribution:

{lang="lisp",linenos=off}
~~~~~~~~
(defun stationary-by-iteration (dist P steps)
  "Approximate the stationary distribution by long-run simulation of v_t P."
  (iterate-chain dist P steps))
~~~~~~~~

The second is by **solving the linear system** pi = pi * P directly. For a two-state chain with transition matrix [[a, 1-a], [b, 1-b]], the stationary distribution is:

    pi_Sunny = b / (1 - a + b)
    pi_Rainy = (1 - a) / (1 - a + b)

{lang="lisp",linenos=off}
~~~~~~~~
(defun stationary-by-linear-system (P)
  "Solve pi = pi P exactly for a 2-state chain."
  (let* ((a (aref (aref P 0) 0))   ; P[S->S]
         (b (aref (aref P 1) 0))   ; P[R->S]
         (denom (+ (- 1 a) b)))
    (vector (/ b denom) (/ (- 1 a) denom))))
~~~~~~~~

## Running the Example

```
=== Markov Chain: Weather Model ===
States: (Sunny Rainy)
Transition matrix P (rows = from-state):
  Sunny -> [ 0.8,  0.2]
  Rainy -> [ 0.4,  0.6]

Start distribution: [1, 0] (certainly Sunny).
Iterating v_{t+1} = v_t P:
  steps=0   : [1.0000, 0.0000]
  steps=1   : [0.8000, 0.2000]
  steps=2   : [0.7200, 0.2800]
  steps=5   : [0.6701, 0.3299]
  steps=10  : [0.6667, 0.3333]
  steps=50  : [0.6667, 0.3333]
  steps=100 : [0.6667, 0.3333]

Stationary distribution (by long iteration, t=10000):
  [0.6668, 0.3334]
Stationary distribution (exact, solving pi = pi P):
  [0.6667, 0.3333]

The chain forgets its initial state and converges to the
unique stationary distribution (ergodic theorem for Markov chains).
```

Watch the distribution evolve. We start certainly Sunny: [1, 0]. After one step, there is an 80% chance of Sunny and 20% of Rainy. After two steps, the distribution has shifted further: [0.72, 0.28]. By step 10, the distribution has settled at approximately [2/3, 1/3], and it stays there forever.

The stationary distribution is [2/3, 1/3]: in the long run, about 67% of days are Sunny and 33% are Rainy. Both methods of finding it agree. The chain has "forgotten" its starting state and settled into equilibrium.

## Why This Matters

Markov chains model a vast range of real-world processes. Google's original PageRank algorithm is a Markov chain on the web graph. Queueing systems in operations research are Markov chains. Biological sequence models (hidden Markov models) are used in gene prediction and speech recognition. Reinforcement learning in AI is built on Markov decision processes, which extend Markov chains with actions and rewards.

The stationary distribution is particularly important. It tells you the long-run behavior of the system: what fraction of time the chain spends in each state, regardless of where it started. This makes Markov chains a powerful tool for analyzing systems that reach equilibrium, from physical systems to computer networks to economic models.
