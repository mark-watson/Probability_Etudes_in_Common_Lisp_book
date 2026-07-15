# Markov Chains {#markov_chains}

Many random processes unfold over time: the weather changes day by day, a stock price moves tick by tick, a customer navigates through a website page by page. A **Markov chain** is the simplest model for such sequential randomness. It captures the idea that the future depends on the present but not on the distant past.

The example program for this chapter is in the file **10_markov_chains.lisp**.

## Andrey Markov and the Origin of Markov Chains

The Russian mathematician Andrey Markov introduced these chains in 1906 as a way to study sequences of random variables that are not independent. Independence had been the dominant assumption in classical probability theory (as in Bernoulli's theorem and the Central Limit Theorem), but many real phenomena are far from independent. Markov wanted to find the weakest form of dependence for which the great theorems of probability could still be proved.

His original motivating example was surprisingly literary: he analyzed the sequence of consonants and vowels in Pushkin's novel-in-verse *Eugene Onegin*. He found that consonants and vowels are not independent: a vowel is more likely to be followed by a consonant than by another vowel. But he could still describe the sequence with a two-state Markov chain and could prove a version of the Law of Large Numbers for it. This was the birth of a new branch of probability theory.

Markov chains have since found applications in essentially every scientific field. Statistical mechanics, queueing theory, population genetics, computer science, natural language processing, and machine learning all use Markov chains as fundamental modeling tools.

## The Markov Property

A Markov chain is a sequence of random variables X0, X1, X2, ... taking values in a **state space**. The defining property is the **Markov property**:

    P(X_{t+1} = j | X_t = i, X_{t-1}, ...) = P(X_{t+1} = j | X_t = i)

The next state depends only on the current state, not on the history of how we got there. This is sometimes called "memorylessness of the future given the present."

The Markov property is a very strong assumption. It says that the current state contains all the information needed to predict the future. In a weather model, this would mean that today's weather (Sunny or Rainy) tells us as much about tomorrow's weather as any longer history. In reality, richer models could take into account humidity, pressure, season, and so on, but the Markov idealization is often a useful approximation.

When the Markov property does not hold in the raw data, a common trick is to **expand the state space**: define the state to include enough recent history that the resulting sequence is Markovian. For example, a language model that depends on the last two words can be seen as a Markov chain whose "state" is a pair of consecutive words.

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

The n-step transition matrix P^n has a direct probabilistic meaning: the (i, j) entry is P(X_n = j | X_0 = i), the probability of being in state j after n steps starting from state i. Powers of the transition matrix are one of the primary computational objects when working with Markov chains.

## Classification of States

Not every state in a Markov chain behaves the same way. A rich vocabulary describes the possible behaviors.

A state i **communicates** with state j if there is a positive-probability path from i to j and from j to i. Communication is an equivalence relation; the state space partitions into **communicating classes**. A chain is **irreducible** if it has just one communicating class, meaning every state can eventually reach every other.

A state is **recurrent** if the chain returns to it with probability 1, and **transient** otherwise. In a finite Markov chain, recurrence is equivalent to the state being reachable from itself in the long run.

A state is **absorbing** if once the chain enters it, it never leaves. In matrix terms, an absorbing state has a 1 on the diagonal of P. Absorbing chains have their own rich theory, useful for modeling processes that eventually terminate.

The **period** of a state is the greatest common divisor of the return times to that state. A state is **aperiodic** if its period is 1. Aperiodicity means the chain does not get locked into a deterministic cycle. A chain is **ergodic** if it is irreducible and aperiodic.

Our weather chain is both irreducible (Sunny and Rainy each reach the other) and aperiodic (the chain can stay in the same state, so no forced cycles). It is therefore ergodic, which as we will see next gives it a unique long-run behavior.

## The Stationary Distribution

A **stationary distribution** pi is a row vector that satisfies:

    pi = pi * P

If the chain starts distributed as pi, it stays distributed as pi forever. The stationary distribution is a fixed point of the chain's evolution.

For an **irreducible** chain (you can get from any state to any other state) that is **aperiodic** (not locked into a deterministic cycle), the stationary distribution exists, is unique, and the chain converges to it from any starting state:

    v_t approaches pi  as  t approaches infinity

This is the **ergodic theorem** for Markov chains. No matter where you start, the chain eventually settles into its equilibrium distribution. It generalizes the Law of Large Numbers: the long-run fraction of time the chain spends in state i is pi_i, regardless of the starting state.

### Detailed Balance and Reversibility

A stronger condition than the stationarity equation pi = pi P is **detailed balance**:

    pi_i * P[i][j] = pi_j * P[j][i] for all i, j

If pi and P satisfy detailed balance, then pi is a stationary distribution (summing over j gives the ordinary stationarity equation). A chain that satisfies detailed balance is called **reversible**: viewed in reverse time, it has the same statistical properties.

Detailed balance is a much easier condition to check than the full stationarity equation. Many designed Markov chains, especially those used in Markov chain Monte Carlo, are constructed to satisfy detailed balance for a target distribution pi.

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

For larger chains, we solve the system (P^T - I) pi^T = 0 subject to the constraint sum pi_i = 1. This is a standard linear algebra problem and can also be posed as finding the left eigenvector of P corresponding to eigenvalue 1.

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

## Applications of Markov Chains

Markov chains model a vast range of real-world processes. A partial list:

**PageRank**. Google's original algorithm for ranking web pages models the web as a huge Markov chain: at each step, a "random surfer" clicks a link on the current page, with occasional teleportation to a random page to guarantee irreducibility. The stationary distribution assigns high weight to pages that many random surfers spend time on. PageRank is a computation of this stationary distribution.

**Queueing systems**. In operations research, the number of customers in a queue often forms a Markov chain. Its stationary distribution tells us the long-run behavior of the queue: average wait times, server utilization, and so on.

**Hidden Markov models**. In speech recognition, bioinformatics, and finance, we often observe some quantity that depends on an underlying Markov chain we cannot directly see. Hidden Markov models let us infer the underlying states from the observations.

**Reinforcement learning**. Modern reinforcement learning is built on **Markov decision processes**, which extend Markov chains with actions and rewards. An agent chooses actions to maximize expected cumulative reward in an environment whose state evolves as a Markov chain conditional on the agent's actions.

**Markov chain Monte Carlo (MCMC)**. To sample from a complicated target distribution pi, we can construct a Markov chain whose stationary distribution is exactly pi and run it for a long time. The Metropolis-Hastings algorithm and Gibbs sampling are the two workhorses of MCMC. This has become the dominant computational tool of Bayesian statistics.

**Statistical mechanics**. The state of a physical system at thermal equilibrium is distributed according to the Boltzmann distribution, which is the stationary distribution of many physical Markov chains. MCMC and physics share deep intellectual roots.

## Beyond Discrete Time

The Markov chains we have discussed evolve at discrete time steps. There is a parallel theory of **continuous-time Markov chains**, in which transitions happen at random times drawn from exponential distributions. Continuous-time Markov chains describe systems like radioactive decay, chemical reactions, and queueing systems.

There is also a theory of **general Markov processes** with continuous state spaces, culminating in stochastic differential equations and Brownian motion. These lie beyond the scope of this book, but they build on the same central intuition: the future depends on the present, not on the past.

## Why This Matters

The stationary distribution is particularly important. It tells you the long-run behavior of the system: what fraction of time the chain spends in each state, regardless of where it started. This makes Markov chains a powerful tool for analyzing systems that reach equilibrium, from physical systems to computer networks to economic models.

The convergence to stationarity is also the foundation of Markov chain Monte Carlo, which has become the computational engine of modern Bayesian statistics and countless applications in machine learning. Understanding Markov chains is the entry point to a huge portion of modern probabilistic modeling.

## Problem Set

**Problem 10.1.** For the weather chain with transition matrix [[0.8, 0.2], [0.4, 0.6]], verify by direct substitution that pi = [2/3, 1/3] satisfies pi = pi P.

**Problem 10.2 (Different starting state).** Rerun the example program starting from [0, 1] (certainly Rainy). Does the chain still converge to the same stationary distribution? By what step is the distribution effectively at equilibrium? Explain in words why the answer does not depend on the starting state.

**Problem 10.3 (A three-state chain).** Consider a chain over states A, B, C with transition matrix:

    A -> [0.5, 0.3, 0.2]
    B -> [0.1, 0.6, 0.3]
    C -> [0.2, 0.2, 0.6]

Extend the example program to handle three states and compute the stationary distribution by long iteration. Confirm your answer by verifying pi = pi P.

**Problem 10.4 (Detailed balance).** Does the two-state weather chain satisfy detailed balance for pi = [2/3, 1/3]? Check the condition pi_S * P[S][R] = pi_R * P[R][S] with numerical values. Every two-state chain with a stationary distribution satisfies detailed balance, so this should hold.

**Problem 10.5 (An absorbing chain).** Consider a chain over states {0, 1, 2, 3} where state 0 and state 3 are absorbing (all transitions from them stay in place) and states 1 and 2 each transition to either neighbor with probability 1/2. Starting in state 1, compute the probability of eventually being absorbed at state 0 versus state 3.

**Problem 10.6 (Periodic chain).** Consider a chain over states {0, 1} where state 0 always transitions to state 1 and state 1 always transitions to state 0. Starting from [1, 0], iterate the chain. What do you observe? Explain why this chain does not converge to a stationary distribution, despite the fact that [0.5, 0.5] satisfies the stationarity equation.

**Problem 10.7 (Expected hitting time).** For the weather chain, starting Sunny, what is the expected number of days until the first Rainy day? Hint: the number of days is a geometric random variable with success probability 0.2 (from the Sunny -> Rainy transition). Confirm this both by formula and by simulation.

**Problem 10.8 (Simulating the chain).** Modify the program to simulate one long trajectory of the weather chain: at each step, draw a uniform random number and use it to determine the next state given the current state. Run for 100,000 steps and compute the empirical fraction of Sunny days. Compare against the stationary probability 2/3.

**Problem 10.9 (Convergence rate).** For the weather chain, plot the total variation distance ||v_t - pi||_TV between the distribution at step t and the stationary distribution, as t grows from 0 to 20. How quickly does the distance shrink? The rate is controlled by the second eigenvalue of P.

**Problem 10.10 (Coding exercise).** Add a function `mixing-time` that takes a Markov chain P, a starting distribution, and a tolerance epsilon, and returns the smallest t such that the total variation distance from stationarity is at most epsilon. Test it on the weather chain with epsilon = 0.01.
