# Introduction {#introduction}

Probability theory is the mathematical framework for reasoning about uncertainty. Whenever we make decisions under incomplete information, build models of noisy data, or try to understand random phenomena, probability is the tool we reach for. This book uses short, runnable Common Lisp programs to make the key ideas concrete.

## Why Study Probability?

Probability theory underlies a remarkable range of fields. Statistical inference, machine learning, physics, finance, biology, and engineering all rely on probabilistic models. When you train a neural network, you are implicitly working with probability distributions. When you build a spam filter, you are applying Bayes' theorem. When you model customer behavior or queueing systems, you are using Markov chains.

A solid grasp of probability theory gives you the foundation to understand all of these applications from first principles rather than treating the math as a black box.

## What This Book Covers

We cover the core concepts of a one-semester probability course, organized into ten example programs:

1. **Basic Probability** - sample spaces, events, and the Kolmogorov axioms
2. **Conditional Probability** - conditioning, independence, and Bayes' theorem
3. **Discrete Random Variables** - PMF, CDF, expectation, and variance
4. **Binomial and Geometric Distributions** - two fundamental discrete distributions
5. **Continuous Distributions** - PDF, CDF, and the uniform, exponential, and normal distributions
6. **Law of Large Numbers** - why sample averages converge to the true mean
7. **Central Limit Theorem** - why sums of random variables tend toward the normal distribution
8. **Monte Carlo Methods** - estimating quantities by random sampling
9. **Bayesian Inference** - updating beliefs as evidence arrives
10. **Markov Chains** - sequences of random states with memoryless transitions

Each chapter presents the theory in plain language, then walks through the corresponding Common Lisp program that demonstrates the concepts in action.

## A Note on the Code

The example programs are written to be readable. They are not optimized for performance. The goal is clarity of the probability concepts, not speed of execution. When we need to simulate a million random trials, the code is straightforward even if it takes a few seconds to run.

The code comments explain the probability theory, not the Common Lisp. If you are new to Common Lisp, you can still follow the mathematical ideas by reading the comments and running the programs to see their output.

## How Probability Connects to the Rest of Mathematics

Probability theory sits at a crossroads. It draws on set theory for its foundations, on calculus for continuous distributions, and on linear algebra for multivariate models. In return, it provides the theoretical basis for statistics, machine learning, information theory, and stochastic processes.

In this book we will see these connections emerge naturally. The Kolmogorov axioms are really statements about sets. The expectation of a continuous random variable is an integral. Markov chains involve matrix multiplication. Each connection becomes clear when we implement it in code.
