# Wrapping Up {#wrapup}

We have covered a lot of ground in this short book. Starting from the three Kolmogorov axioms, we built up to the great theorems of probability and their practical applications. Let us look back at the journey.

## The Arc of the Book

We began with **sample spaces and events**, the set-theoretic foundations of probability. The Kolmogorov axioms gave us three simple rules from which everything else follows.

**Conditional probability** taught us to update our view when new information arrives. Bayes' theorem showed us how to reverse the direction of conditioning, revealing that a positive medical test for a rare disease is far less alarming than most people think.

**Random variables** gave us a numerical language for talking about randomness. We learned about PMFs and CDFs for discrete variables, and PDFs for continuous variables. The expected value and variance became our standard tools for summarizing distributions.

The **binomial and geometric distributions** showed us what happens when we repeat simple experiments. The binomial counts successes in a fixed number of trials, while the geometric measures the wait for the first success. Both emerge naturally from Bernoulli trials, and the geometric distribution has the elegant memoryless property.

**Continuous distributions** introduced the uniform, exponential, and normal distributions. We saw that probability is the area under a density curve, and we used numerical integration to compute means, variances, and probabilities. The 68-95-99.7 rule gave us intuition for the normal distribution's shape.

The **Law of Large Numbers** explained why averages converge to expected values. We watched the sample mean of a die roll settle toward 3.5 as the number of rolls grew from 10 to 1,000,000. This theorem is the justification for using sample averages as estimates of true means.

The **Central Limit Theorem** revealed why the normal distribution appears everywhere. Sums of independent random variables tend toward normality, regardless of the original distribution. We saw a bell-shaped histogram emerge from the highly non-normal Bernoulli distribution, just by averaging 50 draws.

**Monte Carlo methods** showed us how to estimate hard quantities by random sampling. We estimated pi by throwing points at a square and counting how many landed in a quarter disk. The 1/sqrt(n) convergence rate is both the power and the limitation of Monte Carlo.

**Bayesian inference** taught us to learn from data by updating beliefs. Starting from a uniform prior, we watched the posterior distribution concentrate around the true parameter value as evidence accumulated. The Beta-Bernoulli conjugate pair made the updates as simple as counting.

**Markov chains** modeled sequential randomness with the Markov property. We saw a weather model converge to its stationary distribution from any starting state, illustrating the ergodic theorem.

Finally, we built **a probabilistic programming language**. It ties the earlier ideas together: a model is declared with a small macro DSL, and three inference engines fit it from the same unnormalized log posterior. Metropolis-Hastings, Hamiltonian Monte Carlo, and variational inference each return the posterior in a different way. The Metropolis-Hastings sampler is a Markov chain whose stationary distribution is the posterior, so the R-hat and ESS diagnostics measure the same mixing behavior the Markov chain chapter analyzed.

## Where To Go From Here

This book covers the foundations. If you want to go further, here are some directions:

**Statistics** builds directly on probability theory. Statistical inference uses probability distributions to draw conclusions from data. Hypothesis testing, confidence intervals, and regression analysis are all applications of the concepts in this book.

**Machine learning** is deeply probabilistic. Bayesian methods, probabilistic graphical models, and Monte Carlo techniques are central to modern ML. The Bayesian inference and Markov chain chapters are particularly relevant.

**Stochastic processes** extend Markov chains to continuous time and more complex state spaces. Poisson processes, Brownian motion, and martingales are the next steps beyond discrete-time Markov chains.

**Information theory** uses probability to quantify uncertainty and information. Concepts like entropy and mutual information are grounded in the probability theory we have covered.

## The Value of Implementation

Every concept in this book was demonstrated with runnable code. I believe that implementing a mathematical idea in code is one of the best ways to truly understand it. When you write a function to compute a probability, you are forced to make the abstract concrete. When you run a simulation, you see the theory come alive in the output.

I encourage you to modify the example programs. Change the parameters, add new distributions, extend the simulations. The best way to learn probability theory is to play with it.

## Final Thoughts

Probability theory is a beautiful and practical branch of mathematics. It gives us the tools to reason clearly about an uncertain world. I hope this book has given you a solid foundation and the confidence to apply these ideas in your own work.

If you enjoyed this book, please consider reading my other book "Loving Common Lisp" for a broader treatment of Common Lisp programming. You can find more of my work at [markwatson.com](http://markwatson.com).
