# Further Reading {#further_reading}

This book covers the foundations. The references below point to authoritative
texts and papers for readers who want to go deeper. The historical sources
named in the chapters (Kolmogorov, Bernoulli, de Moivre, Laplace, Markov, and
others) are not repeated here; this list gathers modern texts and the specific
papers behind the methods in Chapter 11.

## General Probability and Measure

- William Feller, *An Introduction to Probability Theory and Its Applications*,
  Volumes I and II (Wiley). The classic two-volume treatment.
- Sheldon Ross, *A First Course in Probability* (Pearson). A standard
  undergraduate text.
- Geoffrey Grimmett and David Stirzaker, *Probability and Random Processes*
  (Oxford). Strong on random processes.
- Joseph Blitzstein and Jessica Hwang, *Introduction to Probability* (CRC). Clear
  and example-driven, with free online lectures.
- Rick Durrett, *Probability: Theory and Examples* (Cambridge). A graduate
  measure-theoretic text.
- Patrick Billingsley, *Probability and Measure* (Wiley). The measure-theoretic
  foundations sketched in Chapter 1.
- David Williams, *Probability with Martingales* (Cambridge). A short, elegant
  route through measure theory and martingales.

## Statistics, Bayesian Methods, and Machine Learning

- Larry Wasserman, *All of Statistics* (Springer). A fast, broad bridge from
  probability to statistics and machine learning.
- Andrew Gelman, John Carlin, Hal Stern, David Dunson, Aki Vehtari, and Donald
  Rubin, *Bayesian Data Analysis*, 3rd edition (CRC). The standard reference for
  the Bayesian material in Chapters 9 and 11.
- Christopher Bishop, *Pattern Recognition and Machine Learning* (Springer,
  2006).
- Kevin Murphy, *Probabilistic Machine Learning: An Introduction* (MIT, 2022).
- David MacKay, *Information Theory, Inference, and Learning Algorithms*
  (Cambridge, 2003). Free online, and excellent on the link between probability
  and information.

## Markov Chains and Monte Carlo

- David Levin and Yuval Peres, *Markov Chains and Mixing Times* (AMS). The
  eigenvalue and mixing-time material from Chapter 10.
- James Norris, *Markov Chains* (Cambridge).
- Christian Robert and George Casella, *Monte Carlo Statistical Methods*
  (Springer).
- Steve Brooks, Andrew Gelman, Galin Jones, and Xiao-Li Meng, editors,
  *Handbook of Markov Chain Monte Carlo* (CRC, 2011).

## Papers Behind Chapter 11

- Nicholas Metropolis, Arianna Rosenbluth, Marshall Rosenbluth, Augusta Teller,
  and Edward Teller (1953), "Equation of State Calculations by Fast Computing
  Machines," *Journal of Chemical Physics*.
- W. K. Hastings (1970), "Monte Carlo Sampling Methods Using Markov Chains and
  Their Applications," *Biometrika*.
- Simon Duane, Anthony Kennedy, Brian Pendleton, and Duncan Roweth (1987),
  "Hybrid Monte Carlo," *Physics Letters B*. The origin of Hamiltonian Monte
  Carlo.
- Radford Neal (2011), "MCMC Using Hamiltonian Dynamics," in the *Handbook of
  Markov Chain Monte Carlo*.
- Michael Betancourt (2017), "A Conceptual Introduction to Hamiltonian Monte
  Carlo," arXiv:1701.02434.
- Matthew Hoffman and Andrew Gelman (2014), "The No-U-Turn Sampler: Adaptively
  Setting Path Lengths in Hamiltonian Monte Carlo," *Journal of Machine Learning
  Research*. The step-size adaptation used by the HMC engine.
- Alp Kucukelbir, Dustin Tran, Rajesh Ranganath, Andrew Gelman, and David Blei
  (2017), "Automatic Differentiation Variational Inference," *Journal of Machine
  Learning Research*. The recipe the variational engine follows.
- David Blei, Alp Kucukelbir, and Jon McAuliffe (2017), "Variational Inference:
  A Review for Statisticians," *Journal of the American Statistical Association*.
- Diederik Kingma and Jimmy Ba (2015), "Adam: A Method for Stochastic
  Optimization," *International Conference on Learning Representations*. The
  optimizer used by the variational engine.
- Charles Geyer (1992), "Practical Markov Chain Monte Carlo," *Statistical
  Science*. The initial positive sequence estimator behind the ESS diagnostic.
- Andrew Gelman and Donald Rubin (1992), "Inference from Iterative Simulation
  Using Multiple Sequences," *Statistical Science*. The original R-hat statistic.
- Aki Vehtari, Andrew Gelman, Daniel Simpson, Bob Carpenter, and Paul-Christian
  Bürkner (2021), "Rank-normalization, Folding, and Localization: An Improved
  R-hat for Assessing Convergence of MCMC," *Bayesian Analysis*.
- George Box and Mervin Muller (1958), "A Note on the Generation of Random
  Normal Deviates," *Annals of Mathematical Statistics*. The transform used to
  sample normals.
- Atılım Güneş Baydin, Barak Pearlmutter, Alexey Radul, and Jeffrey Siskind
  (2018), "Automatic Differentiation in Machine Learning: A Survey," *Journal of
  Machine Learning Research*. Background for the dual-number differentiation.

## Special Functions

- The NIST Digital Library of Mathematical Functions, https://dlmf.nist.gov, is
  the modern successor to Abramowitz and Stegun and gives the error-function and
  Gamma-function approximations used in Chapters 5, 9, and 11.
