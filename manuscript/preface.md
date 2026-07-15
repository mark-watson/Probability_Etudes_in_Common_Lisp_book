# Preface

Probability is the mathematical language we use to reason about uncertainty, and it is one of the most practical branches of mathematics you can learn. This short book uses Common Lisp as the vehicle for exploring probability theory. Each chapter is built around a self-contained program that you can load and run immediately. The code comments and the text focus on the probability concepts, not on the details of Common Lisp itself. If you want to learn Common Lisp programming in depth, I recommend my other book "Loving Common Lisp."

## How To Read This Book

The chapters are ordered to build concepts progressively. We start with sample spaces and the basic axioms of probability, then move through conditional probability, random variables, named distributions, the great convergence theorems, Monte Carlo methods, Bayesian inference, and finally Markov chains. Each chapter can also be read on its own if you already have some background.

I encourage you to run the example programs and have the source code open in an editor as you read. Seeing the output of a simulation and having the code in front of you makes the theory tangible in a way that formulas or just reading the text alone cannot.

## Running the Examples

The example programs are in the parent directory of this manuscript. Each file is self-contained and can be run using LispWorks (which I have aliased to `lw`) or SBCL:

```
lw -eval "(progn (load \"01_basic_probability.lisp\") (quit))"
```

Replace the filename with any of the example files. The examples also run under SBCL and other Common Lisp implementations with the only adjustments being the different command line options for **lw** and **sbcl**.

## Acknowledgments

I would like to thank my wife Carol Watson for her work in editing my books.

Please visit the [author's website](http://markwatson.com).
