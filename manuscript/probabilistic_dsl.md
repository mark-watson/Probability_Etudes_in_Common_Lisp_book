# An Embedded Probabilistic Programming Language {#probabilistic_dsl}

A probabilistic program describes a generative model in ordinary code: it draws latent parameters from priors, then ties those parameters to observed data through likelihoods. A **probabilistic programming language** (PPL) separates the model from the inference. You write the model once; the language supplies a generic engine that returns the posterior distribution of the parameters given the data. The same model can be fit by several different engines without rewriting a line.

This chapter builds a small but complete PPL in Common Lisp. It provides a macro DSL for declaring models and three inference engines written from scratch: random-walk Metropolis-Hastings, Hamiltonian Monte Carlo with automatic differentiation, and mean-field variational inference. A terminal REPL and text plots let you inspect the results.

The example program for this chapter is in the file **11_probabilistic_dsl.lisp**.

## What Inference Computes

Chapter 9 introduced Bayesian updating for a single conjugate pair, where the posterior fell out as a closed-form Beta. Real models rarely have that luxury. A probabilistic program encodes the joint distribution `p(\theta, D)`$ over latent parameters `\theta`$ and data `D`$. The program draws each latent from a **prior** with `sample` and ties latents to data through a **likelihood** with `observe`. Inference targets the posterior

```$
p(\theta \mid D) = \frac{p(D \mid \theta)\, p(\theta)}{p(D)},
```

where the marginal likelihood `p(D) = \int p(D \mid \theta)\, p(\theta)\, d\theta`$ is usually an intractable integral over the parameter space.

The key observation that makes all three engines in this chapter work is that we never need `p(D)`$. Every acceptance ratio and every gradient involves a ratio or a difference of log densities, and the constant `\log p(D)`$ cancels. So each engine needs only the **unnormalized log posterior**

```$
\log g(\theta) = \log p(D \mid \theta) + \log p(\theta).
```

The program computes `\log g(\theta)`$ for any `\theta`$. The engine's job is to explore the shape of `g`$ and report where its mass sits. Sampling engines return a set of draws whose empirical distribution approximates the posterior; the variational engine returns a simple approximating distribution tuned to match it.

## Unconstrained Space and the Jacobian

A prior often lives on a bounded set. A standard deviation is positive; a probability lies in `(0, 1)`$. The samplers in this chapter move on the whole real line, so each latent is reparameterized. We work in an **unconstrained** space `u \in \mathbb{R}`$ and map back to the natural scale with a transform `x = t(u)`$: `x = \exp(u)`$ for positive support, `x = \mathrm{sigmoid}(u)`$ for the unit interval, and a scaled sigmoid for a bounded interval.

A change of variables bends the density. For the transformed variable the density picks up the Jacobian:

```$
\log p_u(u) = \log p_x(t(u)) + \log \left| \frac{dt}{du} \right|.
```

For `x = \exp(u)`$ the Jacobian term is `\log|dt/du| = \log(\exp(u)) = u`$. For the sigmoid it is `\log x + \log(1 - x)`$. The DSL applies these corrections automatically, so the model author never writes a transform or a Jacobian by hand. This keeps the model readable and keeps the math correct.

## Three Ways to Find the Posterior

There are three standard computational routes from `\log g`$ to the posterior, and this chapter implements all three.

**Metropolis-Hastings (MCMC).** Start at a point `u`$. Propose a new point `u' = u + \text{normal noise}`$. Accept it with probability

```$
\alpha = \min\!\left(1,\, \frac{g(u')}{g(u)}\right) = \min\!\left(1,\, \exp(\log g(u') - \log g(u))\right).
```

If accepted, the chain moves to `u'`$; otherwise it stays. Under mild conditions the chain's stationary distribution is the posterior. The method is simple and needs no gradients, but it mixes slowly in high dimensions because random proposals wander.

**Hamiltonian Monte Carlo.** Treat `-\log g(u)`$ as a potential energy and add a momentum variable `p`$ with kinetic energy `\tfrac{1}{2}\|p\|^2`$. The total energy is

```$
H(u, p) = -\log g(u) + \tfrac{1}{2}\|p\|^2.
```

Simulate the Hamiltonian dynamics with the leapfrog integrator for several steps, which slides the state along the posterior's level sets. Accept the endpoint with probability `\min(1, \exp(H(u,p) - H(u',p')))`$. Long, low-rejection moves need the **gradient** of `\log g`$, which we get exactly from automatic differentiation rather than finite differences. HMC mixes far better than random-walk Metropolis in moderate dimensions.

**Variational inference.** Instead of simulating the posterior, fit a simple distribution `q(u)`$ to it. We use a factorized Normal `q(u) = \prod_i \mathcal{N}(m_i, s_i^2)`$. Fitting `q`$ means maximizing the Evidence Lower Bound

```$
\mathrm{ELBO} = \mathbb{E}_{q}[\log g(u)] + \mathcal{H}(q),
```

where `\mathcal{H}(q)`$ is the entropy of `q`$. The **reparameterization trick** writes a sample as `u = m + s \odot \epsilon`$ with `\epsilon \sim \mathcal{N}(0, I)`$, so the expectation becomes a differentiable function of `m`$ and `s`$. The same automatic differentiation supplies the gradient, and Adam ascends the ELBO. Variational inference is fast and deterministic, but it approximates the posterior with a shape (independent Gaussians) that may not match the truth.

The next sections show how the file builds the pieces these engines share: automatic differentiation, a distribution library, support transforms, and the model DSL.

## Forward-Mode Automatic Differentiation

HMC and variational inference both need the gradient of `\log g`$. The file implements **forward-mode automatic differentiation** with dual numbers. A dual number carries a value together with the full gradient vector of that value with respect to the model parameters. Every primitive operation propagates the gradient by the chain rule, so a single evaluation of the log density yields both its value and its exact gradient.

```$
\text{dual}(v, g), \qquad v \in \mathbb{R},\; g \in \mathbb{R}^d.
```

For a unary function `f`$ with derivative `f'`$, the rule is `f(\text{dual}(v, g)) = \text{dual}(f(v),\, f'(v)\, g)`$. For a binary function the partials combine by the chain rule. Constants stay plain doubles, so Metropolis-Hastings, which needs no gradient, pays no AD overhead.

{lang="lisp",linenos=off}
~~~~~~~~
(defstruct (dual (:constructor make-dual (v g)))
  (v 0.0d0 :type double-float)   ; the value
  (g))                            ; simple-vector of partials, length *ad-dim*

(defun d-unary (x val dval)
  "Build f(x): VAL is f(value), DVAL is f'(value). Returns a dual iff X is."
  (if (dual-p x)
      (let* ((n *ad-dim*) (gx (dual-g x))
             (g (make-array n :element-type 'double-float)))
        (dotimes (i n) (setf (aref g i) (* dval (aref gx i))))
        (make-dual val g))
      val))

(defun d-binary (a b val da db)
  "Build f(a,b): VAL is the value, DA and DB the partials wrt a and b."
  (if (or (dual-p a) (dual-p b))
      (let* ((n *ad-dim*) (ga (dg a)) (gb (dg b))
             (g (make-array n :element-type 'double-float)))
        (dotimes (i n)
          (setf (aref g i) (+ (* da (aref ga i)) (* db (aref gb i)))))
        (make-dual val g))
      val))
~~~~~~~~

On top of these two combinators, the file defines differentiable versions of the arithmetic operators. Each one passes the right partials through `d-binary` or `d-unary`:

{lang="lisp",linenos=off}
~~~~~~~~
(defun g+ (a b) (d-binary a b (+ (dv a) (dv b)) 1.0d0 1.0d0))
(defun g* (a b) (let ((av (dv a)) (bv (dv b)))
                  (d-binary a b (* av bv) bv av)))
(defun gexp (x)
  ;; Cap the exponent to keep exp from overflowing to infinity.
  (let* ((ax (min 700.0d0 (dv x))) (e (exp ax))) (d-unary x e e)))
(defun glog (x)
  ;; Clamp the argument into the positive domain of log.
  (let ((ax (max 1d-300 (dv x)))) (d-unary x (log ax) (/ 1.0d0 ax))))
~~~~~~~~

The entry point seeds coordinate `i`$ with the unit dual `e_i`$ (value `0`$ everywhere except partial `i`$ set to `1`$), runs the function once, and reads the value and full gradient off the result. One pass gives the whole gradient vector:

{lang="lisp",linenos=off}
~~~~~~~~
(defun ad-gradient (f x)
  "Differentiate scalar F at point X (a vector of doubles).
   Seeds coordinate i with the unit dual e_i, runs F once, and reads the
   value and full gradient off the result. Returns (values value gradient)."
  (let* ((n (length x))
         (*ad-dim* n)
         (duals (make-array n)))
    (dotimes (i n)
      (let ((g (make-array n :element-type 'double-float :initial-element 0.0d0)))
        (setf (aref g i) 1.0d0)
        (setf (aref duals i) (make-dual (coerce (aref x i) 'double-float) g))))
    (let ((r (funcall f duals)))
      (values (dv r) (dg r)))))
~~~~~~~~

The cost is one evaluation of `f`$ per parameter dimension, which is fine for the small models in this chapter. (Reverse-mode AD would be cheaper for functions with many inputs and one output, but forward mode is simple to build and exact.) The numerical guards in `gexp`$ and `glog`$ matter: a divergent HMC trajectory can push a parameter to an extreme value, and clamping keeps the gradient finite so the trajectory can be rejected rather than crashing the sampler.

## The Distribution Library

Each distribution is a small record holding three things: a log-density function written with the differentiable operators, a plain sampler used to initialize chains, and the **support** of the distribution (`:real`, `:positive`, `:unit`, `:bounded`, or `:discrete`). The support fixes the transform to unconstrained space.

{lang="lisp",linenos=off}
~~~~~~~~
(defstruct (dist (:conc-name dist-))
  name logpdf sampler (support :real) (lo nil) (hi nil))

(defun normal (mu sigma)
  "Normal(mu, sigma) on the whole real line."
  (make-dist
   :name 'normal :support :real
   :sampler (lambda () (rnorm (dv mu) (dv sigma)))
   :logpdf (lambda (x)
             (let ((z (g/ (g- x mu) sigma)))
               (g- (g- (gneg (g* 0.5d0 (gsquare z))) (glog sigma))
                   (* 0.5d0 +log2pi+))))))
~~~~~~~~

The log-density is built from `g-`, `g*`, `g/`, `gsquare`, and `glog`, so it carries gradients when its inputs are duals and computes plainly when they are doubles. The parameters `mu`$ and `sigma`$ may themselves be duals, which is what makes hierarchical models work: a likelihood mean can depend on other latents, and the gradient still flows through.

The library supplies Normal, Half-Normal, Exponential, Gamma, Beta, Uniform, Bernoulli, and Poisson. The Beta and Gamma log densities use a Lanczos approximation to `\log \Gamma`$, and the Poisson likelihood uses `\log k! = \log \Gamma(k+1)`$. These are the same special functions that appeared in Chapter 5 (continuous distributions) and Chapter 9 (the Beta-Bernoulli update), restated here so the file stands alone.

## Support Transforms

Two small functions handle the mapping between constrained and unconstrained space. `constrain-support` maps an unconstrained real `u`$ back to the distribution's support, and `log-jac-support`$ returns the log Jacobian of that map. Both are written with the differentiable operators so the Jacobian term enters the gradient correctly.

{lang="lisp",linenos=off}
~~~~~~~~
(defun constrain-support (support lo hi u)
  "Map an unconstrained real U to the distribution's support."
  (case support
    (:real u)
    (:positive (gexp u))
    (:unit (sigmoidg u))
    (:bounded (g+ lo (g* (- hi lo) (sigmoidg u))))
    (t u)))

(defun log-jac-support (support lo hi u)
  "Log absolute Jacobian log|dx/du| of the support transform at U."
  (case support
    (:real 0.0d0)
    (:positive u)                               ; d/du exp(u) = exp(u); log = u
    (:unit (let ((s (sigmoidg u)))
             (g+ (glog s) (glog (g- 1.0d0 s)))))
    (:bounded (let ((s (sigmoidg u)))
                (g+ (log (- hi lo)) (g+ (glog s) (glog (g- 1.0d0 s))))))
    (t 0.0d0)))
~~~~~~~~

For `:positive`$ support the Jacobian term is just `u`$, because `d/du\, \exp(u) = \exp(u)`$ and `\log(\exp(u)) = u`$. For `:unit`$ support it is `\log \sigma(u) + \log(1 - \sigma(u))`$, the log of the sigmoid derivative. The sampler side uses a plain `log`$ and `logit`$ to invert the transform when placing initial values.

## The Model DSL: defmodel, sample, observe

The macro `defmodel` defines a function that, when called with data, returns a **model** object. The model body runs under a dynamic context in one of two modes:

- `:init` discovers the latents in order and picks dispersed starting values.
- `:density` reads a parameter vector and accumulates the unconstrained log posterior (prior plus Jacobian plus likelihood).

The body never mentions transforms or gradients. `sample` and `observe` handle all of it.

{lang="lisp",linenos=off}
~~~~~~~~
(defmacro defmodel (name (&rest params) &body body)
  "Define a probabilistic model NAME with formal data PARAMS. Inside BODY use
   (sample name dist) to declare a latent and (observe dist datum) to score
   data. Calling (NAME data...) returns a model object to hand to `infer`."
  `(defun ,name ,params
     (make-model :name ',name :param-names ',params
                 :thunk (lambda () ,@body))))
~~~~~~~~

`sample` declares a latent. In `:init`$ mode it records the latent's support and returns a starting value. In `:density`$ mode it reads the next parameter from the vector, transforms it to the natural scale, and adds the prior log density plus the Jacobian to the running total:

{lang="lisp",linenos=off}
~~~~~~~~
(defun sample (name dist)
  "Declare a latent variable NAME with prior DIST; return its current value."
  (let ((ctx *ctx*))
    (ecase (ctx-mode ctx)
      (:init
       ;; Start each chain dispersed but sane: uniform in [-2,2] on the
       ;; unconstrained line. This over-disperses relative to typical
       ;; posteriors (good for R-hat) without the extreme values a vague
       ;; prior draw could produce.
       (let* ((u0 (- (* 4.0d0 (runif)) 2.0d0))
              (x0 (constrain-support (dist-support dist)
                                     (dist-lo dist) (dist-hi dist) u0)))
         (push (make-pstate :name name :support (dist-support dist)
                            :lo (dist-lo dist) :hi (dist-hi dist))
               (ctx-supports ctx))
         (push u0 (ctx-inits ctx))
         x0))
      (:density
       (let* ((i (ctx-index ctx))
              (u (aref (ctx-params ctx) i))
              (x (constrain-support (dist-support dist)
                                    (dist-lo dist) (dist-hi dist) u)))
         (setf (ctx-index ctx) (1+ i))
         (setf (ctx-logdens ctx)
               (g+ (ctx-logdens ctx)
                   (g+ (funcall (dist-logpdf dist) x)
                       (log-jac-support (dist-support dist)
                                        (dist-lo dist) (dist-hi dist) u))))
         x)))))

(defun observe (dist value)
  "Score an observed VALUE under likelihood DIST. Ignored while initializing."
  (let ((ctx *ctx*))
    (when (and ctx (eq (ctx-mode ctx) :density))
      (setf (ctx-logdens ctx)
            (g+ (ctx-logdens ctx)
                (funcall (dist-logpdf dist) (coerce value 'double-float)))))
    value))
~~~~~~~~

The same body serves both modes because `sample` and `observe` dispatch on `ctx-mode`. Running the body once in `:init`$ mode builds the latent list and the start vector; running it in `:density`$ mode with a parameter vector builds the log posterior. Two thin wrappers package these:

{lang="lisp",linenos=off}
~~~~~~~~
(defun model-init (model)
  "Run MODEL once in init mode. Returns (values init-u-vector pstates)."
  (let ((*ctx* (make-ctx :mode :init)))
    (funcall (model-thunk model))
    (values (coerce (nreverse (ctx-inits *ctx*)) 'vector)
            (nreverse (ctx-supports *ctx*)))))

(defun model-logdensity (model)
  "Return a closure u-vector -> unconstrained log posterior (dual-aware)."
  (lambda (u)
    (let ((*ctx* (make-ctx :mode :density :params u)))
      (funcall (model-thunk model))
      (ctx-logdens *ctx*))))
~~~~~~~~

The dynamic variable `*ctx*`$ carries the mode and the accumulator, so the model body reads like a direct description of the generative process. This is the whole point of an embedded DSL: the model is just a Lisp function, and the host language's machinery does the rest.

## The Example Data

The file ships three data generators, one per built-in model. Before fitting, it helps to see what the data look like.

{lang="lisp",linenos=off}
~~~~~~~~
(defun example-coin-data ()
  "40 flips with 28 heads (true bias near 0.7)."
  (append (make-list 28 :initial-element 1) (make-list 12 :initial-element 0)))

(defun example-normal-data (&optional (n 60))
  "N draws from Normal(5, 2)."
  (let ((v (make-dvec n)))
    (dotimes (i n v) (setf (aref v i) (rnorm 5.0d0 2.0d0)))))

(defun example-regression-data (&optional (n 40))
  "N points from y = 1 + 2x + Normal(0,1). Returns (values xs ys)."
  (let ((xs (make-dvec n)) (ys (make-dvec n)))
    (dotimes (i n (values xs ys))
      (let ((x (- (* 4.0d0 (runif)) 2.0d0)))
        (setf (aref xs i) x
              (aref ys i) (+ 1.0d0 (* 2.0d0 x) (rnorm 0.0d0 1.0d0)))))))
~~~~~~~~

The three formats are:

- **Coin**: a list of `0`$ and `1`$, here `28`$ ones and `12`$ zeros: `(1 1 1 0 1 ... 1 0 0)`.
- **Normal**: a vector of doubles drawn from `\mathcal{N}(5, 2)`$, for example `#(4.83 6.21 3.40 5.57 ...)`.
- **Regression**: two parallel vectors `xs`$ and `ys`$, where each `y_i = 1 + 2 x_i + \varepsilon_i`$ with `\varepsilon_i \sim \mathcal{N}(0, 1)`$.

Because the global random state is freshly seeded on every load, the exact data (and so the exact posterior) differ slightly each session. The true parameters stay fixed: coin bias near `0.7`$, normal mean `5`$ and standard deviation `2`$, regression intercept `1`$, slope `2`$, noise `1`$.

## The Three Inference Engines

### Metropolis-Hastings

The simplest engine is a random-walk Metropolis sampler in unconstrained space. At each iteration it adds isotropic Gaussian noise to the current point, evaluates `\log g`$ at the proposal, and accepts with the Metropolis ratio. Divergent proposals (infinite or NaN log density) are rejected. The engine records the acceptance rate and returns posterior draws transformed back to the natural scale.

{lang="lisp",linenos=off}
~~~~~~~~
(defun engine-mh (model &rest args)
  "Random-walk Metropolis in unconstrained space. Gradient-free MCMC."
  (let* ((iters (getf args :iters 6000))
         (burn  (getf args :burn 2000))
         (step  (getf args :step 0.4d0))
         (thin  (getf args :thin 1)))
    (multiple-value-bind (u0 pstates) (model-init model)
      (let* ((f (model-logdensity model))
             (dim (length u0))
             (u (to-dvec u0))
             (lp (dv (funcall f u)))
             (rows nil) (nacc 0) (nprop 0))
        (dotimes (it iters)
          (let ((prop (make-dvec dim)))
            (dotimes (i dim)
              (setf (aref prop i) (+ (aref u i) (rnorm 0.0d0 step))))
            (let ((lp2 (dv (funcall f prop))))
              (incf nprop)
              (when (and (finitep lp2)
                         (< (log (max 1d-300 (runif))) (- lp2 lp)))
                (setf u prop lp lp2) (incf nacc))))
          (when (and (>= it burn) (zerop (mod (- it burn) thin)))
            (push (constrain-draw pstates u) rows)))
        (values (nreverse rows) pstates
                (list :accept-rate (/ nacc (max 1 nprop) 1.0d0)))))))
~~~~~~~~

The proposal scale `step`$ controls the tradeoff between acceptance rate and move size. Too small and nearly every proposal is accepted but the chain crawls; too large and nearly every proposal is rejected and the chain stalls. A rate around `0.2`$ to `0.5`$ usually mixes well for random-walk Metropolis.

### Hamiltonian Monte Carlo

HMC treats `-\log g(u)`$ as a potential energy. Each iteration draws a fresh momentum `p \sim \mathcal{N}(0, I)`$, then simulates the Hamiltonian dynamics with the **leapfrog integrator**. The integrator alternates half-step momentum updates with full-step position updates, using the gradient of the potential at each step. The gradient comes from `ad-gradient`, so it is exact.

{lang="lisp",linenos=off}
~~~~~~~~
(labels ((potential-grad (uv)
           ;; Potential U = -log g; its gradient is -grad(log g).
           (multiple-value-bind (val g) (ad-gradient f uv)
             (let ((ng (make-dvec dim)))
               (dotimes (i dim) (setf (aref ng i) (- (aref g i))))
               (values (- val) ng)))))
  (let ((steps (1+ (random steps0)))      ; jitter path length to avoid resonance
        (p (randn-vec dim))
        (cur-pos (copy-seq u)))
    (multiple-value-bind (cur-pot cur-grad) (potential-grad cur-pos)
      (when (finitep cur-pot)
        (let ((cur-kin 0.0d0))
          (dotimes (i dim) (incf cur-kin (* 0.5d0 (aref p i) (aref p i))))
          (let ((pp (copy-seq p)) (uu (copy-seq cur-pos))
                (grad cur-grad) (prop-pot cur-pot))
            (dotimes (i dim)           ; half step for momentum
              (decf (aref pp i) (* 0.5d0 eps (aref grad i))))
            (block leapfrog
              (dotimes (l steps)       ; full leapfrog steps
                (dotimes (i dim) (incf (aref uu i) (* eps (aref pp i))))
                (multiple-value-bind (vv gg) (potential-grad uu)
                  (setf prop-pot vv grad gg)
                  (unless (finitep vv) (return-from leapfrog)) ; diverged
                  (when (< l (1- steps))
                    (dotimes (i dim)
                      (decf (aref pp i) (* eps (aref grad i))))))))
            (dotimes (i dim)           ; final half step
              (decf (aref pp i) (* 0.5d0 eps (aref grad i))))
            (let ((prop-kin 0.0d0))
              (dotimes (i dim)
                (incf prop-kin (* 0.5d0 (aref pp i) (aref pp i))))
              (incf nprop)
              (let ((alpha (if (and (finitep prop-pot) (finitep prop-kin))
                               (min 1.0d0 (exp (- (+ cur-pot cur-kin)
                                                  (+ prop-pot prop-kin))))
                               0.0d0)))
                (when (< (runif) alpha) (setf u uu) (incf nacc))))))))))
~~~~~~~~

Two details matter in practice. The path length is jittered each iteration (`steps = 1 + random(steps0)`) so the integrator cannot resonate with the target's curvature, which would cripple mixing. And a divergent trajectory (the potential becomes infinite) breaks out of the leapfrog loop early and is rejected, the same safety net as in Metropolis. During warmup a Nesterov dual-averaging scheme adapts the step size `eps`$ toward a target acceptance rate of `0.8`$, then freezes it to the running average for the sampling phase. This is the adaptation recipe from Hoffman and Gelman's NUTS paper, applied to plain HMC.

### Variational Inference

The variational engine fits a factorized Normal `q(u) = \prod_i \mathcal{N}(m_i, s_i^2)`$ by maximizing the ELBO. It stores the mean `m`$ and the log standard deviation `\ell = \log s`$ (so the standard deviation stays positive for free). Each iteration draws `mc`$ samples from `q`$, and for each sample it computes the gradient of `\log g`$ at `u = m + s \odot \epsilon`$. The reparameterization trick moves the gradient of the expectation onto `m`$ and `\ell`$:

```$
\frac{\partial\, \mathrm{ELBO}}{\partial m_i} = \frac{\partial \log g}{\partial u_i}, \qquad
\frac{\partial\, \mathrm{ELBO}}{\partial \ell_i} = \frac{\partial \log g}{\partial u_i}\, s_i\, \epsilon_i + 1.
```

The `+1`$ on the log-standard-deviation gradient is the entropy term `\mathcal{H}(q)`$: it widens `q`$ unless the data pulls the mean toward higher `\log g`$. Adam ascends both parameter sets.

{lang="lisp",linenos=off}
~~~~~~~~
(dotimes (s mc)
  (let ((eps (randn-vec dim)) (u (make-dvec dim)))
    (dotimes (i dim)
      (setf (aref u i) (+ (aref m i)
                          (* (exp (aref ls i)) (aref eps i))))) ; u = m + s*eps
    (multiple-value-bind (val g) (ad-gradient f u)
      (when (and (finitep val) (every #'finitep g))
        (dotimes (i dim)
          ;; dELBO/dm = grad ; dELBO/dls = grad*sd*eps + 1 (entropy).
          (incf (aref gm i) (aref g i))
          (incf (aref gl i)
                (+ (* (aref g i) (exp (aref ls i)) (aref eps i))
                   1.0d0)))))))
~~~~~~~~

After Adam converges, the engine draws a large sample from the fitted `q`$ and returns those draws, transformed back to the natural scale, as its posterior summary. Variational inference is fast and deterministic, but its answer is only as good as the factorized Normal approximation. For a unimodal posterior on the unconstrained scale it does well; for a posterior with strong correlations or multiple modes it can miss.

## Diagnostics: ESS and R-hat

Sampling only approximates the posterior, so the file reports two standard diagnostics alongside each fit.

The **effective sample size** (ESS) measures how many independent draws the autocorrelated chain is worth. If consecutive draws are correlated, the chain carries less information than its length suggests. Geyer's initial positive sequence estimator sums the autocorrelations until they turn negative:

```$
n_{\text{eff}} = \frac{n}{1 + 2 \sum_{k=1}^{\infty} \rho_k},
```

where `\rho_k`$ is the lag-`k`$ autocorrelation. An ESS near `n`$ means the draws are nearly independent; an ESS far below `n`$ means the chain mixes poorly and you should run longer or switch engines.

{lang="lisp",linenos=off}
~~~~~~~~
(defun ess (series)
  "Effective sample size from Geyer's initial positive autocorrelation sum."
  (let* ((n (length series)) (m (mean-of series))
         (var (/ (reduce #'+ (map 'list (lambda (x) (expt (- x m) 2)) series))
                 n)))
    (if (<= var 0.0d0)
        (coerce n 'double-float)
        (let ((s 0.0d0))
          (loop for lag from 1 below n
                for rho = (autocorr series m var lag)
                while (> rho 0.0d0)
                do (incf s rho))
          (clampd (/ n (+ 1.0d0 (* 2.0d0 s))) 1.0d0 (coerce n 'double-float))))))
~~~~~~~~

The **Gelman-Rubin** `\hat{R}`$ (R-hat) compares the variance within each chain to the variance between chains. Run several chains from dispersed starts. If they all sample the same posterior, the between-chain variance should match the within-chain variance and `\hat{R} \approx 1`$. If the chains disagree, `\hat{R}`$ rises above `1`$, signalling that the chains have not converged. The statistic is

```$
\hat{R} = \sqrt{\frac{\hat{V}}{W}}, \qquad \hat{V} = \frac{n-1}{n}\, W + \frac{1}{n}\, B,
```

with `B`$ the between-chain variance and `W`$ the average within-chain variance. A common rule of thumb is to trust the fit when `\hat{R} < 1.01`$ for every parameter. R-hat needs at least two chains; with one chain the file reports ESS only.

## Running the Example

Load the file and it drops you into the PPL REPL:

```
rlwrap sbcl --load 11_probabilistic_dsl.lisp
```

Type `demo` at the `ppl>` prompt for a guided tour, or fit a model directly. The `(demo)` function fits the coin and the Normal models, plots their posteriors, and compares all three engines against the exact Beta answer. Because the random state is freshly seeded each session, your numbers will differ slightly from those below, but the patterns hold.

### The Coin Bias

The coin model puts a `\text{Beta}(1, 1)`$ prior (uniform on `[0, 1]`$) on the bias `p`$ and scores `40`$ flips, `28`$ of them heads, with a Bernoulli likelihood:

{lang="lisp",linenos=off}
~~~~~~~~
(defmodel coin-model (flips)
  "Coin bias: p ~ Beta(1,1), each flip ~ Bernoulli(p)."
  (let ((p (sample :p (beta-dist 1.0d0 1.0d0))))
    (dolist (f flips) (observe (bernoulli p) f))))
~~~~~~~~

With a `\text{Beta}(1, 1)`$ prior and `28`$ heads, `12`$ tails, the exact posterior is `\text{Beta}(29, 13)`$ with mean `29/42 \approx 0.6905`$ (Chapter 9's conjugate update). HMC with two chains reproduces it:

```
Posterior summary  [method hmc, 2000 draws, 2 chain(s)]
  param    mean       sd         2.5%       50%        97.5%      ess      R-hat
  p            0.6859     0.0713     0.5321     0.6890     0.8241    2000.   1.000
  acceptance rate: 0.853

posterior of p
      0.484 | ###
      0.506 | ####
      0.527 | ###
      0.549 | ######
      0.570 | #########
      0.592 | #################
      0.613 | ##########################
      0.635 | ####################################
      0.656 | ####################################
      0.678 | #######################################
      0.699 | ############################################
      0.721 | ########################################
      0.742 | #################################
      0.764 | ############################
      0.785 | #################
      0.807 | ########
      0.828 | ######
      0.850 | ####
      0.871 | #
```

The posterior mean `0.6859`$ sits next to the exact `0.6905`$. The 95% credible interval `[0.532, 0.824]`$ is what the conjugate `\text{Beta}(29, 13)`$ gives. The histogram is the shape of that Beta density: roughly symmetric, peaked near `0.69`$, tapering to both sides. ESS equals the full `2000`$ draws and R-hat is `1.000`$, so the two chains agree and mix almost perfectly. The acceptance rate `0.853`$ is right at the HMC target of `0.8`$.

### Three Engines Versus the Exact Answer

The demo fits the same coin with each engine and checks the posterior mean against the exact `29/42`$:

```
  engine         E[p]       abs err
  Metropolis         0.6896     0.0009
  HMC                0.6884     0.0021
  Variational        0.6746     0.0159
```

All three land within a couple of percentage points of `0.6905`$. Metropolis and HMC, both sampling the posterior directly, are essentially exact up to Monte Carlo noise. Variational inference is less accurate here because the factorized Normal `q`$ is fit in unconstrained space and then pushed through the sigmoid; the mean of a transformed Normal is not the transform of the mean, so the bias shifts slightly. This is the usual cost of variational inference: speed and determinism in exchange for a small, controllable bias.

### Unknown Mean and Spread

The mean-var model puts a `\mathcal{N}(0, 10)`$ prior on the mean, a Half-Normal`(5)`$ prior on the standard deviation, and scores `60`$ draws from `\mathcal{N}(5, 2)`$ with a Normal likelihood:

{lang="lisp",linenos=off}
~~~~~~~~
(defmodel mean-var-model (data)
  "Unknown mean and spread: mu ~ Normal(0,10), sigma ~ HalfNormal(5),
   each datum ~ Normal(mu, sigma)."
  (let ((mu (sample :mu (normal 0.0d0 10.0d0)))
        (sigma (sample :sigma (half-normal 5.0d0))))
    (loop for x across data do (observe (normal mu sigma) x))))
~~~~~~~~

Two chains recover both parameters:

```
Posterior summary  [method hmc, 2400 draws, 2 chain(s)]
  param    mean       sd         2.5%       50%        97.5%      ess      R-hat
  mu           5.0625     0.3098     4.4528     5.0655     5.6644    1768.   1.000
  sigma        2.4944     0.2326     2.0893     2.4749     3.0082    1566.   1.000
  acceptance rate: 0.846
```

The posterior mean of `mu`$ is `5.06`$, bracketing the true `5`$. The posterior mean of `sigma`$ is `2.49`$: with only `60`$ draws the sample standard deviation has upward scatter, and the half-normal prior adds a gentle pull, so the estimate sits a little above the true `2`$. The credible intervals cover both true values. R-hat is `1.000`$ for both and ESS is high, so the chains agree and mix well.

### Bayesian Linear Regression

The regression model puts Normal priors on the intercept and slope and a Half-Normal prior on the noise, then scores `40`$ `(x, y)`$ points. The mean `a + b x`$ is built with the differentiable operators `g+`$ and `g*`$, so HMC and VI get gradients through it:

{lang="lisp",linenos=off}
~~~~~~~~
(defmodel regression-model (xs ys)
  "Linear regression y ~ Normal(a + b x, sigma). The mean a + b*x is built
   with the differentiable ops so HMC and VI get gradients through it."
  (let ((a (sample :a (normal 0.0d0 10.0d0)))
        (b (sample :b (normal 0.0d0 10.0d0)))
        (sigma (sample :sigma (half-normal 5.0d0))))
    (loop for x across xs for y across ys
          do (observe (normal (g+ a (g* b x)) sigma) y))))
~~~~~~~~

Two chains recover the true line `y = 1 + 2x`$ with noise `1`$:

```
Posterior summary  [method hmc, 2400 draws, 2 chain(s)]
  param    mean       sd         2.5%       50%        97.5%      ess      R-hat
  a            1.0078     0.1612     0.6836     1.0036     1.3277    1813.   1.000
  b            1.9109     0.1508     1.6102     1.9089     2.2197    2005.   1.000
  sigma        1.0118     0.1236     0.8070     0.9959     1.2813    1283.   1.001
  acceptance rate: 0.837
```

The intercept `1.008`$ and slope `1.911`$ both cover their true values of `1`$ and `2`$, and the noise `1.012`$ covers the true `1`$. This is the first model with a **linear predictor**: the likelihood mean depends on two latents at once. The AD pass threads the gradient through `g+`$ and `g*`$ without the model author doing anything special, which is what makes the DSL usable for models more complex than a single prior.

### Reading the Diagnostics

The same coin fit by random-walk Metropolis shows why diagnostics matter. A standalone run gives:

```
Posterior summary  [method mh, 10000 draws, 2 chain(s)]
  param    mean       sd         2.5%       50%        97.5%      ess      R-hat
  p            0.6903     0.0694     0.5465     0.6930     0.8153    2234.   1.000
  acceptance rate: 0.533
```

The posterior mean is right on the exact `0.6905`$, but note the cost. Metropolis used `10000`$ iterations to get ESS `2234`$, while HMC got ESS `2000`$ from only `2000`$ draws. Each Metropolis draw is worth roughly a fifth of an independent draw; each HMC draw is worth nearly one. The autocorrelation plot makes the cause visible:

```
autocorrelation of p
  lag   0  1.000 |########################################
  lag   1  0.654 |##########################
  lag   2  0.431 |#################
  lag   3  0.294 |############
  lag   4  0.211 |########
  lag   5  0.149 |######
  lag   6  0.110 |####
  lag   7  0.080 |###
  lag   8  0.060 |##
  lag   9  0.041 |##
  lag  10  0.022 |#
  lag  11  0.017 |#
  lag  12  0.017 |#
  lag  13  0.018 |#
  lag  14  0.006 |
  lag  15 -0.008 |
```

The autocorrelation decays slowly: even at lag `10`$ it is still positive. The area under this curve is the `2 \sum \rho_k`$ in the ESS denominator, and that area is large, so ESS is far below `n`$. HMC's leapfrog moves slide along the posterior and decorrelate quickly, which is why its autocorrelation drops to near zero within a few lags and its ESS sits near `n`$. R-hat is `1.000`$ for both engines, so both have converged; the difference is purely efficiency. This is the practical case for HMC over random-walk Metropolis in even modest dimensions.

## Wrap Up

This chapter assembled a working probabilistic programming language from three ingredients: a tiny DSL (`defmodel`, `sample`, `observe`) that lets you state a model as plain Lisp; a layer of forward-mode automatic differentiation that turns a log-density function into its gradient; and three inference engines that share the same model and the same unnormalized log posterior. Metropolis-Hastings is the simple, gradient-free baseline. HMC adds a gradient and a momentum and mixes far better. Variational inference trades exact sampling for a fast, deterministic approximation.

The whole system fits in one file of roughly `1100`$ lines and uses only Common Lisp primitives plus a few SBCL-specific float-trap masks. It connects the earlier chapters into a single tool: the Normal, Beta, Gamma, and Poisson distributions from Chapters 5 and 9, the Monte Carlo idea from Chapter 8, and the Markov chain theory from Chapter 10 all reappear here. The Metropolis-Hastings chain is a Markov chain whose stationary distribution is the posterior, and the convergence diagnostics R-hat and ESS measure exactly the mixing behavior Chapter 10 analyzed.

The DSL is extensible. Adding a new distribution means writing one `make-dist`$ call with a log density and a sampler. Adding a new model means writing one `defmodel`$ form. To go further you might add a Student-t likelihood for robust regression, a log-normal prior, a No-U-Turn Sampler to tune HMC's path length automatically, or a reparameterization that handles posteriors with strong correlations. The framework is already there; the model and the engine are finally separate.

## Problem Set

**Problem 11.1.** For the coin model with a `\text{Beta}(1, 1)`$ prior and `28`$ heads in `40`$ flips, the exact posterior is `\text{Beta}(29, 13)`$ with mean `29/42 \approx 0.6905`$ and variance `29 \cdot 13 / (42^2 \cdot 43)`$. Compute the exact posterior standard deviation and compare it with the `sd`$ column the HMC summary reports (about `0.0713`$).

**Problem 11.2 (Two modes).** Trace through the coin model body in both `:init`$ and `:density`$ modes. In `:init`$ mode, what value does `sample` return, and what gets pushed onto the context? In `:density`$ mode with parameter vector `u = (u_0)`$, write down the three terms that `sample` adds to `ctx-logdens` (the Beta prior log density, the Jacobian, and the Bernoulli likelihood is added later by `observe`). Why does `observe` do nothing in `:init`$ mode?

**Problem 11.3 (ESS and efficiency).** The Metropolis run used `10000`$ iterations to reach ESS `2234`$, while HMC used `2000`$ draws to reach ESS `2000`$. How many Metropolis iterations would you need to match HMC's `2000`$ effective draws, assuming the same efficiency? What does this say about the per-draw cost of random-walk Metropolis versus HMC?

**Problem 11.4 (R-hat).** Run the coin model with `:method :hmc :chains 4`$ and read the R-hat column. Now deliberately break mixing by setting the HMC step size very large (`:step 2.0d0`) with few leapfrog steps. What happens to the acceptance rate, the ESS, and R-hat? Explain how R-hat detects the problem even when ESS alone might look acceptable.

**Problem 11.5 (The Jacobian).** A standard deviation `sigma`$ has `:positive`$ support, so the DSL uses `x = \exp(u)`$. Show that the log Jacobian is `u`$, as the code claims. Now suppose a probability `p`$ has `:unit`$ support with `x = \mathrm{sigmoid}(u) = 1/(1 + e^{-u})`$. Derive the log Jacobian `\log \sigma(u) + \log(1 - \sigma(u))`$ and confirm it matches `log-jac-support`. What would go wrong if the sampler moved in the natural space and ignored the Jacobian?

**Problem 11.6 (Adding a distribution).** Add a `\text{LogNormal}(\mu, \sigma)`$ distribution to the library. Its density on the positive axis is `p(x) = \frac{1}{x \sigma \sqrt{2\pi}} \exp\!\big(-(\log x - \mu)^2 / (2\sigma^2)\big)`$ for `x > 0`$. Write the `make-dist`$ call using the differentiable operators, choosing the right `:support`. What is the Jacobian contribution, and why is it the same as for any `:positive`$ latent?

**Problem 11.7 (A new model).** Using the existing distributions and the DSL, write a `defmodel`$ for a Poisson rate: `lambda ~ \text{Gamma}(2, 1)`$ prior, and `observe`$ a vector of count data with `(poisson lambda)`. Generate `50`$ counts from a Poisson with true rate `3`$, fit the model with HMC and two chains, and check that the posterior mean of `lambda`$ covers `3`$. Compare with the exact Gamma posterior from Chapter 9's conjugacy (Gamma-Poisson).

**Problem 11.8 (Comparing engines).** Fit the regression model with all three engines (`:mh`, `:hmc`, `:vi`), two chains each where supported. For each engine, report the posterior mean of `b`$, the ESS, and the wall-clock time (use `get-internal-real-time`$ around the call). Rank the engines by accuracy and by speed. Which engine gives the best tradeoff for this three-parameter model?

**Problem 11.9 (Variational bias).** The three-engine table shows variational inference with a larger error than the two samplers. The variational guide is a factorized Normal in unconstrained space, transformed by the sigmoid. Explain in your own words why the mean of `p = \mathrm{sigmoid}(u)`$ with `u \sim \mathcal{N}(m, s^2)`$ is not `sigmoid(m)`$, and why this introduces a bias that grows with the variance `s^2`$. Suggest one way to reduce the bias (for example, fitting `q`$ on the natural scale, or using a Beta guide).

**Problem 11.10 (Coding exercise).** The HMC engine uses a fixed, jittered path length. Replace it with a simple doubling scheme: start with one leapfrog step, and double the path length while the proposed Hamiltonian energy decreases, stopping at a maximum or when the energy starts to rise. This is a simplified version of the No-U-Turn condition. Test your modified engine on the regression model and report the effect on ESS and acceptance rate compared with the fixed-path-length version.
