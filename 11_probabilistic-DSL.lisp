;;;; 11_probabilistic-DSL.lisp
;;;; An embedded probabilistic programming language (a "PPL") in Common Lisp.
;;;; Define Bayesian models with a small macro DSL, then fit them with three
;;;; inference engines written from scratch: Metropolis-Hastings, Hamiltonian
;;;; Monte Carlo (using automatic differentiation), and mean-field variational
;;;; inference. A terminal REPL and text plots let you explore posteriors.
;;;;
;;;; Load it and it drops you into the PPL REPL:
;;;;   rlwrap sbcl --load 11_probabilistic-DSL.lisp
;;;;   rlwrap lw    -load 11_probabilistic-DSL.lisp
;;;; (rlwrap gives line editing and history.) Type `help` at the ppl> prompt.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; A PROBABILISTIC PROGRAM describes a joint distribution p(theta, D) over
;; latent parameters theta and data D. The program draws each latent from a
;; PRIOR with `sample`, and ties latents to data through a LIKELIHOOD with
;; `observe`. Inference then targets the POSTERIOR
;;
;;     p(theta | D)  =  p(D | theta) p(theta) / p(D),
;;
;; where p(D) = integral of p(D|theta) p(theta) dtheta is usually intractable.
;; Every engine here needs only the UNNORMALIZED log posterior
;;
;;     log g(theta)  =  log p(D | theta) + log p(theta),
;;
;; because the constant p(D) cancels in every accept ratio and every gradient.
;;
;; UNCONSTRAINED SPACE. A prior may live on a bounded set: a standard deviation
;; is positive, a probability lies in (0,1). We reparameterize each latent to
;; the whole real line so the samplers can move freely. For a transform
;; x = t(u) the density picks up the Jacobian:
;;
;;     log p_u(u)  =  log p_x(t(u)) + log |dt/du|.
;;
;; We use x = exp(u) for positive support (log-Jacobian u) and x = sigmoid(u)
;; for the unit interval (log-Jacobian log x + log(1-x)).
;;
;; THREE ENGINES.
;;   1. METROPOLIS-HASTINGS (MCMC). Propose u' = u + normal noise; accept with
;;      probability min(1, g(u')/g(u)). The chain's stationary distribution is
;;      the posterior. Simple, gradient-free, but mixes slowly in high dimension.
;;   2. HAMILTONIAN MONTE CARLO. Treat -log g(u) as a potential energy, add a
;;      momentum variable, and simulate Hamiltonian dynamics with the leapfrog
;;      integrator. Long, low-rejection moves need the GRADIENT of log g, which
;;      we get exactly from forward-mode AUTOMATIC DIFFERENTIATION (dual numbers)
;;      rather than finite differences.
;;   3. VARIATIONAL INFERENCE. Fit a factorized Normal q(u)=prod N(m_i, s_i^2)
;;      by maximizing the Evidence Lower BOund
;;          ELBO = E_q[log g(u)] + H(q).
;;      The reparameterization u = m + s*eps (eps ~ N(0,1)) makes the ELBO
;;      differentiable through the same AD, and Adam ascends it.
;;
;; DIAGNOSTICS. Sampling only approximates the posterior, so we report the
;; effective sample size (ESS), the Gelman-Rubin R-hat across chains, trace
;; plots, autocorrelation, and posterior histograms, all in the terminal.
;;
;; A note on reuse. The earlier chapters (05 continuous distributions, 07 the
;; histogram, 09 the Beta/Bayes update) supply the same math this file builds
;; on. Each of those files runs its own demo on load, so to keep this REPL
;; clean we restate the few helpers we need here and mark where they came from.
;; ============================================================================

;; ----------------------------------------------------------------------------
;; 0. Globals and scalar math helpers
;; ----------------------------------------------------------------------------

(defvar *ppl-rng* (make-random-state t)
  "Random state for the whole PPL, freshly seeded so each session differs.")

(defvar *ppl-autostart* t
  "When true, loading this file starts the REPL. Set to nil before loading to
   just define the functions (for scripting or testing).")

(defparameter +2pi+ (* 2.0d0 (coerce pi 'double-float)))
(defparameter +log2pi+ (log +2pi+))

(defun runif ()
  "Uniform(0,1) draw as a double-float."
  (random 1.0d0 *ppl-rng*))

(defun clampd (x lo hi)
  "Clamp X into [LO, HI]."
  (max lo (min hi x)))

(defun sigmoid (u)
  "Logistic sigmoid 1/(1+e^-u) on a plain double."
  (/ 1.0d0 (+ 1.0d0 (exp (- u)))))

(defun logit (p)
  "Inverse sigmoid log(p/(1-p))."
  (log (/ p (- 1.0d0 p))))

(defun rnorm (&optional (mean 0.0d0) (sd 1.0d0))
  "One Normal(mean, sd^2) draw via the Box-Muller transform.
   Box-Muller maps two uniforms to an independent standard normal pair; we
   keep one. See chapter 05 for the normal density this samples from."
  (let ((u1 (max 1d-12 (runif)))
        (u2 (runif)))
    (+ mean (* sd (sqrt (* -2.0d0 (log u1))) (cos (* +2pi+ u2))))))

;; Lanczos approximation to log Gamma, needed for the Beta and Gamma log
;; densities and for log-factorial in the Poisson likelihood.
(defparameter +lanczos+
  (vector 0.99999999999980993d0 676.5203681218851d0 -1259.1392167224028d0
          771.32342877765313d0 -176.61502916214059d0 12.507343278686905d0
          -0.13857109526572012d0 9.9843695780195716d-6 1.5056327351493116d-7))

(defun log-gamma (x)
  "Natural log of the Gamma function for real x > 0 (Lanczos, g = 7)."
  (let ((x (coerce x 'double-float)))
    (if (< x 0.5d0)
        ;; Reflection: Gamma(x)Gamma(1-x) = pi/sin(pi x).
        (- (log (/ (coerce pi 'double-float) (sin (* (coerce pi 'double-float) x))))
           (log-gamma (- 1.0d0 x)))
        (let* ((x (- x 1.0d0))
               (a (aref +lanczos+ 0))
               (tt (+ x 7.5d0)))
          (loop for i from 1 below 9
                do (incf a (/ (aref +lanczos+ i) (+ x i))))
          (+ (* 0.5d0 (log +2pi+))
             (* (+ x 0.5d0) (log tt))
             (- tt)
             (log a))))))

(defun log-factorial (k)
  "log(k!) = log Gamma(k+1) for a non-negative integer count k."
  (log-gamma (+ (coerce k 'double-float) 1.0d0)))

(defun rgamma (shape)
  "Draw from Gamma(shape, rate=1) with the Marsaglia and Tsang method.
   For shape < 1 we boost with the standard u^(1/shape) correction."
  (let ((shape (coerce shape 'double-float)))
    (if (< shape 1.0d0)
        (* (rgamma (+ shape 1.0d0))
           (expt (max 1d-12 (runif)) (/ 1.0d0 shape)))
        (let* ((d (- shape (/ 1.0d0 3.0d0)))
               (c (/ 1.0d0 (sqrt (* 9.0d0 d)))))
          (loop
            (let* ((xn (rnorm))
                   (v (expt (+ 1.0d0 (* c xn)) 3)))
              (when (> v 0.0d0)
                (let ((u (runif)))
                  (when (< (log u)
                           (+ (* 0.5d0 xn xn) (* d (+ 1.0d0 (- v) (log v)))))
                    (return (* d v)))))))))))

(defun rpois (lam)
  "Draw from Poisson(lam) with Knuth's product method (fine for small lam)."
  (let ((el (exp (- lam))) (k 0) (p 1.0d0))
    (loop
      (incf k)
      (setf p (* p (runif)))
      (when (<= p el) (return (- k 1))))))

;; ----------------------------------------------------------------------------
;; 1. Forward-mode automatic differentiation (dual numbers)
;; ----------------------------------------------------------------------------
;;
;; A DUAL number carries a value plus the whole gradient of that value with
;; respect to the model parameters. Every primitive op propagates the gradient
;; by the chain rule, so ONE evaluation of the log density yields its exact
;; gradient. This is what turns a plain Lisp log-density into the force field
;; HMC and VI need. Constants stay plain doubles; the ops below accept a mix.

(defvar *ad-dim* 1
  "Length of the gradient vectors during a differentiation pass.")

(defstruct (dual (:constructor make-dual (v g)))
  (v 0.0d0 :type double-float)   ; the value
  (g))                            ; simple-vector of partials, length *ad-dim*

(declaim (inline dv))
(defun dv (x)
  "Value part of X, whether X is a dual or a plain real."
  (if (dual-p x) (dual-v x) (coerce x 'double-float)))

(defun zero-grad ()
  (make-array *ad-dim* :element-type 'double-float :initial-element 0.0d0))

(defun dg (x)
  "Gradient part of X: its stored partials, or zeros if X is a constant."
  (if (dual-p x) (dual-g x) (zero-grad)))

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

;; Differentiable arithmetic. Model code that combines latents (a linear
;; predictor, say) uses these so gradients flow. On plain doubles they just
;; compute, so Metropolis-Hastings pays no AD overhead.
(defun g+ (a b) (d-binary a b (+ (dv a) (dv b)) 1.0d0 1.0d0))
(defun g- (a b) (d-binary a b (- (dv a) (dv b)) 1.0d0 -1.0d0))
(defun g* (a b) (let ((av (dv a)) (bv (dv b)))
                  (d-binary a b (* av bv) bv av)))
(defun g/ (a b)
  ;; Guard the denominator so a divergent trajectory cannot divide by zero.
  (let* ((av (dv a)) (b0 (dv b))
         (bv (if (< (abs b0) 1d-300) (if (minusp b0) -1d-300 1d-300) b0)))
    (d-binary a b (/ av bv) (/ 1.0d0 bv) (/ (- av) (* bv bv)))))
(defun gneg (x) (d-unary x (- (dv x)) -1.0d0))
(defun gexp (x)
  ;; Cap the exponent to keep exp from overflowing to infinity.
  (let* ((ax (min 700.0d0 (dv x))) (e (exp ax))) (d-unary x e e)))
(defun glog (x)
  ;; Clamp the argument into the positive domain of log.
  (let ((ax (max 1d-300 (dv x)))) (d-unary x (log ax) (/ 1.0d0 ax))))
(defun gsqrt (x)
  (let* ((ax (max 0.0d0 (dv x))) (s (sqrt ax)))
    (d-unary x s (/ 0.5d0 (max 1d-300 s)))))
(defun gsquare (x) (let ((ax (dv x))) (d-unary x (* ax ax) (* 2.0d0 ax))))
(defun gexpt (x p) (let ((ax (dv x)))
                     (d-unary x (expt ax p) (* p (expt ax (- p 1))))))

(defun sigmoidg (u)
  "Differentiable logistic sigmoid."
  (g/ 1.0d0 (g+ 1.0d0 (gexp (gneg u)))))

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

;; ----------------------------------------------------------------------------
;; 2. Distribution library
;; ----------------------------------------------------------------------------
;;
;; Each distribution knows its log density (in differentiable ops), a sampler
;; (plain draws, used to initialize chains), and its SUPPORT, which fixes the
;; transform to unconstrained space. Parameters may themselves be duals, so a
;; likelihood mean can depend on other latents (hierarchical models).

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

(defun half-normal (sigma)
  "Half-Normal(sigma): a Normal folded to the positive axis. Good scale prior."
  (make-dist
   :name 'half-normal :support :positive
   :sampler (lambda () (abs (rnorm 0.0d0 (dv sigma))))
   :logpdf (lambda (x)
             (g- (g- (* 0.5d0 (log (/ 2.0d0 (coerce pi 'double-float))))
                     (glog sigma))
                 (g* 0.5d0 (gsquare (g/ x sigma)))))))

(defun exponential (rate)
  "Exponential(rate) on the positive axis. See chapter 05."
  (make-dist
   :name 'exponential :support :positive
   :sampler (lambda () (/ (- (log (max 1d-12 (runif)))) (dv rate)))
   :logpdf (lambda (x) (g- (glog rate) (g* rate x)))))

(defun gamma-dist (shape rate)
  "Gamma(shape, rate) on the positive axis."
  (let ((const (- (* shape (log (coerce rate 'double-float)))
                  (log-gamma shape))))
    (make-dist
     :name 'gamma :support :positive
     :sampler (lambda () (/ (rgamma shape) (dv rate)))
     :logpdf (lambda (x)
               (g+ (g- (g* (- shape 1.0d0) (glog x)) (g* rate x)) const)))))

(defun beta-dist (a b)
  "Beta(a, b) on the unit interval. Conjugate coin-bias prior of chapter 09."
  (let ((logB (- (+ (log-gamma a) (log-gamma b)) (log-gamma (+ a b)))))
    (make-dist
     :name 'beta :support :unit
     :sampler (lambda () (let ((g1 (rgamma a)) (g2 (rgamma b)))
                           (/ g1 (+ g1 g2))))
     :logpdf (lambda (x)
               (g- (g+ (g* (- a 1.0d0) (glog x))
                       (g* (- b 1.0d0) (glog (g- 1.0d0 x))))
                   logB)))))

(defun uniform (lo hi)
  "Uniform(lo, hi) on a bounded interval."
  (make-dist
   :name 'uniform :support :bounded :lo (coerce lo 'double-float)
   :hi (coerce hi 'double-float)
   :sampler (lambda () (+ lo (* (- hi lo) (runif))))
   :logpdf (let ((c (- (log (- hi lo)))))
             (lambda (x) (declare (ignore x)) c))))

(defun bernoulli (p)
  "Bernoulli(p): a 0/1 likelihood for coin-flip style data."
  (make-dist
   :name 'bernoulli :support :discrete
   :sampler (lambda () (if (< (runif) (dv p)) 1.0d0 0.0d0))
   :logpdf (lambda (x)
             (g+ (g* x (glog p)) (g* (- 1.0d0 x) (glog (g- 1.0d0 p)))))))

(defun poisson (lam)
  "Poisson(lam): a count likelihood."
  (make-dist
   :name 'poisson :support :discrete
   :sampler (lambda () (coerce (rpois (dv lam)) 'double-float))
   :logpdf (lambda (k)
             (g- (g- (g* k (glog lam)) lam) (log-factorial (round (dv k)))))))

;; ----------------------------------------------------------------------------
;; 3. Support transforms (constrained <-> unconstrained real line)
;; ----------------------------------------------------------------------------

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

(defun unconstrain-support (support lo hi x)
  "Inverse transform, on plain doubles, used to place initial values."
  (case support
    (:real x)
    (:positive (log (max 1d-12 x)))
    (:unit (logit (clampd x 1d-6 (- 1.0d0 1d-6))))
    (:bounded (logit (clampd (/ (- x lo) (- hi lo)) 1d-6 (- 1.0d0 1d-6))))
    (t x)))

;; ----------------------------------------------------------------------------
;; 4. The model DSL: defmodel / sample / observe
;; ----------------------------------------------------------------------------
;;
;; `defmodel` defines a function that, when given data, returns a MODEL object.
;; The model body runs under a dynamic context in one of two modes:
;;   :init    discover the latents in order and pick starting values.
;;   :density read a parameter vector and accumulate the unconstrained log
;;            posterior (prior + Jacobian + likelihood).
;; The body never mentions transforms or gradients; `sample` and `observe`
;; handle all of it.

(defstruct (model (:conc-name model-)) name param-names thunk)
(defstruct (pstate (:conc-name pstate-)) name support lo hi)

(defmethod print-object ((m model) stream)
  (print-unreadable-object (m stream :type t)
    (format stream "~(~a~) ~a" (model-name m) (model-param-names m))))

(defstruct (ctx (:conc-name ctx-))
  mode                 ; :init or :density
  params               ; vector of unconstrained values (density mode)
  (index 0)            ; next parameter slot (density mode)
  (logdens 0.0d0)      ; accumulator (density mode)
  (supports nil)       ; reversed list of pstate (init mode)
  (inits nil))         ; reversed list of unconstrained starts (init mode)

(defvar *ctx* nil "Dynamic model-execution context.")

(defmacro defmodel (name (&rest params) &body body)
  "Define a probabilistic model NAME with formal data PARAMS. Inside BODY use
   (sample name dist) to declare a latent and (observe dist datum) to score
   data. Calling (NAME data...) returns a model object to hand to `infer`."
  `(defun ,name ,params
     (make-model :name ',name :param-names ',params
                 :thunk (lambda () ,@body))))

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
                (funcall (dist-logpdf dist) (coerce value 'double-float))))))
  value)

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

(defun constrain-draw (pstates u)
  "Transform an unconstrained draw U back to the parameters' natural scale."
  (let* ((n (length u))
         (v (make-array n :element-type 'double-float)))
    (loop for i below n for ps in pstates
          do (setf (aref v i)
                   (coerce (constrain-support (pstate-support ps)
                                              (pstate-lo ps) (pstate-hi ps)
                                              (aref u i))
                           'double-float)))
    v))

;; ----------------------------------------------------------------------------
;; 5. Small vector helpers
;; ----------------------------------------------------------------------------

(defun make-dvec (n &optional (init 0.0d0))
  (make-array n :element-type 'double-float :initial-element init))

(defun randn-vec (n)
  (let ((v (make-dvec n)))
    (dotimes (i n v) (setf (aref v i) (rnorm)))))

(defun to-dvec (seq)
  (map '(vector double-float) (lambda (x) (coerce x 'double-float)) seq))

(defun finitep (x)
  "True when X is a finite double (not an infinity or NaN)."
  (and (= x x) (< (abs x) 1d308)))

(defmacro with-safe-floats (&body body)
  "Run BODY with floating traps masked so a divergent sampler trajectory
   yields infinities or NaNs (which we reject) instead of crashing. On
   SBCL this masks the traps; elsewhere it relies on the primitive clamps."
  #+sbcl `(sb-int:with-float-traps-masked (:invalid :overflow :divide-by-zero
                                           :underflow :inexact)
            ,@body)
  #-sbcl `(progn ,@body))

;; ----------------------------------------------------------------------------
;; 6. Inference engine: Metropolis-Hastings
;; ----------------------------------------------------------------------------

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

;; ----------------------------------------------------------------------------
;; 7. Inference engine: Hamiltonian Monte Carlo
;; ----------------------------------------------------------------------------

(defun engine-hmc (model &rest args)
  "HMC with leapfrog integration. Uses AD for the exact gradient of log g."
  (let* ((iters  (getf args :iters 2500))
         (burn   (getf args :burn 1000))
         (eps0   (getf args :step 0.1d0))
         (steps0 (getf args :leapfrog 20))
         (thin   (getf args :thin 1)))
    (multiple-value-bind (u0 pstates) (model-init model)
      (let* ((f (model-logdensity model))
             (dim (length u0))
             (u (to-dvec u0))
             (rows nil) (nacc 0) (nprop 0)
             (eps eps0)
             ;; Nesterov dual-averaging state to adapt the step size in warmup
             ;; toward a target acceptance rate (Hoffman and Gelman, 2014).
             (mu-da (log (* 10.0d0 eps0)))
             (log-eps-bar 0.0d0) (h-bar 0.0d0)
             (gamma 0.05d0) (t0 10.0d0) (kappa 0.75d0) (target 0.8d0))
        (labels ((potential-grad (uv)
                   ;; Potential U = -log g; its gradient is -grad(log g).
                   (multiple-value-bind (val g) (ad-gradient f uv)
                     (let ((ng (make-dvec dim)))
                       (dotimes (i dim) (setf (aref ng i) (- (aref g i))))
                       (values (- val) ng)))))
          (dotimes (it iters)
            ;; Jitter the path length so the leapfrog cannot resonate with the
            ;; target's curvature (which would cripple mixing).
            (let ((steps (1+ (random steps0)))
                  (p (randn-vec dim))
                  (cur-pos (copy-seq u))
                  (alpha 0.0d0))
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
                        (setf alpha (if (and (finitep prop-pot) (finitep prop-kin))
                                        (min 1.0d0 (exp (- (+ cur-pot cur-kin)
                                                           (+ prop-pot prop-kin))))
                                        0.0d0))
                        (when (< (runif) alpha) (setf u uu) (incf nacc)))))))
              ;; Adapt the step size during warmup, then freeze it to the
              ;; running average for the sampling phase.
              (cond ((and (< it burn) (> burn 0))
                     (let* ((m (+ it 1.0d0)) (w (/ 1.0d0 (+ m t0))))
                       (setf h-bar (+ (* (- 1.0d0 w) h-bar) (* w (- target alpha))))
                       (let ((log-eps (- mu-da (* (/ (sqrt m) gamma) h-bar)))
                             (eta (expt m (- kappa))))
                         (setf log-eps-bar
                               (+ (* eta log-eps) (* (- 1.0d0 eta) log-eps-bar)))
                         (setf eps (exp log-eps)))))
                    ((and (= it burn) (> burn 0)) (setf eps (exp log-eps-bar)))))
            (when (and (>= it burn) (zerop (mod (- it burn) thin)))
              (push (constrain-draw pstates u) rows))))
        (values (nreverse rows) pstates
                (list :accept-rate (/ nacc (max 1 nprop) 1.0d0)
                      :step-size eps))))))

;; ----------------------------------------------------------------------------
;; 8. Inference engine: mean-field variational inference (ADVI-style)
;; ----------------------------------------------------------------------------

(defstruct (adam (:conc-name adam-)) (iter 0) m v (lr 0.05d0))

(defun adam-init (dim lr) (make-adam :m (make-dvec dim) :v (make-dvec dim) :lr lr))

(defun adam-ascend! (state params grad)
  "Take one Adam ascent step: PARAMS += lr * m_hat / (sqrt(v_hat)+eps)."
  (let ((tt (1+ (adam-iter state))) (m (adam-m state)) (v (adam-v state))
        (lr (adam-lr state)) (b1 0.9d0) (b2 0.999d0) (eps 1d-8))
    (setf (adam-iter state) tt)
    (dotimes (i (length params))
      (let ((gi (aref grad i)))
        (setf (aref m i) (+ (* b1 (aref m i)) (* (- 1.0d0 b1) gi)))
        (setf (aref v i) (+ (* b2 (aref v i)) (* (- 1.0d0 b2) gi gi)))
        (let ((mhat (/ (aref m i) (- 1.0d0 (expt b1 tt))))
              (vhat (/ (aref v i) (- 1.0d0 (expt b2 tt)))))
          (incf (aref params i) (* lr (/ mhat (+ (sqrt vhat) eps)))))))))

(defun engine-vi (model &rest args)
  "Fit a factorized Normal q(u) by maximizing the ELBO with reparameterized
   gradients (AD) and Adam. Then draw from q to summarize the posterior."
  (let* ((iters  (getf args :iters 3000))
         (lr     (getf args :step 0.05d0))
         (mc     (getf args :mc 4))
         (ndraws (getf args :draws 4000)))
    (multiple-value-bind (u0 pstates) (model-init model)
      (let* ((f (model-logdensity model))
             (dim (length u0))
             (m (to-dvec u0))
             (ls (make-dvec dim -1.0d0))   ; log-sd, so sd starts near 0.37
             (am (adam-init dim lr))
             (als (adam-init dim lr)))
        (dotimes (it iters)
          (let ((gm (make-dvec dim)) (gl (make-dvec dim)))
            (dotimes (s mc)
              (let ((eps (randn-vec dim)) (u (make-dvec dim)))
                (dotimes (i dim)
                  (setf (aref u i) (+ (aref m i)
                                      (* (exp (aref ls i)) (aref eps i)))))
                (multiple-value-bind (val g) (ad-gradient f u)
                  (when (and (finitep val) (every #'finitep g))
                    (dotimes (i dim)
                      ;; dELBO/dm = grad ; dELBO/dls = grad*sd*eps + 1 (entropy).
                      (incf (aref gm i) (aref g i))
                      (incf (aref gl i)
                            (+ (* (aref g i) (exp (aref ls i)) (aref eps i))
                               1.0d0)))))))
            (dotimes (i dim)
              (setf (aref gm i) (/ (aref gm i) mc)
                    (aref gl i) (/ (aref gl i) mc)))
            (adam-ascend! am m gm)
            (adam-ascend! als ls gl)))
        (let ((rows nil))
          (dotimes (k ndraws)
            (let ((u (make-dvec dim)))
              (dotimes (i dim)
                (setf (aref u i) (+ (aref m i) (* (exp (aref ls i)) (rnorm)))))
              (push (constrain-draw pstates u) rows)))
          (values (nreverse rows) pstates
                  (list :vi-mean (copy-seq m)
                        :vi-sd (map '(vector double-float) #'exp ls))))))))

;; ----------------------------------------------------------------------------
;; 9. Posterior object, the `infer` front door, and diagnostics
;; ----------------------------------------------------------------------------

(defstruct (posterior (:conc-name posterior-))
  method param-names supports draws chains extra)

(defmethod print-object ((p posterior) stream)
  (print-unreadable-object (p stream :type t)
    (format stream "~(~a~) ~d draws ~a"
            (posterior-method p) (length (posterior-draws p))
            (coerce (posterior-param-names p) 'list))))

(defun infer (model &rest args)
  "Fit MODEL and return a posterior. Keys:
     :method  :hmc (default) | :mh | :vi
     :chains  number of independent chains (default 1; needed for R-hat)
     plus engine keys such as :iters :burn :step :leapfrog :thin :draws :mc."
  (let* ((method (getf args :method :hmc))
         (chains (getf args :chains 1))
         (engine (ecase method
                   (:mh #'engine-mh) (:hmc #'engine-hmc) (:vi #'engine-vi)))
         (chain-rows nil) (pstates nil) (extra nil))
    (dotimes (c chains)
      (multiple-value-bind (rows ps ex)
          (with-safe-floats (apply engine model args))
        (setf pstates ps extra ex)
        (push (coerce rows 'vector) chain-rows)))
    (setf chain-rows (nreverse chain-rows))
    (let ((all (apply #'concatenate 'list
                      (mapcar (lambda (v) (coerce v 'list)) chain-rows))))
      (make-posterior :method method
                      :param-names (map 'vector #'pstate-name pstates)
                      :supports pstates
                      :draws (coerce all 'vector)
                      :chains chain-rows
                      :extra extra))))

(defun post-index (post name)
  (position name (posterior-param-names post)))

(defun post-column (post name)
  "Extract parameter NAME across all draws as a double vector."
  (let* ((j (post-index post name))
         (rows (posterior-draws post))
         (n (length rows))
         (v (make-dvec n)))
    (unless j (error "no parameter ~s in posterior" name))
    (dotimes (i n v) (setf (aref v i) (aref (aref rows i) j)))))

(defun quantile (sorted q)
  "Linear-interpolation quantile Q of an ascending double vector."
  (let* ((n (length sorted)))
    (if (= n 1) (aref sorted 0)
        (let* ((h (* (- n 1) q)) (lo (floor h)) (hi (min (- n 1) (ceiling h))))
          (+ (aref sorted lo)
             (* (- h lo) (- (aref sorted hi) (aref sorted lo))))))))

(defun mean-of (series)
  (/ (reduce #'+ series) (length series)))

(defun sd-of (series)
  (let ((m (mean-of series)) (n (length series)))
    (sqrt (/ (reduce #'+ (map 'list (lambda (x) (expt (- x m) 2)) series))
             (max 1 (1- n))))))

(defun autocorr (series mean var lag)
  "Lag-LAG autocorrelation of SERIES given its MEAN and VAR."
  (let ((n (length series)) (s 0.0d0))
    (loop for i from lag below n
          do (incf s (* (- (aref series i) mean)
                        (- (aref series (- i lag)) mean))))
    (/ s (* n var))))

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

(defun gelman-rubin (post name)
  "Gelman-Rubin R-hat for parameter NAME across chains. ~1.0 signals mixing."
  (let* ((j (post-index post name))
         (chs (posterior-chains post))
         (series (mapcar
                  (lambda (rows)
                    (let* ((nn (length rows)) (v (make-dvec nn)))
                      (dotimes (i nn v) (setf (aref v i) (aref (aref rows i) j)))))
                  chs))
         (m (length series))
         (n (reduce #'min (mapcar #'length series))))
    (if (or (< m 2) (< n 2))
        1.0d0
        (let* ((means (mapcar (lambda (s) (/ (loop for i below n sum (aref s i)) n))
                              series))
               (grand (/ (reduce #'+ means) m))
               (b (* (/ n (- m 1))
                     (reduce #'+ (mapcar (lambda (mu) (expt (- mu grand) 2)) means))))
               (w (/ (reduce #'+
                             (mapcar (lambda (s mu)
                                       (/ (loop for i below n
                                                sum (expt (- (aref s i) mu) 2))
                                          (- n 1)))
                                     series means))
                     m))
               (varhat (+ (* (/ (- n 1.0d0) n) w) (/ b n))))
          (if (<= w 0.0d0) 1.0d0 (sqrt (/ varhat w)))))))

(defun summarize (post &key (stream t))
  "Print a per-parameter posterior summary table with ESS and (if available)
   R-hat, plus the method's own diagnostic (acceptance rate or the VI fit)."
  (format stream "~%Posterior summary  [method ~(~a~), ~d draws, ~d chain(s)]~%"
          (posterior-method post) (length (posterior-draws post))
          (length (posterior-chains post)))
  (format stream "  ~8a ~10a ~10a ~10a ~10a ~10a ~8a~@[ ~7a~]~%"
          "param" "mean" "sd" "2.5%" "50%" "97.5%" "ess"
          (when (>= (length (posterior-chains post)) 2) "R-hat"))
  (loop for name across (posterior-param-names post) do
    (let* ((col (post-column post name))
           (sorted (sort (copy-seq col) #'<)))
      (format stream "  ~8a ~10,4f ~10,4f ~10,4f ~10,4f ~10,4f ~8,0f~@[ ~7,3f~]~%"
              name (mean-of col) (sd-of col)
              (quantile sorted 0.025d0) (quantile sorted 0.5d0)
              (quantile sorted 0.975d0) (ess col)
              (when (>= (length (posterior-chains post)) 2)
                (gelman-rubin post name)))))
  (let ((ar (getf (posterior-extra post) :accept-rate))
        (vm (getf (posterior-extra post) :vi-mean)))
    (when ar (format stream "  acceptance rate: ~,3f~%" ar))
    (when vm (format stream "  variational q fitted in unconstrained space~%")))
  post)

;; ----------------------------------------------------------------------------
;; 10. Terminal plots (histogram, trace, density, autocorrelation)
;; ----------------------------------------------------------------------------
;; The histogram binning follows chapter 07; here it also finds its own range.

(defun ascii-hist (series &key (bins 21) (width 48) (title ""))
  "Text histogram of SERIES scaled to WIDTH columns."
  (let* ((mn (reduce #'min series)) (mx (reduce #'max series)))
    (when (= mn mx) (setf mn (- mn 1.0d0) mx (+ mx 1.0d0)))
    (let* ((bw (/ (- mx mn) bins)) (counts (make-array bins :initial-element 0)))
      (loop for x across series do
        (incf (aref counts (min (1- bins) (max 0 (floor (/ (- x mn) bw)))))))
      (let ((mc (reduce #'max counts)))
        (unless (string= title "") (format t "~a~%" title))
        (dotimes (b bins)
          (format t "  ~9,3f | ~a~%"
                  (+ mn (* bw (+ b 0.5d0)))
                  (make-string (if (zerop mc) 0
                                   (round (* (/ (aref counts b) mc) width)))
                               :initial-element #\#)))))))

(defun ascii-line (series &key (rows 11) (cols 70) (title ""))
  "Text line plot of SERIES over its index (a trace plot)."
  (let ((n (length series)))
    (when (< n 2) (format t "~a: too few points~%" title) (return-from ascii-line))
    (let ((mn (reduce #'min series)) (mx (reduce #'max series)))
      (when (= mn mx) (setf mx (+ mn 1d-9)))
      (let ((grid (make-array (list rows cols) :initial-element #\Space)))
        (dotimes (c cols)
          (let* ((idx (min (1- n) (round (* c (/ (- n 1) (float (1- cols) 1.0d0))))))
                 (val (aref series idx))
                 (r (min (1- rows) (max 0 (round (* (/ (- val mn) (- mx mn))
                                                    (1- rows)))))))
            (setf (aref grid (- (1- rows) r) c) #\*)))
        (unless (string= title "") (format t "~a~%" title))
        (dotimes (r rows)
          (let ((line (make-string cols)))
            (dotimes (c cols) (setf (char line c) (aref grid r c)))
            (cond ((= r 0)         (format t "  ~9,3f |~a~%" mx line))
                  ((= r (1- rows)) (format t "  ~9,3f |~a~%" mn line))
                  (t               (format t "            |~a~%" line)))))))))

(defun plot-hist (post name &rest opts)
  "Posterior histogram for parameter NAME."
  (apply #'ascii-hist (post-column post name)
         :title (format nil "posterior of ~a" name) opts))

(defun plot-density (post name &rest opts)
  "Alias for plot-hist: a text density of parameter NAME."
  (apply #'plot-hist post name opts))

(defun plot-trace (post name &rest opts)
  "Trace plot of parameter NAME over sampling iterations."
  (apply #'ascii-line (post-column post name)
         :title (format nil "trace of ~a" name) opts))

(defun plot-autocorr (post name &key (max-lag 30))
  "Autocorrelation bars for parameter NAME up to MAX-LAG."
  (let* ((s (post-column post name)) (m (mean-of s))
         (var (/ (reduce #'+ (map 'list (lambda (x) (expt (- x m) 2)) s))
                 (length s))))
    (format t "autocorrelation of ~a~%" name)
    (loop for lag from 0 to max-lag
          for rho = (if (zerop var) 0.0d0 (autocorr s m var lag))
          do (format t "  lag ~3d ~6,3f |~a~%" lag rho
                     (make-string (max 0 (round (* (abs rho) 40)))
                                  :initial-element (if (>= rho 0) #\# #\-))))))

;; ----------------------------------------------------------------------------
;; 11. Example models and data
;; ----------------------------------------------------------------------------

(defmodel coin-model (flips)
  "Coin bias: p ~ Beta(1,1), each flip ~ Bernoulli(p)."
  (let ((p (sample :p (beta-dist 1.0d0 1.0d0))))
    (dolist (f flips) (observe (bernoulli p) f))))

(defmodel mean-var-model (data)
  "Unknown mean and spread: mu ~ Normal(0,10), sigma ~ HalfNormal(5),
   each datum ~ Normal(mu, sigma)."
  (let ((mu (sample :mu (normal 0.0d0 10.0d0)))
        (sigma (sample :sigma (half-normal 5.0d0))))
    (loop for x across data do (observe (normal mu sigma) x))))

(defmodel regression-model (xs ys)
  "Linear regression y ~ Normal(a + b x, sigma). The mean a + b*x is built
   with the differentiable ops so HMC and VI get gradients through it."
  (let ((a (sample :a (normal 0.0d0 10.0d0)))
        (b (sample :b (normal 0.0d0 10.0d0)))
        (sigma (sample :sigma (half-normal 5.0d0))))
    (loop for x across xs for y across ys
          do (observe (normal (g+ a (g* b x)) sigma) y))))

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

(defparameter *examples*
  '((:coin       "coin bias, Beta-Bernoulli (1 latent)")
    (:mean-var   "mean and sd of a Normal (2 latents)")
    (:regression "Bayesian linear regression (3 latents)"))
  "Registry of built-in example models for the REPL.")

(defun build-example (which)
  (ecase which
    (:coin (coin-model (example-coin-data)))
    (:mean-var (mean-var-model (example-normal-data)))
    (:regression (multiple-value-bind (xs ys) (example-regression-data)
                   (regression-model xs ys)))))

(defun run-example (which &rest infer-args)
  "Fit a built-in example. WHICH is :coin, :mean-var, or :regression;
   INFER-ARGS pass straight to `infer`, for example :method :hmc :chains 2."
  (apply #'infer (build-example which) infer-args))

;; ----------------------------------------------------------------------------
;; 12. Guided demonstration
;; ----------------------------------------------------------------------------

(defun demo ()
  "Run a short tour: fit the coin and mean/var models, plot posteriors, and
   compare the three engines on the coin bias against the exact answer."
  (format t "~%================ PPL DEMONSTRATION ================~%")

  (format t "~%--- 1. Coin bias with HMC (2 chains) ---~%")
  (format t "Data: 28 heads in 40 flips. Exact posterior is Beta(29,13),~%")
  (format t "so the posterior mean of p is 29/42 = ~,4f.~%" (/ 29.0d0 42.0d0))
  (let ((post (run-example :coin :method :hmc :chains 2
                           :iters 1500 :burn 500 :leapfrog 12 :step 0.18d0)))
    (summarize post)
    (plot-hist post :p :bins 19 :width 44)
    (plot-trace post :p :rows 9 :cols 64))

  (format t "~%--- 2. Unknown mean and sd with HMC ---~%")
  (format t "Data: 60 draws from Normal(5, 2).~%")
  (let ((post (run-example :mean-var :method :hmc :chains 2
                           :iters 2000 :burn 800 :leapfrog 18 :step 0.1d0)))
    (summarize post)
    (plot-density post :mu :bins 19 :width 44)
    (plot-density post :sigma :bins 19 :width 44))

  (format t "~%--- 3. Same coin, three engines vs the exact mean ---~%")
  (let ((truth (/ 29.0d0 42.0d0)))
    (format t "  ~14a ~10a ~10a~%" "engine" "E[p]" "abs err")
    (dolist (spec '((:mh  "Metropolis" (:iters 8000 :burn 3000 :step 0.6d0))
                    (:hmc "HMC"        (:iters 1500 :burn 500 :leapfrog 12 :step 0.18d0))
                    (:vi  "Variational" (:iters 2500 :step 0.05d0))))
      (destructuring-bind (method label opts) spec
        (let* ((post (apply #'run-example :coin :method method opts))
               (m (mean-of (post-column post :p))))
          (format t "  ~14a ~10,4f ~10,4f~%" label m (abs (- m truth)))))))

  (format t "~%Try it yourself. Examples:~%")
  (format t "  (defparameter *fit* (run-example :regression :method :hmc :chains 2))~%")
  (format t "  (summarize *fit*)~%")
  (format t "  (plot-hist *fit* :b)~%")
  (format t "===================================================~%")
  (values))

;; ----------------------------------------------------------------------------
;; 13. The command REPL
;; ----------------------------------------------------------------------------

(defun print-banner ()
  (format t "~%  Common Lisp Embedded Probabilistic Programming Language~%")
  (format t "  MCMC . Hamiltonian Monte Carlo . Variational Inference~%"))

(defun print-mini-help ()
  (format t "~%Type a command or any Lisp form. Commands:~%")
  (format t "  help    full help and worked examples~%")
  (format t "  models  list built-in example models~%")
  (format t "  demo    run the guided demonstration~%")
  (format t "  quit    leave the REPL~%")
  (format t "Quick start:  demo    then try  (run-example :coin :method :hmc)~%"))

(defun print-models ()
  (format t "~%Built-in example models (use with run-example):~%")
  (dolist (e *examples*)
    (format t "  ~12s ~a~%" (first e) (second e)))
  (format t "~%Fit one, for example:~%")
  (format t "  (run-example :mean-var :method :hmc :chains 2)~%"))

(defun print-help ()
  (print-banner)
  (format t "~%DEFINE A MODEL~%")
  (format t "  (defmodel my-model (data)~%")
  (format t "    (let ((mu (sample :mu (normal 0 10)))~%")
  (format t "          (s  (sample :s  (half-normal 5))))~%")
  (format t "      (loop for x across data do (observe (normal mu s) x))))~%")
  (format t "  sample declares a latent from a prior; observe scores data.~%")
  (format t "~%DISTRIBUTIONS~%")
  (format t "  normal half-normal exponential gamma-dist beta-dist~%")
  (format t "  uniform bernoulli poisson~%")
  (format t "  Combine latents in a mean with g+ g- g* g/ (they carry gradients).~%")
  (format t "~%FIT  (build a model by CALLING it with data, then infer)~%")
  (format t "  (infer (my-model data) :method :hmc :chains 2)~%")
  (format t "    :method  :hmc (default) | :mh | :vi~%")
  (format t "    :chains :iters :burn :step :leapfrog :thin :draws :mc~%")
  (format t "  (run-example :coin :method :vi)   fit a built-in model~%")
  (format t "~%INSPECT~%")
  (format t "  (summarize post)          table: mean sd quantiles ess R-hat~%")
  (format t "  (plot-hist post :mu)      posterior histogram~%")
  (format t "  (plot-density post :mu)   same, read as a density~%")
  (format t "  (plot-trace post :mu)     trace over iterations~%")
  (format t "  (plot-autocorr post :mu)  autocorrelation~%")
  (format t "  (post-column post :mu)    raw draws as a vector~%")
  (format t "~%EXAMPLE DATA GENERATORS~%")
  (format t "  (example-coin-data)        40 flips, 28 heads (list of 0/1)~%")
  (format t "  (example-normal-data)      60 draws from Normal(5, 2) (vector)~%")
  (format t "  (example-regression-data)  40 (x,y) points, y = 1 + 2x + noise~%")
  (format t "~%END-TO-END EXAMPLES  (copy, paste, run)~%")
  (format t "  ;; 1. Define a model, then fit and summarize in one line:~%")
  (format t "  (defmodel my-model (data)~%")
  (format t "    (let ((mu (sample :mu (normal 0 10)))~%")
  (format t "          (s  (sample :s  (half-normal 5))))~%")
  (format t "      (loop for x across data do (observe (normal mu s) x))))~%")
  (format t "  (summarize (infer (my-model (example-normal-data)) :method :hmc :chains 2))~%")
  (format t "~%  ;; 2. Keep the posterior in a variable, then inspect and plot:~%")
  (format t "  (defparameter *fit*~%")
  (format t "    (infer (my-model (example-normal-data)) :method :hmc :chains 2))~%")
  (format t "  (summarize *fit*)~%")
  (format t "  (plot-hist *fit* :mu)~%")
  (format t "  (plot-trace *fit* :s)~%")
  (format t "~%  ;; 3. Built-in models, one per engine:~%")
  (format t "  (summarize (run-example :coin :method :mh))~%")
  (format t "  (summarize (run-example :mean-var :method :hmc :chains 2))~%")
  (format t "  (summarize (run-example :regression :method :vi))~%")
  (format t "~%  ;; 4. Regression takes two data vectors:~%")
  (format t "  (defmodel line (xs ys)~%")
  (format t "    (let ((a (sample :a (normal 0 10)))~%")
  (format t "          (b (sample :b (normal 0 10)))~%")
  (format t "          (s (sample :s (half-normal 5))))~%")
  (format t "      (loop for x across xs for y across ys~%")
  (format t "            do (observe (normal (g+ a (g* b x)) s) y))))~%")
  (format t "  (multiple-value-bind (xs ys) (example-regression-data)~%")
  (format t "    (summarize (infer (line xs ys) :method :hmc :chains 2)))~%")
  (format t "~%REPL COMMANDS  help  models  demo  quit~%"))

(defun ppl-command-p (form)
  "True when FORM is a bare symbol naming a REPL command."
  (and (symbolp form) form
       (member (string-downcase (symbol-name form))
               '("help" "?" "h" "about" "models" "demo" "quit" "exit" "q")
               :test #'string=)))

(defun run-ppl-command (form)
  "Run the command named by symbol FORM. Returns :quit to leave, else t."
  (let ((head (string-downcase (symbol-name form))))
    (cond ((member head '("help" "?" "h") :test #'string=) (print-help) t)
          ((string= head "about") (print-banner) t)
          ((string= head "models") (print-models) t)
          ((string= head "demo") (demo) t)
          ((member head '("quit" "exit" "q") :test #'string=) :quit)
          (t t))))

(defun ppl-repl ()
  "Read-eval-print loop. It reads one whole Lisp form at a time, so a form
   that spans several lines (a pasted multi-line model, say) reads as a unit.
   A bare command word such as help, models, demo, or quit runs that command."
  (print-banner)
  (print-mini-help)
  (loop
    (format t "~&ppl> ")
    (finish-output)
    (handler-case
        (let ((form (read *standard-input* nil :eof)))
          (cond
            ((eq form :eof) (format t "~%bye~%") (return))
            ((ppl-command-p form)
             (when (eq (run-ppl-command form) :quit)
               (format t "bye~%") (return)))
            ;; A bare symbol that names a function (a model, say) but no
            ;; variable: Common Lisp is a Lisp-2, so guide instead of erroring.
            ((and (symbolp form) form (fboundp form) (not (boundp form)))
             (format t "~(~a~) names a function, not a variable. Call it with~%  arguments, for example  (~(~a~) data)  or  (~(~a~) (example-normal-data)),~%  then pass the result to infer.~%"
                     form form form))
            (t (let ((vals (multiple-value-list (eval form))))
                 (if vals
                     (format t "~{~s~^, ~}~%" vals)
                     (format t "; no value~%"))))))
      ;; End of stream (Ctrl-D or piped EOF) leaves the loop cleanly.
      (end-of-file () (format t "~%bye~%") (return))
      ;; A malformed form: report it and skip the rest of that line to resync.
      (reader-error (e)
        (format t "error: ~a~%" e)
        (ignore-errors (read-line *standard-input* nil nil)))
      (error (e) (format t "error: ~a~%" e)))))

;; ----------------------------------------------------------------------------
;; 14. Autostart
;; ----------------------------------------------------------------------------
;; Loading the file starts the REPL (set *ppl-autostart* to nil first to skip).
;; When the REPL ends we exit the host Lisp so the file behaves like an app;
;; a direct (ppl-repl) call just returns, which keeps it embeddable.

(when *ppl-autostart*
  (ppl-repl)
  #+sbcl (sb-ext:exit)
  #+lispworks (lispworks:quit)
  #+ccl (ccl:quit)
  #-(or sbcl lispworks ccl) (values))
