;;;; 05_continuous_distributions.lisp
;;;; Continuous random variables: PDF, CDF, expectation for
;;;; Uniform, Exponential, and Normal distributions.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; A CONTINUOUS random variable can take uncountably many values (e.g. any
;; real number in an interval). Because P(X = x) = 0 for any single point,
;; we cannot use a PMF. Instead we use the PROBABILITY DENSITY FUNCTION (PDF)
;; f_X(x). Probability is the AREA under the density:
;;
;;     P(a <= X <= b) = ∫_a^b f_X(x) dx.
;;
;; A valid PDF satisfies: f_X(x) >= 0 and ∫_{-∞}^{∞} f_X(x) dx = 1.
;; The PDF itself is not a probability; it can exceed 1. Only the integral
;; over an interval is a probability.
;;
;; The CDF is F_X(x) = P(X <= x) = ∫_{-∞}^{x} f_X(t) dt, and f_X = F'_X.
;;
;; Expectation and variance for continuous variables:
;;     E[X]     = ∫_{-∞}^{∞} x f_X(x) dx
;;     E[g(X)]  = ∫_{-∞}^{∞} g(x) f_X(x) dx
;;     Var(X)   = E[X^2] - (E[X])^2.
;;
;; -- UNIFORM(a, b): constant density 1/(b-a) on [a,b], zero elsewhere.
;;    f(x) = 1/(b-a) for a <= x <= b.
;;    E[X] = (a+b)/2,  Var(X) = (b-a)^2 / 12.
;;
;; -- EXPONENTIAL(lambda): models waiting time in a Poisson process.
;;    f(x) = lambda e^{-lambda x} for x >= 0.
;;    CDF: F(x) = 1 - e^{-lambda x}.
;;    E[X] = 1/lambda,  Var(X) = 1/lambda^2.
;;    It is MEMORYLESS (continuous analogue of the geometric): the waiting
;;    time remaining is always exponential, regardless of how long you waited.
;;
;; -- NORMAL(mu, sigma^2): the bell curve.
;;    f(x) = (1/(sigma sqrt(2π))) exp(-(x-mu)^2 / (2 sigma^2)).
;;    E[X] = mu,  Var(X) = sigma^2.
;;    The STANDARD normal has mu=0, sigma=1. The 68-95-99.7 rule says
;;    ~68% of mass lies within 1 sigma, ~95% within 2 sigma, ~99.7% within 3.
;;
;; Because Lisp's numeric integration is simple, we approximate integrals
;; numerically with the rectangle (midpoint) rule.
;; ============================================================================

(defun integrate-rectangle (f a b &optional (n 100000))
  "Approximate ∫_a^b f(x) dx via the midpoint rule with N subintervals.
   This is the discrete stand-in for the continuous integral. Accuracy grows
   as 1/N; midpoint rule error is O(1/N^2) for smooth functions."
  (let ((h (/ (- b a) n)))
    (* h (loop for i from 0 below n
               sum (funcall f (+ a (* h (+ i 1/2))))))))

(defun uniform-pdf (a b x)
  "PDF of Uniform(a,b): constant 1/(b-a) on [a,b]."
  (if (and (>= x a) (<= x b))
      (/ 1.0 (- b a))
      0.0))

(defun uniform-cdf (a b x)
  "CDF of Uniform(a,b): 0 below a, (x-a)/(b-a) on [a,b], 1 above b."
  (cond ((< x a) 0.0)
        ((> x b) 1.0)
        (t (/ (- x a) (- b a)))))

(defun uniform-mean (a b)
  "E[X] = (a+b)/2 for Uniform(a,b), by symmetry of the constant density."
  (/ (+ a b) 2.0))

(defun uniform-variance (a b)
  "Var(X) = (b-a)^2/12 for Uniform(a,b). The 1/12 factor comes from
   integrating x^2 over [a,b] against the flat density."
  (/ (expt (- b a) 2) 12.0))

(defun exponential-pdf (lambda-rate x)
  "PDF of Exponential(lambda): lambda e^{-lambda x} for x >= 0.
   The density decays geometrically fast; most mass is near 0."
  (if (>= x 0)
      (* lambda-rate (exp (- (* lambda-rate x))))
      0.0))

(defun exponential-cdf (lambda-rate x)
  "CDF of Exponential(lambda): 1 - e^{-lambda x} for x >= 0."
  (if (>= x 0)
      (- 1.0 (exp (- (* lambda-rate x))))
      0.0))

(defun exponential-mean (lambda-rate)
  "E[X] = 1/lambda. Higher rate => shorter expected wait."
  (/ 1.0 lambda-rate))

(defun exponential-variance (lambda-rate)
  "Var(X) = 1/lambda^2."
  (/ 1.0 (* lambda-rate lambda-rate)))

(defun normal-pdf (mu sigma x)
  "PDF of Normal(mu, sigma^2): the bell curve.
   (1/(sigma sqrt(2π))) exp(-(x-mu)^2/(2 sigma^2))."
  (let ((z (/ (- x mu) sigma)))
    (/ (exp (- (/ (* z z) 2.0)))
       (* sigma (sqrt (* 2.0 pi))))))

(defun standard-normal-cdf (x)
  "CDF of the standard normal via a numerical approximation (Abramowitz &
   Stegun 26.2.17). This avoids needing the non-elementary error function;
   the absolute error is < 7.5e-8."
  (let* ((a1 0.254829592d0) (a2 -0.284496736d0) (a3 1.421413741d0)
         (a4 -1.453152027d0) (a5 1.061405429d0) (p 0.3275911d0)
         (sign (if (>= x 0) 1 -1))
         (ax (abs x))
         (t-val (/ 1.0 (+ 1.0 (* p ax))))
         (y (* t-val (+ a1 (* t-val (+ a2 (* t-val (+ a3 (* t-val (+ a4 (* t-val a5)))))))))))
    ;; y = polynomial in t; erf(|x|) ≈ 1 - y*exp(-x^2)
    ;; Phi(x) = 0.5*(1 + sign(x)*erf(|x|))
    (let ((erf (- 1.0 (* y (exp (- (* ax ax)))))))
      (* 0.5 (+ 1.0 (* sign erf))))))

(defun main ()
  (format t "=== Uniform(0, 2) ===~%")
  (let ((a 0.0) (b 2.0))
    (format t "  E[X] (formula) = ~4f~%" (uniform-mean a b))
    (format t "  E[X] (numeric) = ~4f~%"
            (integrate-rectangle (lambda (x) (* x (uniform-pdf a b x))) a b))
    (format t "  Var(X) (formula) = ~4f~%" (uniform-variance a b))
    (format t "  P(0.5 <= X <= 1.5) = ~4f (exact 0.5)~%"
            (integrate-rectangle (lambda (x) (uniform-pdf a b x)) 0.5 1.5)))

  (format t "~%=== Exponential(lambda=2) ===~%")
  (let ((rate 2.0))
    (format t "  E[X] (formula) = ~4f~%" (exponential-mean rate))
    (format t "  E[X] (numeric) = ~4f~%"
            (integrate-rectangle (lambda (x) (* x (exponential-pdf rate x))) 0 10))
    (format t "  Var(X) (formula) = ~4f~%" (exponential-variance rate))
    (format t "  P(X <= 1) = ~4f (CDF) vs ~4f (numeric)~%"
            (exponential-cdf rate 1.0)
            (integrate-rectangle (lambda (x) (exponential-pdf rate x)) 0 1))
    (format t "  Total probability (should be 1.0): ~4f~%"
            (integrate-rectangle (lambda (x) (exponential-pdf rate x)) 0 20)))

  (format t "~%=== Standard Normal (mu=0, sigma=1) ===~%")
  (let ((mu 0.0) (sigma 1.0))
    (format t "  E[X] (numeric) = ~6f (should be 0)~%"
            (integrate-rectangle (lambda (x) (* x (normal-pdf mu sigma x))) -6 6))
    (format t "  E[X^2] (numeric) = ~6f (should be 1 = Var)~%"
            (integrate-rectangle (lambda (x) (* x x (normal-pdf mu sigma x))) -6 6))
    (format t "  Total probability (should be 1.0): ~4f~%"
            (integrate-rectangle (lambda (x) (normal-pdf mu sigma x)) -6 6))
    (format t "  68%% rule P(-1 <= Z <= 1) = ~4f (should be ~~0.6827)~%"
            (integrate-rectangle (lambda (x) (normal-pdf mu sigma x)) -1 1))
    (format t "  95%% rule P(-2 <= Z <= 2) = ~4f (should be ~~0.9545)~%"
            (integrate-rectangle (lambda (x) (normal-pdf mu sigma x)) -2 2))
    (format t "  99.7%% rule P(-3 <= Z <= 3) = ~4f (should be ~~0.9973)~%"
            (integrate-rectangle (lambda (x) (normal-pdf mu sigma x)) -3 3))
    (format t "  CDF(0) = ~4f (should be 0.5 by symmetry)~%"
            (standard-normal-cdf 0.0))))

(main)
