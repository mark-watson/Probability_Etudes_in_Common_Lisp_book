;;;; 06_law_of_large_numbers.lisp
;;;; Empirical demonstration of the (Weak) Law of Large Numbers.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; The LAW OF LARGE NUMBERS (LLN) is one of the cornerstones of probability.
;; Let X1, X2, ... be independent and identically distributed (i.i.d.) random
;; variables with finite mean mu = E[Xi]. Let the SAMPLE MEAN be
;;
;;     M_n = (X1 + ... + Xn) / n.
;;
;; By linearity of expectation, E[M_n] = mu for every n. The LLN says the
;; sample mean actually CONVERGES to mu as n grows:
;;
;;     WEAK LLN: for any epsilon > 0,
;;         P(|M_n - mu| > epsilon) -> 0   as n -> infinity.
;;
;; That is, the sample mean concentrates around the true mean: the probability
;; of a deviation larger than epsilon shrinks to zero.
;;
;;     STRONG LLN: M_n -> mu almost surely (with probability 1).
;;
;; This is why averaging repeated measurements is a sound estimation strategy:
;; noise averages out, and the average homes in on the expected value.
;;
;; Below we simulate rolling a fair die. The true mean is
;;     mu = (1+2+3+4+5+6)/6 = 3.5.
;; We track M_n for increasing n and watch it approach 3.5.
;; ============================================================================

(defvar *rng-state* (make-random-state t)
  "A freshly seeded random state so each run differs. We use the built-in
   pseudo-random number generator as a stand-in for true randomness.")

(defun roll-die ()
  "Simulate one roll of a fair six-sided die: uniform on {1,...,6}.
   Each face has probability 1/6."
  (+ 1 (random 6 *rng-state*)))

(defun sample-mean-of-rolls (n)
  "Compute M_n = average of n die rolls. By the LLN this approaches 3.5."
  (/ (loop for i below n sum (roll-die)) n))

(defun bernoulli-trial (p)
  "Simulate a Bernoulli(p): return 1 (success) with probability p, else 0.
   Uses the uniform draw U~Uniform(0,1); success when U < p."
  (if (< (random 1.0d0 *rng-state*) p) 1 0))

(defun sample-mean-bernoulli (p n)
  "Average of n Bernoulli(p) trials; estimates the success probability p.
   The LLN guarantees convergence to p."
  (/ (loop for i below n sum (bernoulli-trial p)) n))

(defun demonstrate-lln-die ()
  "Track the running sample mean for a fair die as n grows.
   True mean mu = 3.5. Watch M_n tighten around 3.5."
  (let ((mu 3.5))
    (format t "=== Law of Large Numbers: Fair Die (mu = ~a) ===~%" mu)
    (format t "       n     M_n     |M_n - mu|~%")
    (loop for n in '(10 100 1000 10000 100000 1000000) do
      (let ((mn (float (sample-mean-of-rolls n))))
        (format t "  ~8a  ~6,4f    ~6,4f~%" n mn (abs (- mn mu)))))))

(defun demonstrate-lln-bernoulli (p)
  "Track the running sample mean for Bernoulli(p). True mean = p."
  (let ((mu p))
    (format t "~%=== Law of Large Numbers: Bernoulli(p=~a) ===~%" p)
    (format t "       n     M_n     |M_n - p|~%")
    (loop for n in '(10 100 1000 10000 100000 1000000) do
      (let ((mn (float (sample-mean-bernoulli p n))))
        (format t "  ~8a  ~6,4f    ~6,4f~%" n mn (abs (- mn mu)))))))

(defun main ()
  (demonstrate-lln-die)
  (demonstrate-lln-bernoulli 0.3))

(main)
