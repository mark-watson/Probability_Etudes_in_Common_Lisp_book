;;;; 03_discrete_random_variables.lisp
;;;; Discrete random variables: PMF, CDF, expectation, and variance.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; A RANDOM VARIABLE X is a function from the sample space Omega to the real
;; numbers. It assigns a numerical value to each outcome. A DISCRETE random
;; variable takes values in a countable set (often a finite set).
;;
;; The PROBABILITY MASS FUNCTION (PMF) gives the probability that X equals a
;; specific value x:
;;
;;     p_X(x) = P(X = x) = P({ omega : X(omega) = x }).
;;
;; The PMF satisfies:  p_X(x) >= 0  and  Σ_x p_X(x) = 1.
;;
;; The CUMULATIVE DISTRIBUTION FUNCTION (CDF) accumulates probability up to a
;; value:
;;
;;     F_X(x) = P(X <= x) = Σ_{t <= x} p_X(t).
;;
;; The CDF is non-decreasing, right-continuous, with limits 0 at -∞ and 1 at +∞.
;;
;; The EXPECTED VALUE (mean) is the probability-weighted average:
;;
;;     E[X] = Σ_x x * p_X(x).
;;
;; It is the long-run average of X if the experiment is repeated many times.
;; Two key properties (linearity of expectation):
;;     E[aX + b]        = a E[X] + b
;;     E[X + Y]         = E[X] + E[Y]      (even if X and Y are dependent!)
;;
;; The VARIANCE measures spread around the mean:
;;
;;     Var(X) = E[(X - E[X])^2] = E[X^2] - (E[X])^2.
;;
;; The STANDARD DEVIATION is sigma = sqrt(Var(X)).
;;
;; ============================================================================

(defstruct pmf
  "A discrete probability mass function represented as an alist
   mapping value -> probability."
  table)

(defun make-pmf-from-alist (alist)
  "Build a PMF from an alist of (value . probability) pairs."
  (make-pmf :table alist))

(defun pmf-lookup (p x)
  "Return p_X(x), the probability that X = x. Defaults to 0 if x not listed."
  (let ((entry (assoc x (pmf-table p) :test #'equal)))
    (if entry (cdr entry) 0)))

(defun pmf-values (p)
  "Return all values x for which p_X(x) > 0."
  (mapcar #'car (pmf-table p)))

(defun pmf-total-probability (p)
  "Verify the normalization axiom: Σ_x p_X(x) should equal 1."
  (reduce #'+ (mapcar #'cdr (pmf-table p))))

(defun cdf (p x)
  "Cumulative distribution function F_X(x) = P(X <= x) = Σ_{t <= x} p_X(t).
   Values are compared with <=, so this works for numbers."
  (reduce #'+ (mapcar (lambda (entry)
                        (if (<= (car entry) x) (cdr entry) 0))
                      (pmf-table p))))

(defun expectation (p)
  "E[X] = Σ_x x * p_X(x), the mean of the distribution."
  (reduce #'+ (mapcar (lambda (entry) (* (car entry) (cdr entry)))
                      (pmf-table p))))

(defun expectation-of-square (p)
  "E[X^2] = Σ_x x^2 * p_X(x), needed for the variance shortcut."
  (reduce #'+ (mapcar (lambda (entry)
                        (* (car entry) (car entry) (cdr entry)))
                      (pmf-table p))))

(defun variance (p)
  "Var(X) = E[X^2] - (E[X])^2.
   This identity is equivalent to the definition E[(X - mu)^2] but cheaper:
   expanding (X - mu)^2 = X^2 - 2 mu X + mu^2 and using E[X]=mu gives the form."
  (let ((ex (expectation p)))
    (- (expectation-of-square p) (* ex ex))))

(defun standard-deviation (p)
  "sigma = sqrt(Var(X))."
  (sqrt (variance p)))

(defun print-distribution (p &optional (name "X"))
  "Pretty-print the PMF, CDF at each value, mean, and variance."
  (format t "Distribution of ~a:~%" name)
  (format t "  Normalization check Σ p(x) = ~a~%" (pmf-total-probability p))
  (dolist (x (sort (pmf-values p) #'<))
    (format t "  x=~a : P(X=x)=~a  P(X<=x)=~a~%"
            x (pmf-lookup p x) (cdf p x)))
  (format t "  E[~a]   = ~a = ~4f~%" name (expectation p) (float (expectation p)))
  (format t "  Var(~a) = ~a = ~4f~%" name (variance p) (float (variance p)))
  (format t "  sigma   = ~4f~%" (float (standard-deviation p))))

(defun main ()
  ;; A loaded die: the probability of rolling a 6 is 1/2, the rest share 1/10.
  (let ((loaded-die
          (make-pmf-from-alist
            (list (cons 1 1/10) (cons 2 1/10) (cons 3 1/10)
                  (cons 4 1/10) (cons 5 1/10) (cons 6 1/2)))))
    (format t "=== Loaded Die ===~%")
    (print-distribution loaded-die "LoadedDie"))
  (format t "~%")
  ;; A simple Bernoulli distribution: X = 1 with prob p, 0 with prob 1-p.
  (let ((bernoulli-p-half
          (make-pmf-from-alist (list (cons 0 1/2) (cons 1 1/2)))))
    (format t "=== Bernoulli(p=1/2) ===~%")
    (print-distribution bernoulli-p-half "Bernoulli")
    (format t "  (A Bernoulli(p) has mean p and variance p(1-p).)~%"))
  (format t "~%")
  ;; Sum of two fair dice as a random variable.
  (let ((dice-sum
          (make-pmf-from-alist
            (list (cons 2 1/36) (cons 3 2/36) (cons 4 3/36) (cons 5 4/36)
                  (cons 6 5/36) (cons 7 6/36) (cons 8 5/36) (cons 9 4/36)
                  (cons 10 3/36) (cons 11 2/36) (cons 12 1/36)))))
    (format t "=== Sum of Two Fair Dice ===~%")
    (print-distribution dice-sum "DiceSum")))

(main)
