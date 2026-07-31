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

(defun standardize (p)
  "PMF of the standardized variable Z = (X - mu)/sigma, which has mean 0 and
   variance 1. Each value x maps to (x - mu)/sigma; probabilities are unchanged.
   Standardization puts different distributions on a common scale."
  (let ((mu (expectation p))
        (sigma (standard-deviation p)))
    (make-pmf :table (mapcar (lambda (entry)
                               (cons (/ (- (car entry) mu) sigma) (cdr entry)))
                             (pmf-table p)))))

(defun convolve-pmfs (p q)
  "PMF of the sum X + Y for INDEPENDENT X ~ p and Y ~ q:
     P(X + Y = s) = sum_{x + y = s} p(x) q(y).
   Independence is what lets the joint probability factor into p(x) q(y)."
  (let ((table nil))
    (dolist (ex (pmf-table p))
      (dolist (ey (pmf-table q))
        (let* ((s (+ (car ex) (car ey)))
               (pr (* (cdr ex) (cdr ey)))
               (cell (assoc s table :test #'eql)))
          (if cell (incf (cdr cell) pr) (push (cons s pr) table)))))
    (make-pmf :table (sort table #'< :key #'car))))

(defun joint-expectation (joint fn)
  "E[fn(x, y)] for a JOINT PMF represented as an alist mapping (x y) -> prob."
  (reduce #'+ (mapcar (lambda (entry)
                        (destructuring-bind (x y) (car entry)
                          (* (funcall fn x y) (cdr entry))))
                      joint)))

(defun covariance (joint)
  "Cov(X, Y) = E[XY] - E[X] E[Y] for a joint PMF over (x y) pairs.
   Zero for independent variables; nonzero signals linear co-movement."
  (let ((ex  (joint-expectation joint (lambda (x y) (declare (ignore y)) x)))
        (ey  (joint-expectation joint (lambda (x y) (declare (ignore x)) y)))
        (exy (joint-expectation joint (lambda (x y) (* x y)))))
    (- exy (* ex ey))))

(defun correlation (joint)
  "Correlation coefficient rho = Cov(X,Y)/(sigma_X sigma_Y), always in [-1, 1].
   It measures only the LINEAR part of the association between X and Y."
  (let ((varx (- (joint-expectation joint (lambda (x y) (declare (ignore y)) (* x x)))
                 (expt (joint-expectation joint (lambda (x y) (declare (ignore y)) x)) 2)))
        (vary (- (joint-expectation joint (lambda (x y) (declare (ignore x)) (* y y)))
                 (expt (joint-expectation joint (lambda (x y) (declare (ignore x)) y)) 2))))
    (/ (covariance joint) (* (sqrt varx) (sqrt vary)))))

(defun chebyshev-tail (p k)
  "Actual P(|X - mu| >= k sigma) for this PMF, to compare against Chebyshev's
   distribution-free bound 1/k^2."
  (let ((mu (expectation p))
        (sigma (standard-deviation p)))
    (reduce #'+ (mapcar (lambda (entry)
                          (if (>= (abs (- (car entry) mu)) (* k sigma))
                              (cdr entry) 0))
                        (pmf-table p)))))

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
    (print-distribution dice-sum "DiceSum")

    ;; Build the same distribution by CONVOLVING two single dice.
    (format t "~%=== Convolution: single die + single die ===~%")
    (let* ((die (make-pmf-from-alist (loop for i from 1 to 6 collect (cons i 1/6))))
           (sum2 (convolve-pmfs die die)))
      (format t "  E[die+die]   = ~a = ~4f (matches DiceSum mean 7)~%"
              (expectation sum2) (float (expectation sum2)))
      (format t "  Var(die+die) = ~a = ~4f (matches DiceSum var 35/6)~%"
              (variance sum2) (float (variance sum2))))

    ;; STANDARDIZE the dice sum: the result should have mean 0, variance 1.
    (format t "~%=== Standardization: Z = (X - mu)/sigma ===~%")
    (let ((z (standardize dice-sum)))
      (format t "  mean(Z) = ~6,3f (should be 0)~%" (float (expectation z)))
      (format t "  var(Z)  = ~6,3f (should be 1)~%" (float (variance z))))

    ;; CHEBYSHEV: actual tail vs the distribution-free bound 1/k^2.
    (format t "~%=== Chebyshev's Inequality (dice sum) ===~%")
    (format t "  k     actual P(|X-mu|>=k sigma)   bound 1/k^2~%")
    (dolist (k '(1.0 1.5 2.0 2.5))
      (format t "  ~3,1f   ~18,4f        ~10,4f~%"
              k (float (chebyshev-tail dice-sum k)) (/ 1.0 (* k k)))))

  ;; COVARIANCE and CORRELATION need a JOINT distribution over pairs.
  (format t "~%=== Covariance and Correlation ===~%")
  (let ((indep (loop for i from 1 to 6 append
                  (loop for j from 1 to 6 collect (cons (list i j) 1/36))))
        (x-and-sum (loop for i from 1 to 6 append
                     (loop for j from 1 to 6 collect (cons (list i (+ i j)) 1/36)))))
    (format t "  Two independent dice (X, Y):~%")
    (format t "    Cov = ~6,3f  rho = ~6,3f (independent => 0)~%"
            (float (covariance indep)) (float (correlation indep)))
    (format t "  A die and the two-dice sum (X, X+Y):~%")
    (format t "    Cov = ~6,3f  rho = ~6,3f (rho = 1/sqrt(2) ~~ 0.707)~%"
            (float (covariance x-and-sum)) (float (correlation x-and-sum)))))

(main)
