;;;; 08_monte_carlo.lisp
;;;; Monte Carlo estimation: estimating pi by random sampling.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; MONTE CARLO METHODS estimate quantities by repeated random sampling. The
;; theoretical justification is again the Law of Large Numbers: if we can
;; express a quantity as the expectation of a random variable, then the
;; sample average of repeated draws converges to that expectation.
;;
;; Here we estimate pi. Consider a unit square [0,1]x[0,1] (area 1) and the
;; quarter disk of radius 1 centered at the origin (area pi/4). If we throw
;; uniformly random points into the square, the probability a point lands in
;; the quarter disk equals the ratio of areas:
;;
;;     P(point in disk) = (pi/4) / 1 = pi/4.
;;
;; Let X_i = 1 if the i-th point lands in the disk, 0 otherwise. Then the Xi
;; are i.i.d. Bernoulli(pi/4), and by the LLN:
;;
;;     (1/n) Σ X_i  ->  pi/4,   so   4 * (1/n) Σ X_i  ->  pi.
;;
;; The ESTIMATOR is  4 * (fraction of points inside the disk).
;; Its variance is Var(4 * M_n) = 16 * (pi/4)(1 - pi/4) / n, which shrinks
;; as 1/n. Convergence is slow (1/sqrt(n) in standard deviation), which is
;; the classic downside of Monte Carlo: halving the error needs 4x the samples.
;; ============================================================================

(defvar *mc-random-state* (make-random-state t))

(defun random-point-in-unit-square ()
  "Return (x, y) with x,y ~ Uniform(0,1) independently.
   Uniform sampling over the square gives each region a probability equal to
   its area."
  (values (random 1.0d0 *mc-random-state*)
          (random 1.0d0 *mc-random-state*)))

(defun in-quarter-disk-p (x y)
  "Is (x,y) inside the quarter disk of radius 1 centered at the origin?
   The condition is x^2 + y^2 <= 1 (Pythagorean distance)."
  (<= (+ (* x x) (* y y)) 1.0d0))

(defun estimate-pi (n)
  "Monte Carlo estimate of pi using n random points.
   Estimator = 4 * (count inside disk) / n.
   By the LLN this converges to 4 * (pi/4) = pi."
  (let ((inside 0))
    (dotimes (i n)
      (multiple-value-bind (x y) (random-point-in-unit-square)
        (when (in-quarter-disk-p x y)
          (incf inside))))
    (* 4.0d0 (/ inside n 1.0d0))))

(defun estimate-pi-with-se (n)
  "Estimate pi and also report a STANDARD ERROR for the estimate.
   For a Bernoulli(pi/4) variable, the standard error of M_n is
   sqrt(p(1-p)/n). The estimator is 4*M_n, so its SE is 4*sqrt(p(1-p)/n).
   We plug in the sample proportion p_hat for the unknown p."
  (let* ((inside 0)
         (estimate 4.0d0)
         (p-hat 0.0d0)
         (se 0.0d0))
    (dotimes (i n)
      (multiple-value-bind (x y) (random-point-in-unit-square)
        (when (in-quarter-disk-p x y)
          (incf inside))))
    (setf p-hat (/ inside n 1.0d0))
    (setf estimate (* 4.0d0 p-hat))
    ;; SE of 4*M_n = 4 * sqrt( p_hat (1 - p_hat) / n )
    (setf se (* 4.0d0 (sqrt (/ (* p-hat (- 1.0d0 p-hat)) n 1.0d0))))
    (values estimate se)))

(defun estimate-pi-antithetic (n)
  "Estimate pi with ANTITHETIC VARIATES. For each random point (u, v) we also
   score its reflection (1-u, 1-v). The two disk indicators are negatively
   correlated, so their average has lower variance than two independent draws.
   Uses n pairs, i.e. 2 n points, so it matches (estimate-pi (* 2 n)) in cost."
  (let ((inside 0))
    (dotimes (i n)
      (let ((u (random 1.0d0 *mc-random-state*))
            (v (random 1.0d0 *mc-random-state*)))
        (when (in-quarter-disk-p u v) (incf inside))
        (when (in-quarter-disk-p (- 1.0d0 u) (- 1.0d0 v)) (incf inside))))
    (* 4.0d0 (/ inside (* 2 n) 1.0d0))))

(defun estimate-pi-buffon (n &optional (needle 1.0d0) (spacing 1.0d0))
  "Buffon's needle estimate of pi (1777), an older Monte Carlo estimator. Drop
   N needles of length NEEDLE on a floor ruled with lines SPACING apart
   (needle <= spacing). Let y ~ Uniform(0, spacing/2) be the distance from a
   needle's center to the nearest line and theta ~ Uniform(0, pi/2) its angle;
   the needle crosses a line when (needle/2) sin(theta) >= y. Since
   P(cross) = 2 needle / (pi spacing), pi ~ 2 needle n / (spacing crossings)."
  (let ((crossings 0))
    (dotimes (i n)
      (let ((y (* (/ spacing 2.0d0) (random 1.0d0 *mc-random-state*)))
            (theta (* (/ pi 2.0d0) (random 1.0d0 *mc-random-state*))))
        (when (>= (* (/ needle 2.0d0) (sin theta)) y)
          (incf crossings))))
    (if (zerop crossings)
        0.0d0
        (/ (* 2.0d0 needle n) (* spacing crossings)))))

(defun repeat-estimator (thunk reps)
  "Run THUNK REPS times; return (values mean empirical-standard-deviation) of
   its outputs. This measures an estimator's spread directly, which lets us
   compare variance-reduction schemes at equal cost."
  (let* ((xs (loop repeat reps collect (funcall thunk)))
         (m (/ (reduce #'+ xs) reps))
         (v (/ (reduce #'+ (mapcar (lambda (x) (expt (- x m) 2)) xs)) reps)))
    (values m (sqrt v))))

(defun main ()
  (format t "=== Monte Carlo Estimation of pi ===~%")
  (format t "True pi = ~6,6f~%~%" pi)
  (format t "       n     estimate     error    std-error~%")
  (loop for n in '(1000 10000 100000 1000000) do
    (multiple-value-bind (est se) (estimate-pi-with-se n)
      (format t "  ~8a  ~8,6f  ~8,6f  ~8,6f~%"
              n est (abs (- est pi)) se)))
  (format t "~%Note: the error scales roughly as 1/sqrt(n). Doubling n~%")
  (format t "shrinks the standard error by ~~0.707x; halving the error needs 4x n.~%")

  (format t "~%=== Variance Reduction and Buffon's Needle ===~%")
  (let ((reps 200) (n 5000))
    (format t "  Comparing estimators over ~a repetitions of ~a points each:~%" reps n)
    (multiple-value-bind (m sd) (repeat-estimator (lambda () (estimate-pi n)) reps)
      (format t "    plain disk       mean=~7,4f  empirical SD=~7,4f~%" m sd))
    (multiple-value-bind (m sd)
        (repeat-estimator (lambda () (estimate-pi-antithetic (floor n 2))) reps)
      (format t "    antithetic disk  mean=~7,4f  empirical SD=~7,4f~%" m sd))
    (multiple-value-bind (m sd) (repeat-estimator (lambda () (estimate-pi-buffon n)) reps)
      (format t "    Buffon's needle  mean=~7,4f  empirical SD=~7,4f~%" m sd)))
  (format t "  Antithetic variates lower the SD versus plain sampling at equal cost;~%")
  (format t "  Buffon's needle is a valid but higher-variance estimator.~%"))

(main)
