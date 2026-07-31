;;;; 07_central_limit_theorem.lisp
;;;; Empirical demonstration of the Central Limit Theorem (CLT).

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; The CENTRAL LIMIT THEOREM (CLT) explains why the normal distribution
;; appears everywhere. Let X1, X2, ... be i.i.d. with mean mu and finite
;; variance sigma^2, and let S_n = X1 + ... + Xn. Then for large n,
;;
;;     (S_n - n mu) / (sigma sqrt(n))  --d-->  Normal(0, 1).
;;
;; Equivalently, the SAMPLE MEAN M_n = S_n/n is approximately Normal(mu,
;; sigma^2/n) for large n. Notice the variance of M_n shrinks as 1/n.
;;
;; Key points:
;;   - The original Xi can have ANY distribution (with finite variance); the
;;     SUM/MEAN still becomes approximately normal. This is remarkable.
;;   - The standardized sum converges to the STANDARD normal Z.
;;   - This is the theoretical basis for confidence intervals and z-tests.
;;
;; We demonstrate by drawing Bernoulli(0.5) samples (a very non-normal,
;; two-point distribution), forming sample means of size n, and checking that
;; the histogram of standardized means looks bell-shaped and that the
;; empirical mean/variance match the CLT predictions.
;; ============================================================================

(defvar *clt-random-state* (make-random-state t))

(defun bernoulli (p)
  "Draw from Bernoulli(p): 1 with probability p, else 0."
  (if (< (random 1.0d0 *clt-random-state*) p) 1 0))

(defun sample-mean (p n)
  "M_n = average of n i.i.d. Bernoulli(p). Mean mu=p, var sigma^2=p(1-p)."
  (/ (loop for i below n sum (bernoulli p)) n 1.0d0))

(defun collect-sample-means (p n num-samples)
  "Collect NUM-SAMPLES independent sample means, each from n Bernoulli(p)
   trials. This gives us an empirical distribution of M_n to compare with
   the CLT's predicted Normal(mu, sigma^2/n)."
  (loop for s below num-samples collect (sample-mean p n)))

(defun empirical-mean (data)
  "Arithmetic mean of a list of numbers."
  (/ (reduce #'+ data) (length data) 1.0d0))

(defun empirical-variance (data)
  "Sample variance (using division by n) of a list of numbers."
  (let ((mean (empirical-mean data))
        (n (length data)))
    (/ (reduce #'+ (mapcar (lambda (x) (expt (- x mean) 2)) data)) n 1.0d0)))

(defun make-histogram (data num-bins min-val max-val)
  "Count how many data points fall into each of NUM-BINS equal-width bins
   spanning [min-val, max-val]. Returns a list of (bin-center . count)."
  (let ((bin-width (/ (- max-val min-val) num-bins))
        (counts (make-array num-bins :initial-element 0)))
    (dolist (x data)
      (let* ((frac (/ (- x min-val) (- max-val min-val)))
             (bin (max 0 (min (1- num-bins) (floor (* frac num-bins))))))
        (incf (aref counts bin))))
    (loop for b below num-bins
          collect (cons (+ min-val (* bin-width (+ b 0.5)))
                        (aref counts b)))))

(defun print-histogram (bins max-bar-width)
  "Print a simple text histogram scaled to MAX-BAR-WIDTH characters."
  (let ((max-count (reduce #'max (mapcar #'cdr bins))))
    (dolist (b bins)
      (let* ((center (car b))
             (count (cdr b))
             ;; Give any nonzero bin at least one mark so the tapering tails
             ;; stay visible instead of printing as blank lines.
             (bar-len (if (zerop count) 0
                          (max 1 (round (* (/ count max-count) max-bar-width))))))
        (format t "  ~6,3f | ~a~%" center (make-string bar-len
                                                         :initial-element #\*))))))

(defun main ()
  (let ((p 0.5) (n 50) (num-samples 10000))
    (let* ((mu p)
           (sigma-sq (* p (- 1 p)))            ; Var of one Bernoulli
           (mean-of-mean mu)                    ; E[M_n] = mu
           (var-of-mean (/ sigma-sq n 1.0d0))   ; Var(M_n) = sigma^2/n
           (data (collect-sample-means p n num-samples))
           (emp-mean (empirical-mean data))
           (emp-var (empirical-variance data)))
      (format t "=== Central Limit Theorem ===~%")
      (format t "Source distribution: Bernoulli(p=~a), very non-normal.~%" p)
      (format t "Sample size n = ~a, number of sample means = ~a~%~%" n num-samples)
      (format t "CLT predicts for the sample mean M_n:~%")
      (format t "  E[M_n]        = mu        = ~4f~%" mean-of-mean)
      (format t "  Var(M_n)      = sigma^2/n = ~6,4f~%" var-of-mean)
      (format t "Empirical from simulation:~%")
      (format t "  mean of M_n's = ~6,4f~%" emp-mean)
      (format t "  var  of M_n's = ~6,4f~%" emp-var)
      (format t "  (These should be close to the CLT predictions.)~%~%")
      (format t "Histogram of sample means (bell-shaped = CLT at work):~%")
      ;; Range the histogram to the data so the tails do not waste empty bins.
      (let* ((lo (reduce #'min data))
             (hi (reduce #'max data))
             (bins (make-histogram data 20 lo hi)))
        (print-histogram bins 50))
      (format t "~%Notice the bell shape: sums of ANY i.i.d. variables with~%")
      (format t "finite variance tend toward a normal distribution.~%"))))

(main)
