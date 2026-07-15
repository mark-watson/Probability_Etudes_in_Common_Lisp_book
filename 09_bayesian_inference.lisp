;;;; 09_bayesian_inference.lisp
;;;; Bayesian inference: sequential belief updating for a coin bias.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; BAYESIAN INFERENCE treats an unknown parameter as a random variable and
;; updates a belief distribution over it as evidence arrives. The core engine
;; is Bayes' theorem applied to a parameter theta and observed data D:
;;
;;     p(theta | D)  ∝  p(D | theta) * p(theta)
;;        posterior      likelihood     prior
;;
;; The normalizing constant p(D) = ∫ p(D|theta) p(theta) dtheta makes the
;; posterior a proper density. Often we track the UNNORMALIZED posterior
;;   p(D|theta) p(theta)  and normalize at the end.
;;
;; CONJUGATE PRIORS: For certain likelihood/prior pairs, the posterior is in
;; the same family as the prior, which lets us update with simple arithmetic.
;; For a Bernoulli/binomial likelihood with unknown success probability theta,
;; the BETA(a, b) prior is conjugate:
;;
;;   Prior:        theta ~ Beta(a, b)
;;   Likelihood:   data = s successes in (s+f) trials, each ~ Bernoulli(theta)
;;   Posterior:    theta | data ~ Beta(a + s, b + f).
;;
;; The Beta(a,b) density on [0,1] is
;;   f(theta) = theta^(a-1) (1-theta)^(b-1) / B(a,b),
;; where B(a,b) is the Beta function (the normalizing constant). The
;; parameters a,b act as "pseudo-counts": a-1 prior successes, b-1 failures.
;;
;;   Posterior mean = (a + s) / (a + b + s + f).
;;   With a uniform prior Beta(1,1), this reduces to (s+1)/(s+f+2), the
;;   classic Laplace "rule of succession."
;;
;; As data accumulates (s+f -> infinity), the posterior concentrates around
;; the true theta and the prior's influence vanishes: this is BERNOULLI'S
;; THEOREM, a special case of the LLN applied to the posterior.
;; ============================================================================

(defun beta-mean (a b)
  "Posterior/prior mean of Beta(a,b): E[theta] = a / (a + b).
   This is the Bayesian point estimate under squared-error loss."
  (/ a (+ a b)))

(defun beta-variance (a b)
  "Var(theta) for Beta(a,b):
       a b / ((a+b)^2 (a+b+1)).
   As a+b grows (more data/pseudo-counts), the variance shrinks: the
   posterior tightens around the mean."
  (/ (* a b) (* (expt (+ a b) 2) (+ a b 1))))

(defun beta-pdf-unnormalized (a b theta)
  "Unnormalized Beta(a,b) density: theta^(a-1) (1-theta)^(b-1).
   We drop the Beta function B(a,b) since it does not depend on theta; it
   cancels in the posterior-to-prior ratio and is restored by normalizing."
  (* (expt theta (- a 1)) (expt (- 1 theta) (- b 1))))

(defun beta-pdf-normalized (a b theta)
  "Normalized Beta(a,b) PDF using the Gamma-function form of B(a,b):
       B(a,b) = Gamma(a) Gamma(b) / Gamma(a+b).
   We compute Gamma via Lisp's exp/gamma if available; otherwise approximate.
   Here we normalize numerically over [0,1] for portability."
  (let* ((grid 1000)
         (dx (/ 1.0d0 grid))
         (norm (loop for i below grid
                     sum (* dx (beta-pdf-unnormalized a b (* dx (+ i 0.5)))))))
    (/ (beta-pdf-unnormalized a b theta) norm 1.0d0)))

(defun bayesian-update (prior-a prior-b successes failures)
  "Conjugate update for a Beta prior with Bernoulli/binomial data.
   Prior  Beta(a, b)  ->  Posterior Beta(a + s, b + f).
   The data shifts the pseudo-counts: each observed success adds to a,
   each failure adds to b."
  (values (+ prior-a successes) (+ prior-b failures)))

(defun sequential-bayesian-update (prior-a prior-b data)
  "Process a sequence of observations (1 = success, 0 = failure) one at a
   time, updating the Beta posterior after each. Because Beta is conjugate,
   each step is just incrementing a or b. The result is the same as batching
   all data at once -- a key convenience of conjugacy."
  (let ((a prior-a) (b prior-b))
    (dolist (obs data)
      (if (= obs 1) (incf a) (incf b)))
    (values a b)))

(defun print-beta (a b label)
  "Print mean and variance of a Beta(a,b) distribution."
  (format t "  ~a: Beta(~a, ~a)  mean=~4f  var=~6,4f~%"
          label a b (float (beta-mean a b)) (float (beta-variance a b))))

(defun main ()
  (let ((prior-a 1) (prior-b 1))      ; Beta(1,1) = Uniform(0,1): maximum ignorance
    (format t "=== Bayesian Inference for a Coin Bias ===~%")
    (format t "Likelihood: Bernoulli(theta). Prior: Beta(1,1) = Uniform on [0,1].~%")
    (format t "True theta (used only to generate data) = 0.7.~%~%")
    (print-beta prior-a prior-b "Prior         ")
    ;; Simulate data from a coin with true bias 0.7.
    (let* ((rs (make-random-state t))
          (data (loop for i below 100
                      collect (if (< (random 1.0d0 rs) 0.7d0) 1 0))))
      (let ((successes (count 1 data))
            (failures (count 0 data)))
        (format t "~%Observed ~a successes and ~a failures in 100 flips.~%"
                successes failures)
        ;; Batch update.
        (multiple-value-bind (pa pb) (bayesian-update prior-a prior-b
                                                       successes failures)
          (print-beta pa pb "Posterior     "))
        ;; Sequential update in chunks to show belief concentrating.
        (format t "~%Sequential updates (belief concentrates as data arrives):~%")
        (let ((a prior-a) (b prior-b))
          (dolist (obs data)
            (if (= obs 1) (incf a) (incf b)))
          (print-beta a b "After all 100 "))
        ;; Show intermediate stages.
        (let ((a prior-a) (b prior-b) (count 0))
          (dolist (obs data)
            (if (= obs 1) (incf a) (incf b))
            (incf count)
            (when (member count '(10 50 100))
              (print-beta a b (format nil "After ~4a flips" count)))))))
    (format t "~%True theta = 0.7. As n grows the posterior mean converges to the~%")
    (format t "true value and the variance shrinks to 0 (Bernoulli's theorem).~%")))

(main)
