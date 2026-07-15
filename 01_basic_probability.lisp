;;;; 01_basic_probability.lisp
;;;; Sample spaces, events, and the axioms of probability.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; A RANDOM EXPERIMENT is a procedure whose outcome is not deterministic.
;; The SAMPLE SPACE (written Omega) is the set of all possible outcomes.
;; An EVENT is any subset of the sample space (in the discrete case, any
;;   collection of outcomes).
;;
;; A PROBABILITY MEASURE P is a function from events to [0,1] satisfying
;; the three KOLMOGOROV AXIOMS:
;;   1. Non-negativity:    P(A) >= 0 for every event A.
;;   2. Normalization:     P(Omega) = 1.
;;   3. Countable additivity: for pairwise-disjoint events A1, A2, ...
;;      P(A1 u A2 u ...) = sum P(Ai).
;;
;; When outcomes are EQUALLY LIKELY (e.g. a fair die), the probability of an
;; event A reduces to the classical definition:
;;
;;     P(A) = |A| / |Omega|     (favorable outcomes / total outcomes)
;;
;; This file explores the sample space of rolling two dice and asks:
;;   - What is the probability the sum equals 7?
;;   - What is the probability at least one die shows a 6?
;;   - What is the probability the sum is >= 10?
;; ============================================================================

(defun make-two-dice-sample-space ()
  "Build the sample space Omega for rolling two distinguishable dice.
   Omega = { (i,j) : 1 <= i <= 6, 1 <= j <= 6 } has |Omega| = 36 outcomes."
  (loop for i from 1 to 6
        append (loop for j from 1 to 6 collect (list i j))))

(defun event-sum-is (omega target)
  "Event: outcomes whose dice sum equals TARGET.
   This is the subset { (i,j) in Omega : i + j = TARGET }."
  (remove-if-not (lambda (outcome) (= (reduce #'+ outcome) target)) omega))

(defun event-at-least-one-is (omega value)
  "Event: outcomes where at least one die equals VALUE.
   Equivalent to the union of 'first die = VALUE' and 'second die = VALUE'."
  (remove-if-not (lambda (outcome) (member value outcome)) omega))

(defun event-sum-at-least (omega threshold)
  "Event: outcomes whose dice sum is >= THRESHOLD."
  (remove-if-not (lambda (outcome) (>= (reduce #'+ outcome) threshold)) omega))

(defun probability-classical (event omega)
  "Compute P(event) under the equally-likely-outcomes model:
       P(A) = |A| / |Omega|.
   This respects the axioms: probabilities lie in [0,1], P(Omega)=1, and
   the measure is additive on disjoint events."
  (/ (length event) (length omega)))

(defun demonstrate-axiom-normalization (omega)
  "Verify Kolmogorov's second axiom: P(Omega) = 1."
  (probability-classical omega omega))

(defun demonstrate-complement-rule (event omega)
  "The COMPLEMENT RULE is a consequence of the axioms:
       P(not A) = 1 - P(A).
   This follows because A and its complement partition Omega, so by
   additivity P(A) + P(not A) = P(Omega) = 1."
  (- 1 (probability-classical event omega)))

(defun main ()
  (let ((omega (make-two-dice-sample-space)))
    (format t "Sample space Omega: two dice, |Omega| = ~a~%" (length omega))
    (format t "Axiom check P(Omega) = ~a~%" (demonstrate-axiom-normalization omega))

    (let ((e7 (event-sum-is omega 7)))
      (format t "~%Event 'sum = 7': ~a favorable outcome(s)~%" (length e7))
      (format t "  P(sum = 7) = ~a = ~4f~%" (probability-classical e7 omega)
              (float (probability-classical e7 omega))))

    (let ((e6 (event-at-least-one-is omega 6)))
      (format t "~%Event 'at least one 6': ~a favorable outcome(s)~%" (length e6))
      (format t "  P(at least one 6) = ~a = ~4f~%" (probability-classical e6 omega)
              (float (probability-classical e6 omega)))
      (format t "  Complement P(no 6) = 1 - P(at least one 6) = ~4f~%"
              (float (demonstrate-complement-rule e6 omega))))

    (let ((e10 (event-sum-at-least omega 10)))
      (format t "~%Event 'sum >= 10': ~a favorable outcome(s)~%" (length e10))
      (format t "  P(sum >= 10) = ~a = ~4f~%" (probability-classical e10 omega)
              (float (probability-classical e10 omega))))))

(main)
