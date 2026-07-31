;;;; 02_conditional_probability.lisp
;;;; Conditional probability, independence, the law of total probability,
;;;; and Bayes' theorem.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; CONDITIONAL PROBABILITY measures how likely an event A is GIVEN that
;; another event B is known to have occurred. We restrict our attention to
;; the world where B happened:
;;
;;     P(A | B) = P(A ∩ B) / P(B),    provided P(B) > 0.
;;
;; Intuition: once we know B occurred, B becomes our new sample space, and we
;; measure what fraction of B also contains A.
;;
;; INDEPENDENCE: Two events A and B are independent if knowledge of one gives
;; no information about the other. Formally:
;;
;;     A ⊥ B  iff  P(A ∩ B) = P(A) * P(B).
;;
;; Equivalently, P(A | B) = P(A) when A and B are independent.
;;
;; THE LAW OF TOTAL PROBABILITY: If B1, B2, ..., Bn partition the sample
;; space (they are disjoint and their union is Omega), then for any event A:
;;
;;     P(A) = Σ_i P(A | Bi) * P(Bi).
;;
;; This decomposes an hard-to-compute probability into conditional pieces.
;;
;; BAYES' THEOREM flips the direction of conditioning. It tells us how to
;; update a prior belief P(Bi) after observing evidence A:
;;
;;     P(Bi | A) = P(A | Bi) * P(Bi) / P(A)
;;
;; where P(A) is given by the law of total probability. P(Bi) is the PRIOR,
;; P(A | Bi) is the LIKELIHOOD, and P(Bi | A) is the POSTERIOR.
;;
;; ============================================================================

(defun make-two-dice-sample-space ()
  "Omega for rolling two distinguishable dice: 36 ordered pairs."
  (loop for i from 1 to 6
        append (loop for j from 1 to 6 collect (list i j))))

(defun subset-event (omega predicate)
  "Return the event { outcomes in Omega : predicate(outcome) is true }."
  (remove-if-not predicate omega))

(defun probability-classical (event omega)
  "P(A) = |A| / |Omega| for the equally-likely model."
  (/ (length event) (length omega)))

(defun conditional-probability (a b omega)
  "P(A | B) = P(A ∩ B) / P(B).
   The intersection A ∩ B is the set of outcomes in BOTH events."
  (let ((intersection (intersection a b :test #'equal)))
    (/ (length intersection) (length b))))

(defun independent-p (a b omega)
  "Are A and B independent?  Check whether P(A ∩ B) = P(A) * P(B).
   We compare exact rationals to avoid floating-point ambiguity."
  (let ((intersection (intersection a b :test #'equal)))
    (= (probability-classical intersection omega)
       (* (probability-classical a omega)
          (probability-classical b omega)))))

;; ----------------------------------------------------------------------------
;; The famous medical-testing example for Bayes' theorem.
;; ----------------------------------------------------------------------------

(defun bayes-medical-test ()
  "Apply Bayes' theorem to a classic disease-screening problem.

   Setup:
     - Disease prevalence (prior P(D))      = 1/100   (1% of population).
     - Test sensitivity  P(+ | D)           = 99/100  (true-positive rate).
     - Test specificity  P(- | no D)        = 95/100  (true-negative rate),
       so the false-positive rate P(+ | no D) = 1 - 95/100 = 5/100.

   Question: given a POSITIVE test, what is P(D | +)?
   Many people intuit ~99%, but the answer is far lower because the disease
   is rare and false positives from the large healthy population dominate."
  (let* ((p-disease 1/100)            ; prior P(D)
         (p-healthy 99/100)           ; P(no D) = 1 - P(D)
         (p-pos-given-disease 99/100) ; likelihood P(+ | D)
         (p-pos-given-healthy 5/100)) ; likelihood P(+ | no D)
    ;; Total probability of a positive test:
    ;;   P(+) = P(+ | D) P(D) + P(+ | no D) P(no D)
    (let ((p-positive (+ (* p-pos-given-disease p-disease)
                         (* p-pos-given-healthy p-healthy))))
      ;; Bayes: P(D | +) = P(+ | D) P(D) / P(+)
      (let ((posterior (/ (* p-pos-given-disease p-disease) p-positive)))
        (format t "Medical screening example (Bayes' theorem):~%")
        (format t "  Prior P(D)           = ~a = ~4f~%" p-disease (float p-disease))
        (format t "  Sensitivity P(+|D)   = ~a = ~4f~%" p-pos-given-disease
                (float p-pos-given-disease))
        (format t "  False-positive P(+|no D) = ~a = ~4f~%"
                p-pos-given-healthy (float p-pos-given-healthy))
        (format t "  Total P(+)           = ~a = ~4f~%" p-positive (float p-positive))
        (format t "  Posterior P(D | +)   = ~a = ~4f~%" posterior (float posterior))
        (format t "  => A positive test means only ~,1f% chance of disease!~%"
                (* 100 (float posterior)))
        posterior))))

(defun main ()
  (let ((omega (make-two-dice-sample-space)))
    (format t "=== Conditional Probability with Two Dice ===~%~%")
    (let* ((b (subset-event omega (lambda (o) (= (first o) 1))))   ; first die = 1
           (a (subset-event omega (lambda (o) (= (reduce #'+ o) 7))))) ; sum = 7
      (format t "B = 'first die is 1',  |B| = ~a~%" (length b))
      (format t "A = 'sum is 7',        |A| = ~a~%" (length a))
      (format t "P(A)         = ~a = ~4f~%" (probability-classical a omega)
              (float (probability-classical a omega)))
      (format t "P(A | B)     = ~a = ~4f~%" (conditional-probability a b omega)
              (float (conditional-probability a b omega)))
      (format t "Independent? ~a~%" (if (independent-p a b omega) "yes" "no")))
    (format t "~%")
    ;; A pair that IS independent: first die = 1 vs second die = 1.
    (let ((e1 (subset-event omega (lambda (o) (= (first o) 1))))
          (e2 (subset-event omega (lambda (o) (= (second o) 1)))))
      (format t "e1 = 'first die is 1', e2 = 'second die is 1'~%")
      (format t "P(e1) = ~a, P(e2) = ~a, P(e1 ∩ e2) = ~a~%"
              (probability-classical e1 omega)
              (probability-classical e2 omega)
              (probability-classical (intersection e1 e2 :test #'equal) omega))
      (format t "Independent? ~a (rolling separate dice never affects each other)~%"
              (if (independent-p e1 e2 omega) "yes" "no")))
    (format t "~%=== Bayes' Theorem ===~%~%")
    (bayes-medical-test)))

(main)
