;;;; 10_markov_chains.lisp
;;;; Discrete-time Markov chains: transition matrices and stationary distribution.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; A MARKOV CHAIN is a sequence of random variables X0, X1, X2, ... taking
;; values in a STATE SPACE. Its defining property is the MARKOV PROPERTY:
;; the future depends on the present state only, not on the history:
;;
;;     P(X_{t+1} = j | X_t = i, X_{t-1}, ...) = P(X_{t+1} = j | X_t = i).
;;
;; For a finite state space, the chain is described by a TRANSITION MATRIX P,
;; where P[i][j] = P(X_{t+1} = j | X_t = i). Each ROW of P is a probability
;; distribution: P[i][j] >= 0 and Σ_j P[i][j] = 1.
;;
;; If the state distribution at time t is a row vector v_t, then
;;
;;     v_{t+1} = v_t P        (matrix multiplication).
;;
;; After n steps:  v_n = v_0 P^n.
;;
;; A STATIONARY DISTRIBUTION pi is a row vector satisfying
;;
;;     pi = pi P.
;;
;; It is a fixed point: if the chain starts distributed as pi, it stays
;; distributed as pi forever. For an IRREDUCIBLE and APERIODIC chain (one
;; that can reach any state from any other and is not locked into a cycle),
;; the stationary distribution exists, is unique, and
;;
;;     v_t -> pi   as t -> infinity,  regardless of the starting state.
;;
;; Intuitively, the chain "forgets" its start and settles into equilibrium.
;;
;; We demonstrate a simple weather model with states Sunny/Rainy and find its
;; stationary distribution by (a) iterating the chain and (b) solving the
;; linear system pi = pi P directly.
;; ============================================================================

(defun make-transition-matrix (rows)
  "Store the transition matrix as a vector of row vectors."
  (map 'vector #'(lambda (r) (coerce r 'vector)) rows))

(defun matrix-row-times-vector (row vec)
  "Dot product of a probability row and a distribution vector.
   This is one entry of v P."
  (reduce #'+ (map 'vector #'* row vec)))

(defun step-distribution (dist P)
  "Advance one step: v_{t+1} = v_t P. Returns a new distribution vector."
  (let* ((n (length dist))
         (result (make-array n :initial-element 0.0d0)))
    (dotimes (j n)
      (setf (aref result j)
            (reduce #'+ (loop for i below n
                              collect (* (aref dist i) (aref (aref P i) j))))))
    result))

(defun iterate-chain (dist P steps)
  "Compute v_steps = v_0 P^steps by repeated multiplication.
   As steps -> infinity, this converges to the stationary distribution for
   an irreducible aperiodic chain."
  (let ((d (copy-seq dist)))
    (dotimes (s steps)
      (setf d (step-distribution d P)))
    d))

(defun stationary-by-iteration (dist P steps)
  "Approximate the stationary distribution by long-run simulation of v_t P.
   This is the empirical/iterative approach: it works whenever the chain is
   ergodic (irreducible + aperiodic)."
  (iterate-chain dist P steps))

(defun stationary-by-linear-system (P)
  "Solve pi = pi P exactly for a 2-state chain.
   For 2 states with transition matrix [[a, 1-a],[b, 1-b]], the stationary
   distribution is pi = [b/(a+b), a/(a+b)] (derived from pi = pi P plus
   normalization pi1 + pi2 = 1). This closed form avoids iterative error."
  (let* ((a (aref (aref P 0) 0))   ; P[S->S] = probability of staying Sunny
         (b (aref (aref P 1) 0))   ; P[R->S] = probability Rainy -> Sunny
         ;; From pi = pi P: pi_S (1 - a) = pi_R * b, so pi_S / pi_R = b / (1-a).
         ;; With pi_S + pi_R = 1 this gives:
         ;;   pi_S = b / (1 - a + b)
         ;;   pi_R = (1 - a) / (1 - a + b)
         (denom (+ (- 1 a) b)))
    (vector (/ b denom) (/ (- 1 a) denom))))

(defun main ()
  (let* ((P (make-transition-matrix
              '((0.8 0.2)   ; from Sunny: 80% stay Sunny, 20% -> Rainy
                (0.4 0.6)))) ; from Rainy: 40% -> Sunny, 60% stay Rainy
         (states '("Sunny" "Rainy"))
         (start #(1.0d0 0.0d0)))  ; start certainly Sunny
    (format t "=== Markov Chain: Weather Model ===~%")
    (format t "States: ~a~%" states)
    (format t "Transition matrix P (rows = from-state):~%")
    (format t "  Sunny -> [~4f, ~4f]~%" (aref (aref P 0) 0) (aref (aref P 0) 1))
    (format t "  Rainy -> [~4f, ~4f]~%" (aref (aref P 1) 0) (aref (aref P 1) 1))
    (format t "~%Start distribution: [1, 0] (certainly Sunny).~%")
    (format t "Iterating v_{t+1} = v_t P:~%")
    (loop for steps in '(0 1 2 5 10 50 100) do
      (let ((d (iterate-chain start P steps)))
        (format t "  steps=~4a: [~6,4f, ~6,4f]~%" steps (aref d 0) (aref d 1))))
    (format t "~%Stationary distribution (by long iteration, t=10000):~%")
    (let ((pi-iter (stationary-by-iteration start P 10000)))
      (format t "  [~6,4f, ~6,4f]~%" (aref pi-iter 0) (aref pi-iter 1)))
    (format t "Stationary distribution (exact, solving pi = pi P):~%")
    (let ((pi-exact (stationary-by-linear-system P)))
      (format t "  [~6,4f, ~6,4f]~%" (aref pi-exact 0) (aref pi-exact 1)))
    (format t "~%The chain forgets its initial state and converges to the~%")
    (format t "unique stationary distribution (ergodic theorem for Markov chains).~%")))

(main)
