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
   distribution is pi = [b/(1-a+b), (1-a)/(1-a+b)] (derived from pi = pi P
   plus normalization pi1 + pi2 = 1). This closed form avoids iterative error."
  (let* ((a (aref (aref P 0) 0))   ; P[S->S] = probability of staying Sunny
         (b (aref (aref P 1) 0))   ; P[R->S] = probability Rainy -> Sunny
         ;; From pi = pi P: pi_S (1 - a) = pi_R * b, so pi_S / pi_R = b / (1-a).
         ;; With pi_S + pi_R = 1 this gives:
         ;;   pi_S = b / (1 - a + b)
         ;;   pi_R = (1 - a) / (1 - a + b)
         (denom (+ (- 1 a) b)))
    (vector (/ b denom) (/ (- 1 a) denom))))

(defun solve-linear-system (matrix rhs)
  "Solve the n-by-n system A x = b by Gaussian elimination with partial
   pivoting. MATRIX is a vector of row vectors, RHS a vector; both are copied,
   not modified. Returns the solution vector x."
  (let* ((n (length rhs))
         (a (map 'vector #'copy-seq matrix))
         (x (copy-seq rhs)))
    (dotimes (col n)
      (let ((pivot col))                       ; partial pivot for stability
        (loop for r from (1+ col) below n
              when (> (abs (aref (aref a r) col)) (abs (aref (aref a pivot) col)))
                do (setf pivot r))
        (rotatef (aref a col) (aref a pivot))
        (rotatef (aref x col) (aref x pivot)))
      (let ((d (aref (aref a col) col)))
        (loop for r from (1+ col) below n
              for factor = (/ (aref (aref a r) col) d)
              do (loop for c from col below n
                       do (decf (aref (aref a r) c)
                                (* factor (aref (aref a col) c))))
                 (decf (aref x r) (* factor (aref x col))))))
    (loop for r from (1- n) downto 0           ; back substitution
          do (loop for c from (1+ r) below n
                   do (decf (aref x r) (* (aref (aref a r) c) (aref x c))))
             (setf (aref x r) (/ (aref x r) (aref (aref a r) r))))
    x))

(defun stationary-by-solve (P)
  "Stationary distribution of an n-state chain by solving the linear system
   directly, so it works for any number of states. We want pi P = pi and
   sum(pi) = 1. As a system that is (P^T - I) pi = 0, which is singular, so we
   overwrite its last equation with the normalization sum(pi) = 1 and solve."
  (let* ((n (length P))
         (a (make-array n)))
    (dotimes (i n)
      (let ((row (make-array n :initial-element 0.0d0)))
        (dotimes (j n)
          (setf (aref row j) (- (aref (aref P j) i) (if (= i j) 1.0d0 0.0d0))))
        (setf (aref a i) row)))
    (let ((b (make-array n :initial-element 0.0d0)))
      (dotimes (j n) (setf (aref (aref a (1- n)) j) 1.0d0))  ; last row = all 1s
      (setf (aref b (1- n)) 1.0d0)                            ; sum(pi) = 1
      (solve-linear-system a b))))

(defun tv-distance (u v)
  "Total variation distance between two distributions:
     ||u - v||_TV = (1/2) sum_i |u_i - v_i|.
   It is the largest gap in probability the two distributions assign to any
   event, and it is how we measure the remaining distance to stationarity."
  (* 0.5d0 (loop for i below (length u) sum (abs (- (aref u i) (aref v i))))))

(defun mixing-time (dist P pi-star epsilon)
  "Smallest number of steps t with ||v_t - pi||_TV <= EPSILON, starting from
   DIST. The distance decays geometrically at the rate set by the chain's
   second-largest eigenvalue modulus, so this t grows as epsilon shrinks."
  (let ((v (copy-seq dist)) (tt 0))
    (loop until (or (<= (tv-distance v pi-star) epsilon) (> tt 100000))
          do (setf v (step-distribution v P)) (incf tt))
    tt))

(defvar *rng-state* (make-random-state t)
  "Random state for the single-trajectory simulation.")

(defun simulate-chain (P start steps)
  "Simulate one trajectory: from state index START, take STEPS random
   transitions and return a vector counting visits to each state. Dividing the
   counts by STEPS gives the empirical state distribution, which the ergodic
   theorem says approaches the stationary distribution (an LLN for the chain)."
  (let* ((n (length P))
         (counts (make-array n :initial-element 0))
         (state start))
    (dotimes (s steps counts)
      (incf (aref counts state))
      (let ((u (random 1.0d0 *rng-state*)) (cum 0.0d0) (next (1- n)))
        (dotimes (j n)
          (incf cum (aref (aref P state) j))
          (when (<= u cum) (setf next j) (return)))
        (setf state next)))))

(defun main ()
  (let* ((P (make-transition-matrix
              '((0.8d0 0.2d0)   ; from Sunny: 80% stay Sunny, 20% -> Rainy
                (0.4d0 0.6d0)))) ; from Rainy: 40% -> Sunny, 60% stay Rainy
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
    (format t "Stationary distribution (exact, 2-state closed form):~%")
    (let ((pi-exact (stationary-by-linear-system P)))
      (format t "  [~6,4f, ~6,4f]~%" (aref pi-exact 0) (aref pi-exact 1)))
    (format t "Stationary distribution (general linear solve, any size):~%")
    (let ((pi-solve (stationary-by-solve P)))
      (format t "  [~6,4f, ~6,4f]~%" (aref pi-solve 0) (aref pi-solve 1)))

    (let ((pi-star (stationary-by-linear-system P)))
      (format t "~%Convergence to stationarity (total variation distance to pi):~%")
      (loop for steps in '(0 1 2 3 4 5 10) do
        (format t "  t=~3a  TV=~6,4f~%" steps
                (tv-distance (iterate-chain start P steps) pi-star)))
      (format t "Mixing time (TV <= 0.01) starting from [1,0]: ~a steps~%"
              (mixing-time start P pi-star 0.01d0))

      (format t "~%Simulating one 100000-step trajectory (start Sunny):~%")
      (let* ((n 100000)
             (counts (simulate-chain P 0 n)))
        (format t "  empirical P(Sunny) = ~6,4f (stationary ~6,4f)~%"
                (/ (aref counts 0) n 1.0d0) (aref pi-star 0))))

    (format t "~%The chain forgets its initial state and converges to the~%")
    (format t "unique stationary distribution (ergodic theorem for Markov chains).~%")))

(main)
