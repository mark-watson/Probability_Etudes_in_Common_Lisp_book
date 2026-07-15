;;;; 04_binomial_geometric.lisp
;;;; Two important discrete distributions: Binomial and Geometric.

;; ============================================================================
;; PROBABILITY THEORY BACKGROUND
;; ============================================================================
;;
;; BERNOULLI TRIAL: a single experiment with two outcomes, "success" (prob p)
;; and "failure" (prob 1-p). Let X = 1 on success, 0 on failure; then
;; X ~ Bernoulli(p) with E[X] = p and Var(X) = p(1-p).
;;
;; BINOMIAL DISTRIBUTION: if we run n INDEPENDENT Bernoulli(p) trials and
;; count the number of successes, the count X has the Binomial(n, p)
;; distribution. Its PMF counts the number of sequences with exactly k
;; successes:
;;
;;     P(X = k) = C(n,k) p^k (1-p)^(n-k),   k = 0,1,...,n
;;
;; where C(n,k) = n!/(k!(n-k)!) is the binomial coefficient (the number of
;; ways to choose which k of the n trials are the successes).
;;
;;   E[X] = n p,     Var(X) = n p (1-p).
;;
;; Intuition for the mean: each trial contributes p to the expected count on
;; average, and there are n of them (linearity of expectation).
;;
;; GEOMETRIC DISTRIBUTION: run independent Bernoulli(p) trials until the FIRST
;; success occurs. Let Y = the number of trials needed (including the
;; success). Then Y has the Geometric(p) distribution:
;;
;;     P(Y = k) = (1-p)^(k-1) p,   k = 1,2,3,...
;;
;; (k-1 failures followed by one success). This is a valid PMF because it
;; sums to a geometric series: Σ_{k>=1} (1-p)^(k-1) p = p * 1/(1-(1-p)) = 1.
;;
;;   E[Y] = 1/p,     Var(Y) = (1-p)/p^2.
;;
;; The mean 1/p makes sense: if success probability is p, you wait about 1/p
;; trials on average. The MEMORYLESS PROPERTY says
;;   P(Y > m+n | Y > m) = P(Y > n):
;; past waiting does not change future chances.
;;
;; ============================================================================

(defun binomial-coefficient (n k)
  "C(n,k) = n!/(k!(n-k)!). Computed iteratively to avoid huge intermediate
   factorials. This counts the number of ways to choose k items from n."
  (if (or (< k 0) (> k n))
      0
      (loop with result = 1
            for i from 1 to k
            do (setf result (* result (/ (+ n (- k) i) i)))
            finally (return (round result)))))

(defun binomial-pmf (n p k)
  "P(X = k) for X ~ Binomial(n, p):
       C(n,k) p^k (1-p)^(n-k).
   The binomial coefficient counts WHICH trials succeed; p^k (1-p)^(n-k) is
   the probability of any ONE specific sequence with k successes."
  (* (binomial-coefficient n k)
     (expt p k)
     (expt (- 1 p) (- n k))))

(defun binomial-cdf (n p k)
  "F(k) = P(X <= k) = Σ_{i=0}^{k} P(X = i) for X ~ Binomial(n,p)."
  (loop for i from 0 to k sum (binomial-pmf n p i)))

(defun geometric-pmf (p k)
  "P(Y = k) for Y ~ Geometric(p): (1-p)^(k-1) * p.
   Reads as '(k-1) failures, then a success'."
  (* (expt (- 1 p) (- k 1)) p))

(defun geometric-tail (p k)
  "P(Y > k) = (1-p)^k for Y ~ Geometric(p).
   Surviving past k means the first k trials were ALL failures."
  (expt (- 1 p) k))

(defun demonstrate-memoryless-property (p m n)
  "Show P(Y > m+n | Y > m) = P(Y > n) for a geometric random variable Y.
   By definition of conditional probability:
     P(Y > m+n | Y > m) = P(Y > m+n) / P(Y > m) = (1-p)^(m+n) / (1-p)^m
                        = (1-p)^n = P(Y > n)."
  (let ((conditional (/ (geometric-tail p (+ m n)) (geometric-tail p m))))
    (format t "  Memoryless property check (p=~a, m=~a, n=~a):~%" p m n)
    (format t "    P(Y > m+n | Y > m) = ~a = ~4f~%" conditional (float conditional))
    (format t "    P(Y > n)          = ~a = ~4f~%" (geometric-tail p n)
            (float (geometric-tail p n)))))

(defun main ()
  (format t "=== Binomial Distribution: n=10, p=0.3 ===~%")
  (let ((n 10) (p 3/10))
    (format t "  E[X] = n p = ~a,  Var(X) = n p (1-p) = ~a~%"
            (* n p) (* n p (- 1 p))))
  (let ((n 10) (p 3/10))
    (loop for k from 0 to n do
      (format t "  P(X=~a) = ~a = ~4f   CDF F(~a) = ~a~%"
              k (binomial-pmf n p k) (float (binomial-pmf n p k))
              k (float (binomial-cdf n p k)))))
  (format t "~%=== Geometric Distribution: p=0.2 ===~%")
  (let ((p 1/5))
    (format t "  E[Y] = 1/p = ~a,  Var(Y) = (1-p)/p^2 = ~a~%"
            (/ 1 p) (/ (- 1 p) (* p p))))
  (let ((p 1/5))
    (loop for k from 1 to 10 do
      (format t "  P(Y=~a) = ~a = ~4f   P(Y>~a) = ~a = ~4f~%"
              k (geometric-pmf p k) (float (geometric-pmf p k))
              k (geometric-tail p k) (float (geometric-tail p k)))))
  (format t "~%=== Memoryless Property ===~%")
  (demonstrate-memoryless-property 1/5 3 2))

(main)
