#lang racket

(require "declarations.rkt")
(provide buildTree calcForces moveparticles)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(define (without-null l)
  (remove '() (remove '() (remove '() (remove '() l)))))


(define (buildTree initialArea particles)
  (let* ((x1 (bbox-llx initialArea))
         (x2 (bbox-rux initialArea))
         (x3 (/ (+ x2  x1) 2))
         (y1 (bbox-lly initialArea))
         (y2 (bbox-ruy initialArea))
         (y3 (/ (+ y1 y2) 2)))



(define (quadrant1 particles)
    (define (check str)
       (if (and (< (vec-x (particle-posn str)) x3)
                (>= (vec-y (particle-posn str)) y3)) #t
                                                                   #f))
    (filter check particles))

(define (quadrant2 particles)
    (define (check str)
       (if (and (>= (vec-x (particle-posn str)) x3)
                (>= (vec-y (particle-posn str)) y3)) #t
                                                                   #f))
    (filter check particles))


(define (quadrant3 particles)
    (define (check str)
       (if (and (< (vec-x (particle-posn str)) x3)
                (< (vec-y (particle-posn str)) y3)) #t
                                                                   #f))
    (filter check particles))

(define (quadrant4 particles)
    (define (check str)
       (if (and (>= (vec-x (particle-posn str)) x3)
                (< (vec-y (particle-posn str)) y3)) #t
                                                                   #f))
    (filter check particles))

    
  (cond [(null? particles)'()]
        [(singleton particles) (gnode (particle-mass (car particles)) (particle-posn (car particles)) '())]
        [#t (gnode (sum (masses particles))
              (vec (vec-x (com particles))
                   (vec-y (com particles)))
              (without-null (list (buildTree (bbox x1 y3 x3 y2) (quadrant1 particles))
                    (buildTree (bbox x3 y3 x2 y2) (quadrant2 particles))
                    (buildTree (bbox x1 y1 x3 y3) (quadrant3 particles))
                    (buildTree (bbox x3 y1 x2 y3) (quadrant4 particles)))))])))


(define (masses particles)
  (if(null? particles) '()
     (cons (particle-mass (car particles)) (masses (cdr particles)))))


(define (com particles)
  (let* ((added-mass (sum (masses particles))))
  (define (com-helper particles x-coor y-coor)
    (if (null? particles) (vec (/ x-coor added-mass)
                               (/ y-coor added-mass))
        (com-helper (cdr particles) (+ x-coor (* (particle-mass (car particles)) (vec-x (particle-posn (car particles)))))
                                    (+ y-coor (* (particle-mass (car particles)) (vec-y (particle-posn (car particles))))))))
  (com-helper particles 0 0)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(define (distance p q)
  (sqrt (+ (expt (- (vec-x p) (vec-x q)) 2)
           (expt (- (vec-y p) (vec-y q)) 2))))

(define (cube d)
  (* d d d))

(define (sum* l)
  (define (xy l1 l2 l)
    (if(null? l) (vec (sum l1)
                      (sum l2))
       (xy (cons (vec-x (car l)) l1)
           (cons (vec-y (car l)) l2)
           (cdr l))))
  (xy '() '() l))
  

(define (calcForces initialArea tree particles)
  (let* ((length (abs (- (bbox-rux initialArea)
                        (bbox-llx initialArea)))))
    (define (calc-force-helper len tree p)
      (let* ((mass-particle (particle-mass p))
             (mass-gnode (gnode-mass tree))
             (x-diff (- (vec-x (gnode-posn tree)) (vec-x (particle-posn p))))
             (y-diff (- (vec-y (gnode-posn tree)) (vec-y (particle-posn p))))
             (distance-between (distance (particle-posn p)
                                         (gnode-posn tree))))

        (cond [(or (> distance-between (* theta len)) (null? (gnode-subtrees tree)))
                    (if (= distance-between 0) (vec 0 0)
                        (vec (/ (* g mass-particle mass-gnode x-diff) (cube distance-between))
                             (/ (* g mass-particle mass-gnode y-diff) (cube distance-between))))]
              [#t (sum* (map (lambda (x)(calc-force-helper (/ len 2) x p)) (gnode-subtrees tree)))])))
    (map (lambda(p) (calc-force-helper length tree p)) particles)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


(define (moveparticles particles forces)

  (define (move-helper p force)
    (let* ((mass (particle-mass p))
           (accn-x (/ (vec-x force) mass))
           (accn-y (/ (vec-y force) mass))
           (accn (sqrt (+ (* accn-x accn-x) (* accn-y accn-y))))
           (initialx (vec-x (particle-posn p)))
           (initialy (vec-y (particle-posn p)))
           (velx (vec-x (particle-velocity p)))
           (vely (vec-y (particle-velocity p))))
      (particle mass
                (vec (+ initialx (* velx timeslice) (* 0.5 accn-x timeslice timeslice))
                     (+ initialy (* vely timeslice) (* 0.5 accn-y timeslice timeslice)))
                (vec (+ velx (* accn-x timeslice))
                     (+ vely (* accn-y timeslice))))))

  (zipwith move-helper particles forces))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;