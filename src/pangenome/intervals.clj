(ns pangenome.intervals)

(defrecord GenomeMatch
  [ref_coords 
   query_coords
   pct_id
   ref_landmark
   query_landmark
   ref_length
   query_length])

(defrecord Interval
  [id x y])

(defrecord Genome
  [id
   contigs])

(defn within?
  [i q]
  (and
    (=  (:id i) (:id q))
    (<= (:x i)  (:x q) (:y i))
    (<= (:x i)  (:y q) (:y i))))

(defn intersects?
  [i q]
  (and 
    (= (:id i) (:id q))
    (or
      (<= (:x i) (:x q) (:y i))
      (<= (:x i) (:y q) (:y i))
      (<= (:x q) (:x i) (:y q))
      (<= (:x q) (:y i) (:y q))
      )))


(defn merge-and-enlarge
  "Merges 2 intervals at their further boundaries"
  ([i] (merge-and-enlarge i i))
  ([i q]
  (when (= (:id i) (:id q))
    (->Interval 
      (:id i)
      (min (:x i) (:x q))
      (max (:y i) (:y q))))))

(defn intersection
  [i q]
  (cond
    (not (intersects? i q)) nil
    (within? i q) q
    (within? q i) i

    (<= (:x i) (:x q) (:y i))
         (->Interval (:id i) (:x q) (:y i))

    (<= (:x i) (:y q) (:y i))
         (->Interval (:id i) (:x i) (:y q))
         
    (<= (:x q) (:x i) (:y q))
         (->Interval (:id i) (:x i) (:y q))

    (<= (:x q) (:y i) (:y q))
         (->Interval (:id i) (:x q) (:y i))

    :else (throw (Exception. "Invalid intersection call"))))

(defn length?
  [q]
  (Math/abs
    (- (:y q) (:x q))))

(defn distance
  [i q]
  (if (= (:id i) (:id q))
    (min
      (if (intersects? i q) 0 Double/POSITIVE_INFINITY)
      (Math/abs (- (:x q) (:y i)))
      (Math/abs (- (:x q) (:x i)))
      (Math/abs (- (:y q) (:y i)))
      (Math/abs (- (:y q) (:x i))))
    Double/POSITIVE_INFINITY))

(defn subtract
  [i q]
  (if 
    (within? i q)
      (remove 
        (fn [x]
          (<= (length? x) 5))
        [(->Interval (:id i) (:x i) (:x q))
         (->Interval (:id i) (:y q) (:y i))])
      (throw (Exception. "Invalid subtract for intervals"))))

(defn subtract-from-intervals
  [blocks q]
  (doall
    (flatten
      (doall
        (for [block blocks]
          (cond 
            (within? block q) (subtract block q)
            (intersects? block q) (subtract block (intersection block q))
            :else block))))))

(defn split-interval 
  [i q]
  (if
    (or
      (within? i q)
      (intersects? i q))
        (let [sorted (sort 
                       [(:x i) (:x q) (:y i) (:y q)])
              x (fn [y] (nth sorted y))]
          [(->Interval (:id i) (x 0) (x 1))
           (->Interval (:id i) (inc (x 1)) (x 2))
           (->Interval (:id i) (inc (x 2)) (x 3))])))
        
(defn split-intervals
  [blocks q]
  (flatten
    (for [block blocks]
      (cond 
        (within? block q) (split-interval block q)
        (intersects? block q) (subtract block (intersection block q))
        :else block))))


(defn ->Gene-to-interval [gene]
  (->Interval 
    (:landmark gene) 
    (:start gene) 
    (:end gene)))
