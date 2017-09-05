(ns pangenome.merge
  (:use [pangenome.intervals])
  (:require [pangenome.mummer :as mummer]
            [pangenome.order :as order]
            [clojure.java.shell :as shell]
            [biotools.fasta :as fasta]
            [pangenome.fasta :as pfasta]
            [biotools.gff :as gff]))

(def de-novo-contig-number (atom 0))

(defn find-nearby-intervals
  [max-distance intervals]
  (filter
    (fn [[a b]]
      (<=
        (distance a b)
        max-distance))
        (partition
          2
          1
          (sort-by
            :id
            (sort-by
              :x
              intervals)))))

(defn -merge-nearby-intervals
  [max-distance intervals]
  (loop [final-list []
         unprocessed (sort-by :id (sort-by :x intervals))]
    
    (cond
      (zero? (count unprocessed)) (apply conj final-list unprocessed)
      (= 1 (count unprocessed))   (apply conj final-list unprocessed)
      (<= 2 (count unprocessed))
        (let [[i q] (take 2 unprocessed)]
          (if
            (<=
              (distance i q)
              max-distance)
            (recur final-list
                   (cons (merge-and-enlarge i q) (drop 2 unprocessed)))
            (recur (conj final-list (first unprocessed))
                   (rest unprocessed)))))))

(defn merge-nearby-intervals
  [max-distance intervals]
  (loop [previous-run intervals
         current-run (-merge-nearby-intervals max-distance intervals)]
    (if (=
          previous-run
          current-run)
      current-run
      (recur current-run (-merge-nearby-intervals max-distance intervals)))))

(defn -expand-to-genes
  [intervals query-genes]
  (for [interval intervals]
    (let [intersecting-genes (filter
                               (fn [gene]
                                 (intersects?
                                   interval
                                   (->Gene-to-interval gene)))
                               query-genes)]
      (if (< 0 (count intersecting-genes))
        (do
          (reduce
            merge-and-enlarge
            interval
            (map ->Gene-to-interval intersecting-genes)))
        interval))))

(defn expand-to-genes
  [intervals query-genes]
  (loop [previous-run intervals
         current-run (-expand-to-genes intervals query-genes)]
    (if
      (= previous-run current-run)
      current-run
      (recur current-run (-expand-to-genes intervals query-genes)))))


; Handle all gene expansions and merging
(defn maximize-gene-coverage-and-merge
  [max-distance ; To merge nearby intervals
   query-genes ; Genes found in the query assembly
   intervals ; Current list of intervals
   ]
  
   (loop [previous-set intervals
          current-set (merge-nearby-intervals 
                        max-distance
                        (expand-to-genes intervals query-genes))]
          
          (if 
            (= previous-set current-set) current-set
            (recur 
              current-set
              (merge-nearby-intervals 
                max-distance
                (expand-to-genes current-set query-genes))))))

(defn read-nucmer-output [file]
  (pangenome.mummer/read-nucmer-file file))

(defn number-blocks [gm]
  (apply
    merge
    (map hash-map
         (range (count gm))
         gm)))

(defn runs-of-synteny
  [q-order]
  (loop [p (partition 2 1 q-order)
                        temp []
                        vals []]
                   (if (empty? p)
                     vals
                       (if (or
                             (= (inc (ffirst p)) (second (first p)))
                             (= (dec (ffirst p)) (second (first p))))
                         (recur (rest p)
                                (conj temp (ffirst p))
                                vals)
                         (recur (rest p)
                                []
                                (conj vals (conj temp (ffirst p))))))))

(defn merge-genomes
  [ref-file ref-gff query-file query-gff output results]
  
  (let [config {:pctid 90
                :max-dist-merge 500
                :min-size 500}
        
        ref-genome  (pfasta/load-genome "Reference" ref-file)
        query-genome (pfasta/load-genome "Query" query-file)
        
        ; Process Mummer output
        unfiltered-gm (map mummer/convert-to-genomematch results)
        
        ; Filter on PctID
        filtered-gm   (filter 
                        (fn [x]
                          (<= (:pctid config) (:pct_id x))) 
                        unfiltered-gm)
        
        ; Number blocks based on reference
        numbered-blocks  (number-blocks filtered-gm)
        
        ; Get the reference order, and query-order
        r-order (keys (sort-by (comp :x :ref_coords   val) numbered-blocks))
        q-order (keys (sort-by (comp :x :query_coords val) numbered-blocks))
        
        ; Count runs of synteny, not needed but should be exported
        ; synteny (runs-of-synteny q-order)
       
        novel-sequence (reduce 
                         subtract-from-intervals 
                         (:contigs query-genome)
                         (reverse
                           (sort-by
                             length?
                             (map 
                               (comp :query_coords val) 
                               numbered-blocks))))
        
        query-genes (with-open [rdr (clojure.java.io/reader query-gff)]
                      (doall
                        (remove 
                          (fn [x] (= (:type x) "databank_entry"))
                          (gff/parse-reader rdr))))
        
        novel-sequence-processed (maximize-gene-coverage-and-merge 500 query-genes novel-sequence)
        
        ; Remove novel sequences less than 500bp
        novel-sequence-filtered (remove 
                                  (fn [x] (< (length? x) (:min-size config)))
                                  novel-sequence-processed)

        novel-genes (filter
                      identity
                      (for [gi novel-sequence-filtered
                            gene query-genes]
                        (do
                          (if 
                            (intersects? 
                              gi
                              (->Interval (:landmark gene) (:start gene) (:end gene)))
                            [gi gene]))))
        
        novel-intervals (apply merge 
                               (map hash-map 
                                    (range
                                      (inc @de-novo-contig-number)
                                      (+
                                        @de-novo-contig-number
                                        (inc (count novel-sequence-filtered)))) 
                                    novel-sequence-filtered))
        
        novel-intervals-invert (clojure.set/map-invert (into {} novel-intervals))
        
        ; TODO List:
        ; Expand intervals to encompass gene boundaries...
        ; Merge intervals within 500bp or so...
        ; Order "pr" contigs by size (largest to smallest)
        ; Identify when a pr contig is wholly contained inside of another (newer) contig
        ]
    
    (def novel-intervals novel-intervals)
    (def novel-sequence-filtered novel-sequence-filtered)
    (def query-genes query-genes)
    
    (println "Novel genes:" (count novel-genes))
    (println "Novel Intervals:" (count novel-intervals))
    ;(println (take 5 novel-sequence-filtered))
    ;(println (take 5 novel-sequence-expanded))
    ;(println (take 5 novel-intervals))
        
    (reset! de-novo-contig-number (+ (count novel-sequence-filtered) 1 1 @de-novo-contig-number))
    
    ; Create the new reference sequence file
    
    (shell/sh "cp" ref-gff (str output ".gff3"))
    
    (with-open [
                panref-wrtr (clojure.java.io/writer  (str output ".fasta"))
                panref-gff-wrtr (clojure.java.io/writer (str output ".gff3") :append true)
                ref-rdr     (clojure.java.io/reader ref-file)
                query-rdr   (clojure.java.io/reader query-file)
                ]
      
      ; Output the entire reference sequence in the new file
      (doseq [r (fasta/parse ref-rdr)]
        (.write panref-wrtr (str ">" (:id r) "\n"))
        (.write panref-wrtr (apply str (:seq r)))
        (.write panref-wrtr "\n"))
      
      ; Output only the appropriate parts of the query sequence into the new file
      (doseq [r (fasta/parse query-rdr)]
        (let [novel-contigs (apply hash-set (map :id novel-sequence-filtered))]
          (when (novel-contigs (:id r))
            (doseq [[b i] novel-intervals
                    :when (= (:id i) (:id r))]
              (try
                (.write panref-wrtr (str ">pr." b "\n"))
                
                ; Gotta figure this out
                ; Getting null pointer exception from the following code
                ; Sometimes (:y i) is longer than the length of the sequence(!)
                ; Trying to shorten it to the length of the sequence but getting an unexplained null pointer exception
                ; from this let/if statement
                (let [final-y (if
                                (>= (:count (:seq r)) (:y i))
                                (:y i)
                                (count (:seq r)))]
                  (.write panref-wrtr 
                    (subs 
                      (apply 
                        str 
                        (:seq r)) 
                      (dec (:x i))
                      final-y)))
                (.write panref-wrtr "\n")
                (catch Exception e
                  (println
                    (str "caught exception: " (.getMessage e)
                         " "
                         e
                         " "
                         b
                         " More Data: "
                         (str (:id r) " " (count (:seq r)))
                         " "
                         (:id i)
                         " "
                         (:x i)
                         " "
                         (dec (:x i))
                         " "
                         (:y i)))
                       
                       ))))))
      
      ; Output the gene definitions into a GFF file
      (doseq [[i gene] novel-genes]
        (let [landmark  (str "pr." (get novel-intervals-invert i))
              new-start (- (inc (:start gene)) (:x i))
              new-end   (- (:end gene)   (:x i))]
          
          (when (<= new-start 0)
            (println new-start)
            (println gene)
            (println i))
          
          (.write
            panref-gff-wrtr
            (clojure.string/join
              "\t"
              [landmark
               (:source gene)
               (:type gene)
               new-start
               new-end
               (:score gene)
               (:strand gene)
               (:phase gene)
               (clojure.string/join 
                 ";"
                 (map 
                   (fn [x]
                     ; Need to do URI Encoding here...
                     (str x "=" (get gene x)))
                   (remove 
                     #{:landmark :strand :start :end :phase :type :source :score} 
                     (keys gene))))
               ]))
          
          (.write panref-gff-wrtr "\n")
          )
          
        )
      )
    )
  )


