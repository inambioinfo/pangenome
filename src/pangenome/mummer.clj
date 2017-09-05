(ns pangenome.mummer
  (:use [pangenome.intervals]))

(defn generate-commands
  [file1 file2]
  (println file1)
  (println file2)
  
  [["nucmer" "-p" "out" "--mum" file1 file2]
   ["show-coords" "-cdTH" "out.delta"]])

(defrecord Nucmer-record 
  [ref_start
   ref_end
   query_start
   query_end
   ref_length
   query_length
   pctid
   coverage_ref
   coverage_query
   frm
   tags
   ref_landmark
   query_landmark])

(defn read-nucmer-file
  [filename]
  (let [rdr (clojure.java.io/reader filename)]
    (map
      (fn [x]
        (let [[ref_start ref_end 
               query_start query_end
               ref_length query_length
               pctid
               coverage_ref
               coverage_query
               frm
               tags
               ref_landmark
               query_landmark] (clojure.string/split x #"\t")]
        (->Nucmer-record
          (read-string ref_start)
          (read-string ref_end)
          (read-string query_start)
          (read-string query_end)
          (read-string ref_length)
          (read-string query_length)
          (read-string pctid)
          (read-string coverage_ref)
          (read-string coverage_query)
          frm
          tags
          ref_landmark
          query_landmark)))
      (drop 4 (line-seq rdr)))))

(defn read-nucmer-from-string
  [output]
  (map
    (fn [x]
      (let [[ref_start ref_end 
             query_start query_end
             ref_length query_length
             pctid
             coverage_ref
             coverage_query
             frm
             tags
             ref_landmark
             query_landmark] (clojure.string/split x #"\t")]
      (->Nucmer-record
        (read-string ref_start)
        (read-string ref_end)
        (read-string query_start)
        (read-string query_end)
        (read-string ref_length)
        (read-string query_length)
        (read-string pctid)
        (read-string coverage_ref)
        (read-string coverage_query)
        frm
        tags
        ref_landmark
        query_landmark)))
    (drop 4 (clojure.string/split (:out output) #"\n"))))

(defn convert-to-genomematch [x]
  (let [ref-coords   (sort [(:ref_start x) (:ref_end x)])
        query-coords (sort [(:query_start x) (:query_end x)])]
    (->GenomeMatch
      (->Interval (:ref_landmark x)   (first ref-coords)   (second ref-coords))
      (->Interval (:query_landmark x) (first query-coords) (second query-coords))
      (:pctid x)
      (:ref_landmark x)
      (:query_landmark x)
      (:ref_length x)
      (:query_length x)
      )))

   