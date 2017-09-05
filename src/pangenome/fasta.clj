(ns pangenome.fasta
  (:use [pangenome.intervals])
  (:require clojure.java.io
            [biotools.fasta :as fasta]))

(defn -create-contig
  [record]
  (->Interval
    (:id record)
    1 
    (count (:seq record))))

(defn load-genome
  [id file]

  (with-open [rdr (clojure.java.io/reader file)]
    (->Genome
      id
      (doall
        (map 
          -create-contig
          (fasta/parse rdr))))))
          
      