(ns pangenome.order
  "Determine the order of FASTA files to process"
  (:require clojure.java.io
            [biotools.fasta :as fasta]
            [clojure.core.reducers :as r]
            [biotools.gff :as gff]
            [me.raynes.fs :as fs]))

(defn find-files [dir]
  (fs/find-files dir #".*\.fa(sta)?"))

(defn find-size-of-files [dir]
  (map
    (fn [x]
      (with-open [rdr (clojure.java.io/reader (.toString x))]
        (doall
          {:file (.toString x)
           :size (reduce +
                         (for [y (fasta/parse rdr)]
                            (count (:seq y))))})))
  (find-files dir)))

(defn find-size-of-contigs [dir]
  (map
    (fn [x]
      (with-open [rdr (clojure.java.io/reader (.toString x))]
        (doall
          (for [y (fasta/parse rdr)]
            {:file (.toString x)
             :id (:id y)
             :size (count (:seq y))}))))
  (find-files dir)))


(defn determine-order-of-files [dir]
  (map (juxt :file :size)
       (reverse (sort-by :size (find-size-of-files dir)))))


