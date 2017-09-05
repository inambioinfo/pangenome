(ns pangenome.core
  (:use [pangenome.intervals])
  (:require clojure.java.io
            [biotools.fasta :as fasta]
            [biotools.gff :as gff]
            [clojure.core.reducers :as r]
            [pangenome.mummer :as mummer]
            [pangenome.order :as order]
            [pangenome.merge :as merge]
            [clojure.java.shell :as shell]
            [biotools.gff :as gff]))

; Still TODO:
; Sometimes there are "databank_entry" entries that span an entire interval -- Ignore them
; Sometimes gene start/ends end up negative

; Changelog
; v0.01 -- Initial
; v0.02 -- Merge nearby intervals (within 500bp of intervening sequence)
;       -- Expand de novo sequence intervals to gene boundaries when intersecting
; v0.03 -- Fixed issue of genes with negative coordinates in final gff3
; v0.04 -- Fixed issue with multiple contigs
;
; [Planned Changes]
; v0.05 -- TBD
;       -- If smaller contigs from previous iterations are found 100% in newer intervals,
;          delete previous interval (and genes?) and only incorporate new interval
;          rather than splicing up new interval into 2 additional pieces
;       -- Print out config to a stats or log file
;
;

; (process-pangenome "Smel" "../pangenome-files-mainchr")

; Need to print out config used for panref creation in stats
; or something...
(defn process-pangenome
  [name dir]

  ; Reset the contig number (mostly relevant for development)
  (reset! merge/de-novo-contig-number 0)

  (let [order (pangenome.order/determine-order-of-files dir)
        ; Kick start the pangenome with the first 2 files
        initial-files (map first (take 2 order))
        continuation (map first (drop 2 order))
        ]

    (loop [ref (first initial-files)
           query (second initial-files)
           pipeline continuation
           iteration 1]

      (let [commands (pangenome.mummer/generate-commands ref query)]
        (apply shell/sh (first commands))
        (let [results (mummer/read-nucmer-from-string (apply shell/sh (second commands)))]
          (merge/merge-genomes
            ref
            (clojure.string/replace ref #"fasta" "gff3")
            query
            (clojure.string/replace query #"fasta" "gff3")
            (str "pr" iteration)
            results)))

      (if (not (zero? (count pipeline)))
        (recur (str "pr" iteration ".fasta")
               (first pipeline)
               (rest pipeline)
               (inc iteration))))))

