ontogwas
========

Ontology-based analysis of genomic variants.

## Installation

```r
devtools::install_github("SimonLarsen/ontogwas")
```

## Usage

### Preparing the ontology

```r
library(magrittr)
library(ontogwas)

go <- read_obo("go.obo") # http://geneontology.org/page/download-ontology
anno <- read_gaf("goa_human.gaf") # http://geneontology.org/page/download-go-annotations

go <- go %>%
  filter_obsolete() %>%
  unite_roots(c("GO:0008150","GO:0003674","GO:0005575")) %>%
  annotate(anno) %>%
  propagate_annotations() %>%
  collapse_redundant_terms()

saveRDS(go, "ontology.rds")
```

### Computing principal components

```r
library(ontogwas)
library(snpStats)

go <- readRDS("ontology.rds")

snps <- read.plink("geno.bed", "geno.bim", "geno.fam")

genemap <- make_genemap(snps$map, "hg19", maxgap=10e3)

pc <- compute_term_pcs(go, snps$genotypes, genemap, npcs=4, num_threads=8)

saveRDS(pc, "term_pcs.rds") # warning: large file
```

### Performing association tests

```r
library(ontogwas)

go <- readRDS("ontology.rds")
pc <- readRDS("term_pcs.rds")

```
