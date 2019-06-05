# ontogwas <img src="man/figures/logo.png" align="right" width="100">

Ontology-based analysis of genomic variants.

## Installation

```r
source("https://bioconductor.org/biocLite.R")
biocLite("GenomicRanges")
devtools::install_github("SimonLarsen/ontogwas")
```

## Usage

### Preparing the ontology

```r
library(magrittr)
library(ontogwas)

go <- read_obo("go.obo", c("is_a","part_of")) # 'go.obo' from http://geneontology.org/page/download-ontology
anno <- read_gaf("goa_human.gaf", exclude_evidence="IEA") # 'goa_human.gaf' from http://geneontology.org/page/download-go-annotations

go <- go %>%
  filter_obsolete() %>%
  unite_roots() %>%
  annotate(anno) %>%
  filter_unannotated() %>%
  propagate_annotations() %>%
  collapse_redundant_terms(threshold=0.9)

saveRDS(go, "ontology.rds")
```

### Computing principal components

```r
library(snpStats)

go <- readRDS("ontology.rds")
snps <- read.plink("geno.bed", "geno.bim", "geno.fam")
genemap <- make_genemap(snps$map, "hg19", maxgap=10e3)
pc <- compute_term_pcs(go, snps$genotypes, genemap)
```

### Performing association tests

```r
covars <- read.table("covars.tsv")

results <- test_terms(class ~ sex+age+PC1+PC2+PC3+PC4, covars, go, pc, family=binomial("logit"))
pvalues <- sapply(results$test, "[", 2, 5)
```

### Visualizing results

```r
p.adj <- p.adjust(pvalues, "bonferroni")
sig.terms <- names(which(p.adj < 0.05))
plot_subtree(go, "tree.dot", sig.terms, pvalues=p.adj)
DiagrammeR::grViz("tree.dot")
```
