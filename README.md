# 木槌 : kizuchi

<img src='docs/tap-appear-mallet.jpeg' width='100' />

A `snakemake` workflow for building gene trees from HMM profiles.

---

`kizuchi` is a [`snakemake`](https://snakemake.readthedocs.io/en/stable/)
workflow for building phylogenetic trees, starting from a collection of genomes
and HMM profiles, producing necessary diagnostics along the way in the form of
a Jupyter notebook. The aim of this workflow is to automate and document the
tree generating procedures for phylogenetic analysis in a reproducible way.

- gene prediction using [`prodigal-gv`](https://github.com/apcamargo/prodigal-gv)
- gene annotation using [`hmmer`](http://hmmer.org/)
- ortholog scoring
- amino acid alignment using [`clustalo`](http://www.clustal.org/omega/)
- tree inference using [`fasttree`](http://www.microbesonline.org/fasttree/)

![Rule Graph](docs/rg.svg)
