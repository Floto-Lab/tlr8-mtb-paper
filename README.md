# CRISPR screen analysis for TLR8 paper
Codes for the Masserumule and Passemar et al. paper.
# Instructions

## Install CB2 and other dependencies
```install.packages("CB2")```

```install.packages("tidyverse")```

## Guide sequences
First download guide annotation from  the Moffat lab website http://tko.ccbr.utoronto.ca/ for the base and the supplementary library. Then run `annot2fasta.py` for both the libraries to extract FASTA files.

## Outputs
Start R and run `source(run_cb2.R)`. This will produce a number of descriptive plots.
`guides_w_pval_fc_min.txt` will have a fold change estimate per guide and an associated p-value based on the permutation test. (The results may vary across runs because of the stochastic nature of the permutation test employed.)

## Citation
Maserumule, Charlotte, Charlotte Passemar, et al. 2025. “Phagosomal RNA Sensing through TLR8 Controls Susceptibility to Tuberculosis.” Cell Reports (in press).
