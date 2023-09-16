# Fully synthetic survey data

The code implements the simulation studies in the paper "Fully Synthetic Data for Complex Survey." The code repeatedly simulates PPS survey data, implements the proposed methods (Synrep-R and SynRep-1) and alternatives (generating pseudo-population, pseudo-simple random samples, Horvitz Thompson estimators, and two estimates ignoring the design).

The file "syn.R" considers a probability proportional to size (PPS) design treating the weights from the American Community Survey (ACS2021_w.RDS) as the measure of size. 

The file "syn_srs.R" assumes that original survey is a simple random sample.

The file "summary.Rmd" summarizes the results and generates plots used in the paper.
