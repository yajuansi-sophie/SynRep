# Fully synthetic survey data

The code is used for the simulation and application studies in the paper "Fully Synthetic Data for Complex Survey." The code implements the proposed methods (Synrep-R and SynRep-1) and alternatives (generating pseudo-population, pseudo-simple random samples, Horvitz Thompson estimators, psudo-likehood methods and two estimates ignoring the design).

The file "pps-synrep.R" considers a probability proportional to size (PPS) design treating the weights from the American Community Survey (ACS2021_w.RDS) as the measure of size. The code implements the simulation studies with repeated PPS sampling.

The file "srs-synrep.R" considers a simple random sample (SRS). The code implements the simulation studies with repeated SRSs.

The file "ACS_example.R" applies the methods to the ACS as an illustrative application study.



