# Dipole flux transport

This repository provides code used for the dipole flux transport (DFT) simulations in the paper

**Tähtinen, Asikainen & Mursula (2026)**  
*Ultra-fast simulations of the solar dipole and open flux*  
Submitted to Astronomy & Astrophysics

`CreateFullPropagator.m`: Run this to create full propagator matrices needed for DFT simulations

`DFTExamplesTimeit.m`: Produces the results shown in Table 1 of the paper.

`DFTSingleMap.m`: Runs DFT for a single map or a set of active regions using a precomputed propagator matrix, producing time series of dipole vectors.

`DFTSim.m`: Runs DFT for multiple active regions emerging at arbitrary times

`Calc3DVectorSum.m`: Computes the 3‑component solar dipole vector from a stack of synoptic magnetograms using the vector sum method of Tähtinen et al. (2024,2026)

`SumDipoleVectorsRad.m`: Utility that sums dipole vectors.

`sft_sim_lin.m`: Reference SFT implementation used for validation and baseline comparisons in the paper.
