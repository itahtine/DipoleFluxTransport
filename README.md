# Dipole flux transport

This repository provides code used for the dipole flux transport (DFT) simulations in the paper

**Tähtinen, Asikainen & Mursula (2026)**  
*Ultra-fast simulations of the solar dipole and open flux*  
Submitted to Astronomy & Astrophysics

`Calc3DVectorSum.m`: Computes the 3‑component solar dipole vector (magnitude, θ, φ) from a stack of synoptic magnetograms using the vector‑sum method (vector directions from pixel locations; magnitude equals flux along axis).

`DFTSingleMap.m`: Runs DFT for a single map or a set of active regions using a precomputed propagator (`VectorMap` / `VM(t)`), producing time series of dipole magnitude and angles. Accepts either a `struct regs(i).vals/inds` or a 3‑D map stack.

`DFTSim.m`: Multi‑AR DFT simulation driver

`DFTExamplesTimeit.m`: MATLAB script used to produce the results shown in Table 1 of the paper.

`SumDipoleVectorsRad.m`: Utility that sums dipole vectors across sources/time in spherical coordinates: input is a 3‑D array `(time, [mag θ φ], N)` and it converts to/from Cartesian to return the net dipole as `(mag, θ, φ)` in radians.

`sft_sim_lin.m`: Reference SFT implementation used for validation and baseline comparisons in the paper.
