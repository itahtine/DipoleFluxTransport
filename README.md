# Dipole flux transport

Companion repository for the paper:

**Tähtinen, Asikainen & Mursula (2026)**  
*Ultra-fast simulations of the solar dipole and open flux*  
Submitted to Astronomy & Astrophysics

---

## Overview

This repository provides MATLAB code used in the development and testing of the **Dipole Flux Transport (DFT)** method introduced in Tähtinen et al. (2026).  
DFT combines:

- the **Surface Flux Transport (SFT)** model  
- the **vector sum dipole representation**

to simulate the solar dipole and open solar flux (OSF) with considerable speed.

Compared to classical SFT simulations, DFT achieves speedups from:

- **100×** (daily resolution)  
- up to **50,000×** (individual active region evolution)

while reproducing Solar Cycle 24 dipole evolution within **≈1% accuracy** (see paper).

---

## Repository Contents

```
.
├── Calc3DVectorSum.m                                         # Vector-sum dipole from maps
├── DFTSim.m                                                  # DFT for multiple regions
├── DFTSingleMap.m                                            # DFT for a single map
├── SumDipoleVectorsRad.m                                     # Auxilar file for summing dipole vectors
├── DFTExamplesTimeit.m                                       # Runtime tests
├── sft_sim_lin.m                                             # SFT reference model
└── README.md
```

---

## Main Codes

### `Calc3DVectorSum.m`
Computes the 3‑component solar dipole vector (magnitude, θ, φ) from a stack of synoptic magnetograms using the vector‑sum method (vector directions from pixel locations; magnitude equals flux along axis).  
**Output:** dipole magnitude and orientation (degrees).

### `DFTSingleMap.m`
Runs DFT for a **single map or a set of active regions** using a precomputed propagator (`VectorMap` / `VM(t)`), producing time series of dipole magnitude and angles. Accepts either a `struct regs(i).vals/inds` or a 3‑D map stack.

### `DFTSim.m`
**Multi‑AR DFT simulation driver**. Given:
- `regs` (active-region pixels/indices),
- `tvec` (emergence times),
- `map0` (initial synoptic map),
- `G` (propagator matrix function/array),
- total time `T` and `timestep`,

it outputs `v, t, p` (dipole magnitude and angles) and assembles the per‑source vector matrix `VMat`. Handles edge cases at start times and aligns daily/Carrington steps.

### `DFTExamplesTimeit.m`
Example MATLAB script used in the paper to demonstrate:
- Single AR evolution (4 years)
- Solar Cycle 24 dipole evolution
- DFT vs SFT runtime comparisons
- Daily vs Carrington step effects
- 180‑day ensemble hindcast timing

These examples reproduce performance results shown in Table 1 of the paper.

### `SumDipoleVectorsRad.m`
Utility that **sums dipole vectors across sources/time** in spherical coordinates: input is a 3‑D array `(time, [mag θ φ], N)` and it converts to/from Cartesian to return the net dipole as `(mag, θ, φ)` in radians.

### `sft_sim_lin.m`
Reference **SFT** implementation used for validation and baseline comparisons in the paper.

---

## Propagator Matrices (DFT)

The core of DFT is the compact propagator matrix:

**Download links (insert later):**
- #
- #
- #

---

## Required External Data

To fully reproduce Solar Cycle 24 experiments, you need:

- **HMI synoptic magnetograms**: `hmi.synoptic_mr_polfil_720s`  
- **Active region catalog derived from SHARPs (Yeates 2020 code)`

Because redistribution is restricted, only scripts / links will be provided.

**Download links (insert later):**
- #
- #

---

## Example Usage

```matlab
% Run timing tests and comparisons (DFT vs SFT)
run DFTExamplesTimeit.m

% Compute dipole vector from a magnetogram cube
[v, theta, phi] = Calc3DVectorSum(maps);   % outputs in degrees

% DFT on a single map or AR set with a propagator matrix
[vsum, thetas, phis] = DFTSingleMap(regs, VectorMap, simLength);  % angles in degrees

% Multi-AR DFT with explicit timing and initial map
[v, t, p, VMat] = DFTSim(regs, tvec, map0, G, T, timestep);       % t,p returned in radians inside VMat
```

---

## Citation

If you use this code or reproduce results from the paper, please cite:

```
Tähtinen, I., Asikainen, T., & Mursula, K. (2026).
Ultra-fast simulations of the solar dipole and open flux.
Astronomy & Astrophysics, in press.
```

---

## Contact

Ismo Tähtinen  
Space Physics and Astronomy Research Unit  
University of Oulu, Finland  
Email: ismo.tahtinen@oulu.fi

---

# MIT License

Copyright (c) 2026 Ismo Tähtinen

Permission is hereby granted, free of charge, to any person obtaining a copy  
of this software and associated documentation files (the “Software”), to deal  
in the Software without restriction, including without limitation the rights  
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell  
copies of the Software, and to permit persons to whom the Software is  
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be  
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND,  
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES  
OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND  
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT  
HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,  
WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING  
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR  
OTHER DEALINGS IN THE SOFTWARE.
