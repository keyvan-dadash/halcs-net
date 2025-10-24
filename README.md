# HALCS-Net (Ka) — MATLAB Toolkit

Minimal MATLAB toolkit for planning-grade **lunar communications**:
- L2 **halo-relay** coverage (no SPICE)
- **Ka-band** link budgets: Relay↔Earth (RE), UE↔Earth (LE)
- **ITU-R P.618** atmosphere
- **LDPC-QPSK** Eb/N0 vs availability + margin plots
- Ground-station diversity across multiple dish sizes

---

## Requirements
- **MATLAB R2025a+** (for `p618PropagationLosses`)
- ITU digital map MAT-files  
  *(first run auto-downloads the MathWorks ITU map bundle if missing)*

---

## Quick Start
```matlab
% 1) Configure shared parameters (frequencies, sites, powers, dishes, etc.)
edit constant.m

% 2) Eb/N0 vs antenna size and transmitting power (RL)
run l2.m

% 3) Eb/N0 vs availability (RE & LE) + plot margins at p = 0.001%
run eb_n0_itu.m

% 4) L2 ring coverage map (Earth only / Relay only / Both / None)
run coverage_l2.m
