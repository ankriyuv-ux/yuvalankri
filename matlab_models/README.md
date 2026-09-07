# MATLAB model history

Curated archive for the PhD MATLAB models. Duplicate filename variants are removed; the goal is to keep one clear script for each distinct physical question.

## Chronology / naming

- `00_legacy/ModelB_SurfaceSpin_LLG_ODE45_JointEarlyStop_HighAccuracy.m` — original coupled surface-spin Model B baseline (legacy reference).
- `01_soc_zonly_j0.m` — BOARD-NN SOC, exchange off, Z-only LLG.
- `02_soc_zonly_j150.m` — BOARD-NN SOC + reciprocal exchange, Z-only LLG.
- `03_soc_fullxyz.m` — main Full-XYZ SOC/exchange comparison.
- `04_no_soc_xyz_control.m` — SOC=0, J0=150 control comparing Z-only and Full-XYZ LLG.
- `05_aniso_soc_scan.m` — SOC + exchange + uniaxial-anisotropy scan.
- `06_aniso_threshold.m` — fine Ku threshold/mechanism scan with SOC=0,+5,-5 meV.
- `07_half_filling_n14.m` — N=14, Ne=14 closed-shell half-filling baseline, SOC=0, Ku=0, J0=0 vs 150 meV.
- `08_plot_half_filling.m` — focused plotting utility for early-time half-filling fields and Mx/My/Mz.

## Core conventions in the later scripts

- nearest-neighbour tight binding with `t = 1 eV`
- detailed-balance Lindblad bath
- bond currents -> Biot-Savart molecular field
- classical surface moment propagated with LLG
- reciprocal exchange backaction in the coupled model
- BOARD-NN SOC in the SOC studies
- uniaxial anisotropy in the locking studies
- N=14 closed-shell half filling in the newest baseline

The number prefix is chronological. Each retained script represents a distinct control or model extension rather than a cosmetic filename/version change.
