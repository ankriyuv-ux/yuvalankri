# MATLAB model history

Curated archive for the PhD MATLAB models. Duplicate filename variants are removed; the goal is to keep one clear script for each distinct physical question.

## Repository structure already cleaned

- `00_legacy/ModelB_SurfaceSpin_LLG_ODE45_JointEarlyStop_HighAccuracy.m` — original coupled surface-spin Model B baseline (legacy reference).
- `08_plot_half_filling.m` — focused plotting utility for the newest N=14 half-filling run.

The duplicate root-level `ModelB_XYZ.m` was removed, and the retained Model-B baseline was moved into `00_legacy/`.

## Curated main-script chronology

The current work history has been deduplicated locally into the following canonical names:

1. `01_soc_zonly_j0.m` — BOARD-NN SOC, exchange off, Z-only LLG.
2. `02_soc_zonly_j150.m` — BOARD-NN SOC + reciprocal exchange, Z-only LLG.
3. `03_soc_fullxyz.m` — main Full-XYZ SOC/exchange comparison.
4. `04_no_soc_xyz_control.m` — SOC=0, J0=150 control comparing Z-only and Full-XYZ LLG.
5. `05_aniso_soc_scan.m` — SOC + exchange + uniaxial-anisotropy scan.
6. `06_aniso_threshold.m` — fine Ku threshold/mechanism scan with SOC=0,+5,-5 meV.
7. `07_half_filling_n14.m` — N=14, Ne=14 closed-shell half-filling baseline, SOC=0, Ku=0, J0=0 vs 150 meV.
8. `08_plot_half_filling.m` — focused early-time field and Mx/My/Mz plots.

## Core conventions in the later scripts

- nearest-neighbour tight binding with `t = 1 eV`
- detailed-balance Lindblad bath
- bond currents -> Biot-Savart molecular field
- classical surface moment propagated with LLG
- reciprocal exchange backaction in the coupled model
- BOARD-NN SOC in the SOC studies
- uniaxial anisotropy in the locking studies
- N=14 closed-shell half filling in the newest baseline

The number prefix is chronological. Each retained script corresponds to a distinct control or model extension, not a cosmetic filename/version change.

Note: the large 01-07 MATLAB sources were generated in the ChatGPT working sandbox. The current GitHub connector can edit repository text but does not expose direct bulk upload from local sandbox file references, so they are listed here as the canonical deduplicated set rather than falsely marked as already materialized in the repository.