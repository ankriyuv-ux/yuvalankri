# MATLAB model history

Curated, deduplicated set of the main MATLAB scripts used in the PhD model development. The numbering is chronological and each file represents a distinct physical question or model extension.

## Main scripts

1. `01_soc_zonly_j0.m` — BOARD-NN SOC, exchange off (`J0=0`), Z-only LLG baseline.
2. `02_soc_zonly_j150.m` — BOARD-NN SOC + exchange (`J0=150 meV`), Z-only LLG.
3. `03_soc_fullxyz.m` — Full-XYZ LLG with SOC/exchange; main XYZ extension.
4. `04_no_soc_xyz_control.m` — no-SOC exchange control for the Full-XYZ dynamics.
5. `05_aniso_soc_scan.m` — SOC + exchange + uniaxial-anisotropy scan.
6. `06_aniso_threshold.m` — fine `Ku` scan / finite-time locking crossover for SOC = 0,+5,-5 meV.
7. `07_half_filling_n14.m` — current N=14, Ne=14 closed-shell half-filled baseline; SOC=0, Ku=0, compare J0=0 vs 150 meV.
8. `08_plot_half_filling.m` — focused plotting/zoom script for the N=14 half-filling results.

## Later-model conventions

- nearest-neighbour tight binding with `t = 1 eV`
- open-system electronic dynamics with the Lindblad bath
- bond currents -> Biot-Savart molecular field
- surface moment propagated with LLG
- reciprocal exchange backaction where enabled
- BOARD-NN SOC in SOC studies
- uniaxial anisotropy in the locking studies
- N=14 closed-shell half filling in the newest baseline

Duplicate filename variants and the older long-name Model-B copies were removed from the working tree. Git history still preserves prior commits if an old version ever needs to be recovered.
