# Inbreeding and demographic history of caribou in western Canada

Code accompanying:

> Bourbon, C., Deakin, S., Michalak, A., Hughes, M. M., Cavedon, M., Neufeld, L.,
> Pelletier, A., Polfus, J., Schwantje, H., Thacker, C., Musiani, M., & Poissant, J.
> (2026). Inbreeding and demographic history of caribou (*Rangifer tarandus*) in
> western Canada inferred from genome-wide SNP data. *Evolutionary Applications*,
> 19(8), e70311. https://doi.org/10.1111/eva.70311

This repository contains the analysis code used to estimate inbreeding (runs of
homozygosity, FROH, and kinship), contemporary and historical effective
population size (Ne), and their relationships with census size (Nc) across
caribou subpopulations, pooled genetic units, and metapopulations.

## Data

Genotype data are **not** included in this repository. The SNP genotypes are
publicly archived on Dryad:

> https://doi.org/10.5061/dryad.8pk0p2p3x

The analyses start from the quality-controlled dataset produced by
`01_quality_control.R`: **33,346 autosomal SNPs for 759 individuals** across 45
subpopulations. Download the genotypes from Dryad and place them under `data/`
(see [Repository layout](#repository-layout)) before running any script.

## Pipeline overview

The scripts are numbered in the order they are meant to be run. R scripts were
run under R 4.2.3; shell scripts are SLURM batch scripts for a Linux cluster and
call external population-genetics software (paths are set at the top of each
script).

| Step | Script | Language | What it does |
|------|--------|----------|--------------|
| 1 | `01_quality_control.R` | R | Quality control and filtering from the merged genotypes to the final 759-individual / 33,346-SNP autosomal dataset (missingness, duplicates, herd assignment, HWE, sex-chromosome removal). Also computes heterozygosity, missingness, and an empirical genotyping-error rate from duplicate samples. |
| 2 | `02_roh_inbreeding.R` | R | ROH detection (detectRUNS, consecutive method), per-individual NROH/LROH/SROH, genome-wide FROH, FROH by length class, and group comparisons (Kruskal–Wallis + ε², Dunn's, one-sample Wilcoxon) at the subpopulation and metapopulation levels. |
| — | `froh_inbreeding_class.R` | R | Helper sourced by step 2: FROH partitioned into custom ROH length classes (0.3–2, 2–4, 4–6, 6–8, >8 Mb). |
| 3 | `03_kinship_comparison.R` | R | Individual kinship and the kinship-based inbreeding coefficient Fkin (popkin); comparison of Fkin with FROH via SMA regression (tested against 1:1), per-subpopulation OLS slopes, and Bland–Altman agreement. |
| 4 | `04_currentNe2.sh` | bash | Contemporary Ne with currentNe2: builds pooled genetic units (PLINK) and estimates Ne at the subpopulation, genetic-unit, and metapopulation levels, comparing the panmictic and structure-aware (`-x`) models at the metapopulation scale. Recombination rates r = 0.9/1.0/1.1. |
| 5 | `05_froh_ne_nc_models.R` | R | Relationship between FROH and population size: AICc model selection and conditional model averaging (MuMIn) of log(FROH) on log(Nc) and log(Ne), at the subpopulation and genetic-unit scales, including the percent change in FROH per 10-fold change in predictor. |
| 6 | `06_GONE2_historical.sh` | bash | Historical Ne over ~100 generations with GONE2 (panmictic model, 100 seeds, chromosome 33 removed) at the subpopulation, genetic-unit, and metapopulation levels. |
| 7 | `07_colony_Ne.sh` | bash | Structure-insensitive metapopulation Ne cross-check with COLONY (sibship reconstruction on 500 informative LD-pruned SNPs, Full Likelihood). |
| — | `sensitivity_roh.R` | R | Supporting simulation assessing ROH robustness to genotyping error at the chosen detection parameters. |

### Analyses not scripted here

Two analyses reported in the paper are not included as code:

- **Demographic-scenario testing (DIYABC-RF).** The four demographic scenarios
  (Results 3.9, Fig. S18) were run through the DIYABC-RF graphical interface, so
  there is no script to share. Parameters and priors are given in the Methods
  and Supporting Information S3.
- **Ne/Nc ratios (Tables 1 and 2).** These were computed directly from the
  currentNe2 Ne point estimates and the census sizes (Nc, Appendix Table A1) as
  simple ratios, with 90% CIs obtained by dividing the Ne CI bounds by the Nc
  point estimate. No separate script is required.

## Software requirements

**R (4.2.3)** packages: `detectRUNS`, `popkin`, `BEDMatrix`, `MuMIn`, `lmodel2`,
`FSA`, `dplyr`, `tidyr`, `readxl`, `plyr`, `data.table`.

**External tools** (called by the shell scripts; set the executable paths at the
top of each script):

- PLINK v1.9 and PLINK2 — https://www.cog-genomics.org/plink/
- currentNe2 v2.0 — Santiago et al. (2025)
- GONE2 v2.0 — Santiago et al. (2025)
- COLONY v2.0.7.2 — Jones & Wang (2010)

The shell scripts are written for a SLURM cluster. Adjust the `#SBATCH` headers
(memory, time, partition) for your environment, or run the underlying commands
directly outside SLURM.

## Repository layout

```
.
├── 01_quality_control.R
├── 02_roh_inbreeding.R
├── froh_inbreeding_class.R
├── 03_kinship_comparison.R
├── 04_currentNe2.sh
├── 05_froh_ne_nc_models.R
├── 06_GONE2_historical.sh
├── 07_colony_Ne.sh
├── sensitivity_roh.R
├── data/                     # inputs (not tracked; see Data section)
│   ├── caribou_759_autosomes.{bed,bim,fam,ped,map}   # final dataset (from step 1 / Dryad)
│   ├── subpopulations/       # one PED/MAP per subpopulation
│   ├── metapopulations/      # one PED/MAP per metapopulation
│   ├── genetic_units/        # built by steps 4 and 6 (PLINK merges)
│   ├── population_clusters.txt   # population -> metapopulation
│   ├── genetic_unit_map.txt      # subpopulation -> pooled genetic unit
│   ├── sample_sizes.txt          # genomic n per subpopulation
│   ├── census_size.xlsx          # Nc per subpopulation
│   ├── ne_subpopulations.xlsx    # currentNe2 subpopulation Ne (from step 4)
│   └── ne_genetic_units.xlsx     # currentNe2 genetic-unit Ne (from step 4)
├── output/                   # R script outputs (tables)
├── results/                  # shell script outputs (Ne estimates, logs)
└── logs/                     # SLURM logs
```

Each script defines `data_dir` / `out_dir` (R) or `$BASE` and related paths
(shell) at the top; edit these to match where you place the data. Input file
names in the scripts assume the layout above.

## Notes on reproducibility

- Marker filtering deliberately excludes MAF and LD pruning (following Meyermans
  et al. 2020), as such filtering can bias ROH detection and LD-based Ne
  estimation. The COLONY step is the exception: it uses an LD-pruned,
  high-MAF subset chosen specifically for sibship reconstruction.
- ROH detection parameters (minSNP = 20, maxOppRun = 2, maxMissRun = 2,
  minimum length 300 kb) were selected after the sensitivity assessment
  summarised in `sensitivity_roh.R` and Supporting Information S1.
- currentNe2 and GONE2 were run at r = 1.0 cM/Mb (primary) with sensitivity
  analyses at r = 0.9 and 1.1 cM/Mb. A generation time of 8 years was used to
  convert generations to calendar years.
- Steps 4 and 6 build the pooled genetic units by merging subpopulations with
  PLINK; the membership of each unit is listed in those scripts and in
  `genetic_unit_map.txt`.

## Contact

Questions about the code can be directed to the corresponding author via the
paper, or raised as an issue on this repository.
