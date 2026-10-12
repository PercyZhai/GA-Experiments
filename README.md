# GA-Experiments

Code for the numerical experiments in

> Percy S. Zhai, Mladen Kolar, and Wei Biao Wu.
> **Simultaneous Inference for Covariance and Precision Matrices of Long-Range Dependent Time Series.**
> *IEEE Transactions on Information Theory*, 72(6):4246–4296, 2026.
> [doi:10.1109/TIT.2026.3685246](https://doi.org/10.1109/TIT.2026.3685246)

The repository contains code only. Simulation results, figures, tables and logs are created
by the scripts in `results/`, `output/` and `logs/`; the ABIDE data are downloaded into `data/`
(see Section 3 below). All four directories are git-ignored.

## Repository layout

```
GA-Experiments/
├── functions/                 shared helpers (sourced by the scripts below)
│   ├── bootstrap_sample.R     block bootstrap draws for covariance / precision (simulations)
│   ├── bootstrap_ci.R         entrywise block bootstrap CIs for the precision matrix (ABIDE)
│   ├── load_abide.R           loads the ABIDE NYU time series used in Sec. VIII
│   ├── aal_info.R             AAL-116 labels and plotting coordinates
│   └── graph_plotter.R        plots a conditional independence graph on the AAL layout
├── simulation/                data generation (Sec. VII, Appendix F, Appendix I)
│   ├── lowdim_sim.R           low-dimensional regime, covariance and precision
│   ├── highdim_sim.R          high-dimensional regime, covariance only
│   ├── run_grid.sh            runs the whole simulation grid (SLURM or sequential)
│   └── submit.sh              SLURM job script for one cell of the grid
├── analysis/                  figures and tables from the simulation results
│   ├── lowdim/                qq_{cov,prec}.R, ecdf_{cov,prec}.R, ks_{cov,prec}.R, w1_{cov,prec}.R
│   └── highdim/               qq_cov.R, ecdf_cov.R, ks_cov.R, w1_cov.R
└── real-data/                 ABIDE fMRI analysis (Sec. VIII)
    ├── 01_read_abide.R        reads the ROI time series into data/time.series.ABIDE.RData
    ├── 02_long_memory.R       Hurst exponents and ACF (long-range dependence)
    └── 03_cig_bootstrap.R     brain connectome by block bootstrap
```

## Requirements

- R. The simulations for the paper were run with R 3.6.2 and the figures were drawn with
  R 4.5; the code has also been checked with R 4.6.
- R packages: `foreach` and `doParallel` (high-dimensional simulations), `transport`
  (Wasserstein-1 tables), `pracma` (Hurst exponents) and `igraph` (connectome graphs).
  Everything else uses base R.
- An optimized BLAS (OpenBLAS, MKL, Accelerate, ...) is strongly recommended for the
  high-dimensional simulations.

**All commands are run from the root of the repository.**

## Reproducing the results

### 1. Simulations (Sec. VII and Appendix I)

The data follow the Gaussian linear process $X_i = \sum_{t=0}^{N-1} A_t \epsilon_{i-t}$ with
$A_t = (t+1)^{-\beta} A$, truncation length $N = n^2$, and 200 copies of $\{X_i\}_{i=1}^n$ for every
setting. Two choices of the $p \times p$ matrix $A$ are used:

| setup      | $A_{jk}$                                         | paper        |
|------------|--------------------------------------------------|--------------|
| `toeplitz` | $(\lvert j-k\rvert+1)^{-2}$                       | Sec. VII     |
| `banded`   | $(\lvert j-k\rvert+1)^{-2}\,1\{\lvert j-k\rvert \le 2\}$ | Appendix I |

One cell of the grid:

```bash
Rscript simulation/lowdim_sim.R 1000 10 0.9            # n, p, beta (setup defaults to toeplitz)
Rscript simulation/lowdim_sim.R 1000 10 0.9 banded
Rscript simulation/highdim_sim.R 400 5000 2 toeplitz
```

The whole grid (low-dimensional: $n \in \{200, 500, 1000, 2000, 4000, 10000\}$,
$p \in \{2, 10, 30, 100\}$; high-dimensional: $n \in \{200, 400, 800\}$, $p \in \{2000, 5000, 8000\}$;
$\beta \in \{2, 0.9, 0.55\}$):

```bash
bash simulation/run_grid.sh toeplitz all     # Sec. VII
bash simulation/run_grid.sh banded all       # Appendix I
bash simulation/run_grid.sh toeplitz lowdim  # only the low-dimensional grid
```

(usage: `bash simulation/run_grid.sh [toeplitz|banded] [lowdim|highdim|all]`, defaults
`toeplitz all`). `run_grid.sh` submits one SLURM job per cell through `simulation/submit.sh` when
`sbatch` is available and otherwise runs the cells one after another. The results are written to
`results/lowdim/E2quick_<n>_<p>_<beta>_toeplitz.RData`, `results/lowdim/E2Bquick_<n>_<p>_<beta>_banded.RData`,
`results/highdim/E3quick_<n>_<p>_<beta>_toeplitz.RData` and
`results/highdim/E3Bquick_<n>_<p>_<beta>_banded.RData`, and the R output of every cell to
`logs/<regime>_<n>_<p>_<beta>_<setup>.Rout`.

The largest cells are expensive: the low-dimensional cells with $n = 10000$ use FFTs of length
$10^8$, and the high-dimensional cells with $p = 8000$ hold several matrices of 10 GB (data) and
0.5 GB ($p \times p$, per parallel worker). The high-dimensional script runs its Gaussian
approximation and bootstrap loops in parallel on `SLURM_CPUS_PER_TASK` cores (all cores but one
outside SLURM). The `#SBATCH` defaults in `simulation/submit.sh` (1 CPU, 32 GB, 1 day) are used
for the low-dimensional cells; for the high-dimensional cells `run_grid.sh` requests 28 CPUs,
200 GB and 3 days, as in the published runs (`HIGHDIM_SBATCH_OPTS` in `run_grid.sh`). Adjust
these to your cluster.

### 2. Figures and tables

Once the results of a setup are available:

```bash
for f in analysis/lowdim/*.R analysis/highdim/*.R; do
  Rscript "$f" toeplitz      # use "banded" for Appendix I
done
```

Figures are written to `output/figures/` and tables (as CSV) to `output/tables/`. The file names
follow the paper's figure files: prefix `E2` (low-dimensional), `E3` (high-dimensional), with an
extra `B` for the banded setup (`E2B`, `E3B`).

### 3. ABIDE fMRI data (Sec. VIII)

1. Download the ABIDE I data preprocessed by the
   [Preprocessed Connectomes Project](http://preprocessed-connectomes-project.org/abide/)
   (pipeline `cpac`, strategy `nofilt_noglobal`, derivative `rois_aal`, site NYU) and the
   phenotypic file, and place them as

   ```
   data/abide/Phenotypic_V1_0b_preprocessed1.csv
   data/abide/rois_aal/<FILE_ID>_rois_aal.1D      e.g. NYU_0050952_rois_aal.1D
   ```

   The download links (or PCP's `download_abide_preproc.py`) are given in the header of
   `real-data/01_read_abide.R`; there are 175 NYU files. The data are not redistributed here;
   please follow the ABIDE / PCP terms of use.
2. Run the three scripts in order:

   ```bash
   Rscript real-data/01_read_abide.R      # -> data/time.series.ABIDE.RData
   Rscript real-data/02_long_memory.R     # Fig. 10, Fig. 11, ACF statistic
   Rscript real-data/03_cig_bootstrap.R   # Fig. 12, changes in the number of connections per ROI
   ```

Of the 175 NYU subjects, `01_read_abide.R` skips four that were not part of the paper's data set
(leaving 73 ASD subjects and 98 controls), and `functions/load_abide.R` removes five more whose
data contain a constant ROI series. The analysis therefore uses 71 ASD subjects and 95 controls
($n = 176$ time points, $p = 116$ ROIs).

## Paper figures and tables

| Paper | Content | Script | Output |
|---|---|---|---|
| Fig. 1–3 | QQ-plots, covariance, low-dim, $\beta = 2, 0.9, 0.55$ | `analysis/lowdim/qq_cov.R` | `E2.QQ.cov.<beta>.pdf` |
| Fig. 4–6 | QQ-plots, precision, low-dim | `analysis/lowdim/qq_prec.R` | `E2.QQ.prec.<beta>.pdf` |
| Fig. 7–9 | QQ-plots, covariance, high-dim | `analysis/highdim/qq_cov.R` | `E3.QQ.cov.<beta>.pdf` |
| Fig. 10 | Hurst exponents, ABIDE | `real-data/02_long_memory.R` | `Hurst.ASD.pdf`, `Hurst.ctrl.pdf` |
| Fig. 11 | ACF plots, ABIDE | `real-data/02_long_memory.R` | `acf.1.pdf` – `acf.4.pdf` |
| Fig. 12 | Connectome by block bootstrap, ABIDE | `real-data/03_cig_bootstrap.R` | `Autism_02_95.pdf`, `Control_02_95.pdf` |
| Fig. 13 | Connectome by FGLasso, ABIDE | not included, see Notes | |
| Fig. 14 | ECDF, covariance, low-dim | `analysis/lowdim/ecdf_cov.R` | `E2.ECDF.cov.pdf` |
| Fig. 15 | ECDF, precision, low-dim | `analysis/lowdim/ecdf_prec.R` | `E2.ECDF.prec.pdf` |
| Fig. 16 | ECDF, covariance, high-dim | `analysis/highdim/ecdf_cov.R` | `E3.ECDF.cov.pdf` |
| Fig. 17–19 | QQ-plots, covariance, low-dim, banded | `analysis/lowdim/qq_cov.R banded` | `E2B.QQ.cov.<beta>.pdf` |
| Fig. 20 | ECDF, covariance, low-dim, banded | `analysis/lowdim/ecdf_cov.R banded` | `E2B.ECDF.cov.pdf` |
| Fig. 21–23 | QQ-plots, precision, low-dim, banded | `analysis/lowdim/qq_prec.R banded` | `E2B.QQ.prec.<beta>.pdf` |
| Fig. 24 | ECDF, precision, low-dim, banded | `analysis/lowdim/ecdf_prec.R banded` | `E2B.ECDF.prec.pdf` |
| Fig. 25–27 | QQ-plots, covariance, high-dim, banded | `analysis/highdim/qq_cov.R banded` | `E3B.QQ.cov.<beta>.pdf` |
| Fig. 28 | ECDF, covariance, high-dim, banded | `analysis/highdim/ecdf_cov.R banded` | `E3B.ECDF.cov.pdf` |
| Table I, II | Kolmogorov distances, low-dim, covariance / precision | `analysis/lowdim/ks_{cov,prec}.R` | `E2.KS.{cov,prec}.matrix.csv` |
| Table III | Kolmogorov distances, high-dim | `analysis/highdim/ks_cov.R` | `E3.KS.cov.matrix.csv` |
| Table IV, V | Wasserstein-1 distances, low-dim, covariance / precision | `analysis/lowdim/w1_{cov,prec}.R` | `E2.W1.{cov,prec}.matrix.csv` |
| Table VI | Wasserstein-1 distances, high-dim | `analysis/highdim/w1_cov.R` | `E3.W1.cov.matrix.csv` |
| Table VII | AAL-116 labels | `functions/aal_info.R` (label codes) | |
| Table VIII, X | Kolmogorov distances, low-dim, banded | `analysis/lowdim/ks_{cov,prec}.R banded` | `E2B.KS.{cov,prec}.matrix.csv` |
| Table IX, XI | Wasserstein-1 distances, low-dim, banded | `analysis/lowdim/w1_{cov,prec}.R banded` | `E2B.W1.{cov,prec}.matrix.csv` |
| Table XII, XIII | Kolmogorov / Wasserstein-1 distances, high-dim, banded | `analysis/highdim/{ks,w1}_cov.R banded` | `E3B.{KS,W1}.cov.matrix.csv` |
| Appendix F | FFT simulation of the linear process | `simulation/lowdim_sim.R`, `simulation/highdim_sim.R` | |

The KS / W1 CSV files contain, for every $n$, the rows "Err vs Z" (GA section of the tables)
and "Err vs Boot" (Bootstrap section), plus an additional "Boot vs Z" row that is not reported
in the paper. Their columns are ordered by $p$, then $\beta$ (the paper groups them by $\beta$,
then $p$), and the values are not rounded.

## Notes

- **Simulation method.** Since $A_t = (t+1)^{-\beta} A$, the process is generated as $p$
  independent scalar series with the FFT algorithm of Appendix F, multiplied by $A$. The 200
  copies of $\{X_i\}_{i=1}^n$ are consecutive segments of length $n$ of one simulated series of
  length $N = n^2$, so $n \ge 200$ is required. For every copy, one block bootstrap draw
  $l^{-1/2}|\check B_{i,l}|_\infty$ is computed at a randomly chosen $i$, with $l = \lfloor n^{2/3} \rfloor$.
- **Reproducibility.** All simulations use `set.seed(42)`. Rerunning `simulation/lowdim_sim.R`
  reproduces the published low-dimensional results up to floating-point differences (about
  $10^{-12}$) caused by different R / FFT / BLAS versions. `simulation/highdim_sim.R` reproduces
  the published covariance errors $\sqrt{n}\,|\hat\Sigma_n - \Sigma|_\infty$ in the same way; its
  Gaussian approximation and bootstrap draws are computed on parallel workers whose random
  numbers are not seeded, so they differ between runs (and agree with the published ones in
  distribution). With $p = 2$ the banded and toeplitz setups coincide.
- **Fig. 13.** The comparison with the Functional Graphical Lasso (FGLasso) of
  X. Qiao, S. Guo and G. M. James, "Functional graphical models", *J. Amer. Statist. Assoc.*
  114(525):211–222, 2019, is not part of this repository.
- **Table VII** lists the AAL-116 regions; `functions/aal_info.R` contains the label codes and
  the plotting coordinates used in Fig. 12.

## License

The code is released under the [MIT License](LICENSE).

## Citation

```bibtex
@article{zhai2026simultaneous,
  title   = {Simultaneous Inference for Covariance and Precision Matrices of Long-Range Dependent Time Series},
  author  = {Zhai, Percy S. and Kolar, Mladen and Wu, Wei Biao},
  journal = {IEEE Transactions on Information Theory},
  volume  = {72},
  number  = {6},
  pages   = {4246--4296},
  year    = {2026},
  doi     = {10.1109/TIT.2026.3685246}
}
```
