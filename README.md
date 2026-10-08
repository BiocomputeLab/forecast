# forecast

Simulate and analyse Flow-seq / massively parallel reporter assay (MPRA) experiments.
A Rust port of the core of FORECAST (Gilliot & Gorochowski, *Methods in Molecular Biology* 2553, https://doi.org/10.1007/978-1-0716-2617-7_3):

- **`simulate`** generates synthetic Flow-seq data (cell sorting, PCR, sequencing).
- **`infer`** estimates each construct's expression level and variation (gamma or lognormal
  fluorescence distribution) from Flow-seq data.

## Install

```bash
cargo build --release          # binary: target/release/forecast
cargo test --release           # optional
```

## Quick start

```bash
forecast simulate --reads 1e8 --bins 8 --out-path out/sim
forecast infer --data out/sim --f-max 1e5 --count all --out-path out/inf
```

`simulate` defaults to 12 bins up to `--f-max 1e5`; `infer` takes the bin count from the data and
defaults to `--f-max 1e5`. If you simulate with a different `--f-max`, `--upper-bounds` or
`--bin-ranges`, pass the same layout to `infer`.

Progress bars are shown on a terminal (`--quiet` hides them). Errors print `error: ...` and exit
with code 2. `forecast <command> --help` lists every option.

## Defining bins

Use one of three options (the same ones work in `simulate` and `infer`):

| Option | Meaning | Example |
| --- | --- | --- |
| `--bins N --f-max X` | `N` log-spaced bins from 0 to `X`; the last bin is open-ended (`[.., inf)`) | `--bins 12 --f-max 1e5` |
| `--upper-bounds A,B,C` | Contiguous bins `0:A, A:B, B:C, C:inf` (one more bin than values). Must increase; first must be > 1 | `--upper-bounds 10,100,1000` |
| `--bin-ranges LO:HI,...` | An explicit range for every bin | `--bin-ranges 0:10,10:100,100:1000` |

`--bin-ranges` is the most flexible:

- Bins must be in increasing order and must not overlap. Gaps are allowed (cells falling in a
  gap are discarded in simulation).
- `inf` is allowed for the last upper bound only. A finite last bound models a hard sorter
  ceiling: brighter cells are discarded.
- The number of bins must equal the number of columns in the data (`infer`).

## `simulate`

```bash
# default built-in gamma library, 8 automatic bins, 100M reads
forecast simulate --reads 1e8 --bins 8

# custom bins, reproducible, 4 threads
forecast simulate --upper-bounds 10,100,1000,10000 --seed 42 --workers 4

# explicit ranges with a finite top bin
forecast simulate --bin-ranges 0:5,5:20,20:100,100:1000,1000:20000

# your own library, lognormal model
forecast simulate --distribution lognormal --library my_library.csv --bins 6 --f-max 1e4
```

| Option | Default | Description |
| --- | --- | --- |
| `--distribution` | `gamma` | `gamma` or `lognormal` |
| `--library FILE` | built-in gamma library | Per-construct parameters (see below). Required for `lognormal` |
| `--size N` | `1e6` | Cells sorted through the FACS |
| `--reads N` | `1e5` | Total sequencing reads |
| `--ratio-amplification X` | `100` | PCR amplification ratio |
| `--bias-library` | off | Unequal construct abundances (Dirichlet 0.5) |
| `--f-amp X` | `1` | Fluorescence per protein ratio (scales gamma scale / shifts lognormal mu by ln X) |
| `--bins`, `--f-max`, `--upper-bounds`, `--bin-ranges` | 12 bins, 1e5 | Bin layout, see above |
| `--seed N` | clock | RNG seed; results are identical for any `--workers` |
| `--workers N` | all cores | Threads (`-1` = all) |
| `--out-path DIR` | `out/simulation_<timestamp>` | Output folder |
| `--quiet` | off | Hide progress bars |

Numbers accept scientific notation (`1e8`). Output files: `sequencing.csv`, `cells_bins.csv`
(formats below) and `metadata_simulation.csv` (the parameters used, including the seed and bin ranges).

## `infer`

```bash
# first 100 constructs of the newest out/simulation_* folder
forecast infer

# everything, from a given folder, lognormal model, custom bins
forecast infer --data my_experiment --distribution lognormal \
    --upper-bounds 10,100,1000,10000 --count all

# constructs 500-599
forecast infer --data out/sim --first-index 500 --count 100
```

| Option | Default | Description |
| --- | --- | --- |
| `--data DIR` | newest `out/simulation_*` | Folder containing `sequencing.csv` and `cells_bins.csv` |
| `--distribution` | `gamma` | `gamma` or `lognormal` |
| `--f-max`, `--upper-bounds`, `--bin-ranges` | `--f-max 1e5` | Bin layout (bin count comes from the data) |
| `--first-index N` | `0` | First construct (0-based row) to infer |
| `--count N\|all` | `100` | How many constructs to infer |
| `--workers N` | all cores | Threads |
| `--out-path DIR` | `out/inference_<timestamp>` | Output folder |
| `--quiet` | off | Hide progress bars |

## Input file formats

All files are plain CSV. Whitespace around values is ignored, and counts may be written in
scientific notation (`1.0e+02`).

### Library (`--library`, simulate only)

Two numeric columns per construct, one row each. A header row is optional.

```csv
a,b
4.935,3.401
6.02,3.514
```

| Distribution | Column 1 | Column 2 |
| --- | --- | --- |
| `gamma` | shape (> 0) | scale (> 0) |
| `lognormal` | mu, mean of log fluorescence (any finite value) | sigma, sd of log fluorescence (> 0) |

A built-in 1018-construct gamma library is used when `--library` is omitted.

### `sequencing.csv` (infer input / simulate output)

Read counts per construct (rows) and bin (columns). The first row is a header and the first
column is the construct ID:

```csv
ID,bin_0,bin_1,bin_2
construct_a,0,12,40
construct_b,3,25,0
```

A file with no header row (all numeric) is also accepted; constructs are then numbered from 0.
Every row must have the same number of columns as `cells_bins.csv`; counts must be finite and
non-negative.

### `cells_bins.csv` (infer input / simulate output)

One row giving the number of cells sorted into each bin, in the same order as the columns of
`sequencing.csv`. A header row is optional (a leading `#` is tolerated).

```csv
bin_0,bin_1,bin_2
100000,250000,180000
```

## Output files

### `results.csv` (infer)

One row per construct. Column names use `a`/`b` (gamma shape/scale) or `mu`/`sigma` (lognormal).

| Column | Meaning |
| --- | --- |
| `ID` | Construct ID from `sequencing.csv` |
| `a_mle`, `b_mle` | Maximum-likelihood estimates (0 if no fit) |
| `a_se`, `b_se` | Standard errors of the estimates (0 if unavailable) |
| `mean`, `st_dev` | Method-of-moments mean and standard deviation of fluorescence (a quick estimate; for lognormal, in log space) |
| `inference_grade` | `1` ML fit with standard errors; `2` ML fit but the curvature was unusable (no SEs; likely a poor model fit); `3` naive: reads in a single bin only, so only `mean`/`st_dev` are given; `4` no usable reads |
| `score` | Fraction of the construct's cells in the first and last bins; high values mean the distribution is cut off by the bin range, so treat estimates with care |
| `nll` | Negative log-likelihood at the optimum |
| `a_mom`, `b_mom` | Method-of-moments parameters, on the same scale as the MLE |

Values are rounded to 3 decimals.

### Other files

- `metadata_inference.csv` / `metadata_simulation.csv`: key,value record of the run settings.
- `infer` also copies `sequencing.csv` and `cells_bins.csv` into its output folder.

## Differences from the Python version

- Flat options replace the `auto_bin` / `custom_bin` sub-commands, and `--count N|all`
  replaces `--last_index sample|all`.
- With `--upper-bounds`, `n` values give `n+1` bins (the Python custom-bin mode gave an
  inconsistent bin count).
- Standard errors use a finite-difference Hessian and the optimiser is `argmin`'s Nelder-Mead,
  so estimates agree closely with Python's but not bit-for-bit.

Dependencies: `rand`/`rand_distr`, `statrs`, `argmin`, `csv`, `indicatif`.
