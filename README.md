# forecast

**Simulate and analyse Flow-seq / massively parallel reporter assay (MPRA) experiments.**

`forecast` is a fast, single-binary Rust port of the core of FORECAST
(Gilliot & Gorochowski, *Methods in Molecular Biology* 2553,
<https://doi.org/10.1007/978-1-0716-2617-7_3>). It does two things:

- **`simulate`** generates synthetic Flow-seq data, so you can plan an experiment (how many cells,
  reads and bins do I need?) or test an analysis.
- **`infer`** takes real or simulated Flow-seq data and estimates, for every genetic construct,
  its **expression level** and its **cell-to-cell variation**, with standard errors.

## Contents

- [Background](#background)
- [Install](#install)
- [Quick start: a worked example](#quick-start-a-worked-example)
- [Defining bins](#defining-bins)
- [`simulate`](#simulate)
- [`infer`](#infer)
- [Using your own data](#using-your-own-data)
- [Input file formats](#input-file-formats)
- [Output files and how to read them](#output-files-and-how-to-read-them)
- [Tips and troubleshooting](#tips-and-troubleshooting)
- [Differences from the Python version](#differences-from-the-python-version)

## Background

In a Flow-seq experiment a library of genetic constructs (for example promoters driving a
fluorescent reporter) is put into cells. The cells are sorted by a FACS machine into several
**bins** by fluorescence, the DNA from each bin is amplified by PCR, and each bin is sequenced.
For every construct you end up with a row of **read counts, one per bin**. If a construct is
bright, its reads pile up in the high-fluorescence bins; if it is variable, its reads spread
across several bins.

Together with the number of cells sorted into each bin, those counts are enough to estimate the
distribution of fluorescence of cells carrying each construct. `forecast` models that
distribution as either:

| Model | Parameters | Mean | Standard deviation |
| --- | --- | --- | --- |
| `gamma` (default) | shape `a`, scale `b` | `a·b` | `√a · b` |
| `lognormal` | `mu`, `sigma` of log fluorescence | `exp(mu + sigma²/2)` | `exp(mu + sigma²/2) · √(exp(sigma²) − 1)` |

**How `simulate` works.** It draws how many cells of each construct are sampled, draws each
cell's fluorescence from that construct's distribution, drops the cell into the matching bin,
allocates sequencing reads to bins in proportion to the cells sorted there, and samples reads
for each construct from the PCR-amplified pool. This reproduces both the sorting noise and the
sequencing noise of a real experiment.

**How `infer` works.** For each construct it builds a likelihood of the observed read counts
across bins given a candidate distribution, finds the maximum-likelihood parameters
(Nelder-Mead optimisation, started from a method-of-moments estimate), and estimates standard
errors from the curvature (Hessian) of the likelihood at the optimum. Constructs that cannot be
fitted properly are flagged with an *inference grade* instead of being silently trusted.

## Install

You need a recent [Rust toolchain](https://rustup.rs).

```bash
cargo build --release          # binary: target/release/forecast
cargo test --release           # optional: run the test suite
```
 
Copy `target/release/forecast` anywhere on your `PATH`, or run it from there. The examples below
assume it is called `forecast`.

General behaviour:

- `forecast --help` and `forecast <command> --help` list every option.
- Option names accept `--out-path` or `--out_path`; values can be `--seed 42` or `--seed=42`.
- Numbers can be written in scientific notation (`1e8`).
- Progress bars show on a terminal; `--quiet` hides them.
- Errors print `error: <message>` to stderr and exit with code 2. Nothing is silently ignored:
  bad bins, malformed CSVs, negative counts and so on are all reported.

## Quick start: a worked example

Simulate an experiment with the built-in library of 1018 gamma-distributed constructs, 2 million
sorted cells, 10 million reads and 8 bins:

```bash
forecast simulate --size 2e6 --reads 1e7 --bins 8 --seed 1 --out-path out/sim
```

```
simulated 1018 constructs over 8 bins -> out/sim
```

This writes three files to `out/sim`. `sequencing.csv` holds the read counts per construct and
bin, and `cells_bins.csv` the cells sorted per bin:

```
ID,bin_1,bin_2,bin_3,bin_4,bin_5,bin_6,bin_7,bin_8
0,233,8777,1052,0,0,0,0,0
1,18,7265,2091,0,0,0,0,0
2,0,0,0,22,8186,1488,0,0
...
```

Now infer the distribution of the first few constructs (`--count` defaults to 100; use
`--count all` for every construct):

```bash
forecast infer --data out/sim --count 4 --out-path out/inf
```

`out/inf/results.csv`:

```
ID,a_mle,b_mle,a_se,b_se,mean,st_dev,inference_grade,score,nll,a_mom,b_mom
0,4.791,3.509,0.212,0.15,22.681,20.661,1.0,0.023,-12395.652,1.205,18.822
1,6.987,3.045,0.668,0.27,30.901,27.868,1.0,0.002,-11228.309,1.229,25.133
...
```

The library's true values for construct 0 are `a = 4.935, b = 3.401`; the maximum-likelihood
estimate (`a_mle`, `b_mle`) recovers `4.79 ± 0.21` and `3.51 ± 0.15`. The expression level is
then `a_mle · b_mle` (≈ 16.8) and the variation (standard deviation) is `√a_mle · b_mle`
(≈ 7.7). Prefer the `_mle` columns to the quicker `mean`/`st_dev`/`_mom` columns (see
[Output files](#output-files-and-how-to-read-them)).

Omit `--data` and `infer` uses the newest `out/simulation_*` folder, so the shortest possible
workflow is just:

```bash
forecast simulate
forecast infer
```

## Defining bins

The bin layout must describe how your cells were sorted. Choose **one** of three options; the
same options work in `simulate` and `infer`, and the layout you pass to `infer` must match the
one used to produce the data.

| Option | Meaning | Example |
| --- | --- | --- |
| `--bins N --f-max X` | `N` log-spaced bins from 0 up to `X`; the last bin is open-ended | `--bins 12 --f-max 1e5` |
| `--upper-bounds A,B,C` | Contiguous bins `0:A`, `A:B`, `B:C`, `C:inf` (one more bin than values). Values must increase and the first must be > 1 | `--upper-bounds 10,100,1000` |
| `--bin-ranges LO:HI,...` | An explicit lower and upper fluorescence for every bin | `--bin-ranges 0:10,10:100,100:1000` |

**Automatic bins** (`--bins`/`--f-max`, the default of 12 bins up to 1e5) are log-spaced, which
suits the wide range of fluorescence in most experiments. In `infer` the number of bins is read
from the data, so only `--f-max` is needed.

**`--upper-bounds`** is the simplest way to enter the gate boundaries you used on the sorter.

**`--bin-ranges`** is the most flexible and the only way to describe unusual gating:

- Bins must be in increasing order and must not overlap. Giving `0:10,5:100` is an error.
- **Gaps are allowed**, e.g. `0:5,5:20,30:100`. Cells falling in the 20–30 gap are not collected
  (in simulation they are discarded, as with a real sorter).
- `inf` is accepted for the last upper bound only: `0:10,10:100,100:inf`.
- **A finite last bound is allowed**, e.g. `0:10,10:100,100:1000`. This models a hard ceiling:
  cells brighter than 1000 are discarded, and inference accounts for that truncation.
- The number of ranges must equal the number of bins (`bin_1`..`bin_N` columns in `sequencing.csv`).

For the `lognormal` model, bounds are still given in fluorescence units; `forecast` converts them
to log space internally.

## `simulate`

```bash
forecast simulate [OPTIONS]
```

Examples:

```bash
# 100 million reads, 8 automatic bins, built-in gamma library
forecast simulate --reads 1e8 --bins 8

# your sorter's gates, reproducible, using 4 threads
forecast simulate --upper-bounds 10,100,1000,10000 --seed 42 --workers 4

# explicit bin ranges with a finite top bin (cells above 20000 are lost)
forecast simulate --bin-ranges 0:5,5:20,20:100,100:1000,1000:20000

# unequal construct abundances in the library, as in a real library
forecast simulate --bias-library --size 5e6

# your own construct library, lognormal model
forecast simulate --distribution lognormal --library my_library.csv --bins 6 --f-max 1e4
```

| Option | Default | Description |
| --- | --- | --- |
| `--distribution` | `gamma` | `gamma` or `lognormal` |
| `--library FILE` | built-in gamma library | Per-construct distribution parameters ([format](#library---library)). Required for `lognormal` |
| `--size N` | `1e6` | Number of cells sorted through the FACS |
| `--reads N` | `1e5` | Total sequencing reads shared across all bins |
| `--ratio-amplification X` | `100` | PCR amplification ratio applied after sorting |
| `--bias-library` | off | Give constructs unequal abundance in the starting library (Dirichlet with concentration 0.5) instead of equal abundance |
| `--f-amp X` | `1` | Fluorescence-per-protein ratio. Multiplies the gamma scale, or adds `ln X` to the lognormal `mu` |
| `--bins`, `--f-max`, `--upper-bounds`, `--bin-ranges` | 12 bins, `1e5` | Bin layout, see [Defining bins](#defining-bins) |
| `--seed N` | clock | Random seed. The same seed gives identical output regardless of `--workers` |
| `--workers N` | all cores | Number of threads (`-1` = all cores) |
| `--out-path DIR` | `out/simulation_<timestamp>` | Output folder (created if needed) |
| `--quiet` | off | Hide progress bars |

Planning tip: run `simulate` with several values of `--size`, `--reads` and bin layouts, run
`infer` on each, and compare the standard errors (`a_se`, `b_se`) and the recovery of the known
library values to choose a design.

Output: `sequencing.csv`, `cells_bins.csv` and `metadata_simulation.csv` (all the settings used,
including the seed and the exact bin ranges, so a run can be reproduced).

## `infer`

```bash
forecast infer [OPTIONS]
```

Examples:

```bash
# first 100 constructs of the newest out/simulation_* folder
forecast infer

# every construct in a given folder
forecast infer --data out/sim --count all

# a particular slice: constructs 500-599
forecast infer --data out/sim --first-index 500 --count 100

# lognormal model with the gates used in the experiment
forecast infer --data my_experiment --distribution lognormal \
    --upper-bounds 10,100,1000,10000 --count all

# explicit bin ranges, 8 threads, no progress bars (e.g. in a script)
forecast infer --data my_experiment --bin-ranges 0:10,10:100,100:1000,1000:inf \
    --workers 8 --quiet --count all --out-path results/run1
```

| Option | Default | Description |
| --- | --- | --- |
| `--data DIR` | newest `out/simulation_*` | Folder containing `sequencing.csv` and `cells_bins.csv` |
| `--distribution` | `gamma` | Model to fit: `gamma` or `lognormal` |
| `--f-max`, `--upper-bounds`, `--bin-ranges` | `--f-max 1e5` | Bin layout. The number of bins comes from the data and must agree |
| `--first-index N` | `0` | First construct to infer (0-based row of `sequencing.csv`) |
| `--count N\|all` | `100` | How many constructs to infer |
| `--workers N` | all cores | Number of threads |
| `--out-path DIR` | `out/inference_<timestamp>` | Output folder |
| `--quiet` | off | Hide progress bars |

Inference is independent per construct, so it scales across cores. Output: `results.csv`,
`metadata_inference.csv`, and copies of the input `sequencing.csv` and `cells_bins.csv` so each
result folder is self-contained.

## Using your own data

To analyse a real experiment you need two files in one folder, plus the bin definition:

1. **`sequencing.csv`**: read counts for each construct (rows) in each sorted bin (columns
   `bin_1`..`bin_N`; other columns are ignored).
2. **`cells_bins.csv`**: the number of cells the FACS sorted into each bin (the sorter's event
   counts), in the same bin order as the columns above (see the
   [`cells_bins.csv` format](#cells_binscsv-infer-input-simulate-output)).
3. The **fluorescence range of each bin**, from your gating, passed as `--upper-bounds` or
   `--bin-ranges`. Order bins from lowest to highest fluorescence.

Then run:

```bash
forecast infer --data my_experiment \
    --upper-bounds 10,100,1000,10000 \
    --count all --out-path my_experiment/results
```

Choose `--distribution lognormal` if your reporter's fluorescence is better described on a log
scale, otherwise leave the default `gamma`. Check the grades and scores in `results.csv` (below)
before relying on individual constructs, and consider fitting both models and comparing `nll`
values on constructs you care about.

## Input file formats

All files are plain comma-separated text. Whitespace around values is ignored, and counts may
use scientific notation (`1.0e+02`). Malformed files produce a clear error naming the file and
the problem.

### Library (`--library`, simulate only)

Two numeric columns, one row per construct. A header row is optional.

```csv
a,b
4.935,3.401
6.02,3.514
```

| Model | Column 1 | Column 2 |
| --- | --- | --- |
| `gamma` | shape (> 0) | scale (> 0) |
| `lognormal` | `mu`, mean of log fluorescence (any finite value) | `sigma`, standard deviation of log fluorescence (> 0) |

If `--library` is omitted the built-in 1018-construct gamma library is used. Extra columns are
ignored; rows with missing, non-numeric, non-finite or non-positive parameters are rejected.

### `sequencing.csv` (infer input, simulate output)

Read counts per construct (rows) and bin (columns). The first row is a header and a column
named `ID` holds the construct ID (any text, no commas; if there is no `ID` column the first
column is used). Bins are numbered **from 1** and the count
columns must be named `bin_1`, `bin_2`, ... `bin_N`, lowest to highest fluorescence:

```csv
ID,bin_1,bin_2,bin_3
construct_a,0,12,40
construct_b,3,25,0
```

- `N` is the number of bins in `cells_bins.csv`. The counts are looked up **by column name**, so
  **any additional columns** (annotations, other samples, totals...) and any column order are
  fine: everything except the `ID` column and `bin_1`..`bin_N` is ignored.
- A missing `bin_k` column, or a `bin_{N+1}` column when `cells_bins.csv` has `N` values, is an
  error.
- A file with **no header** (all-numeric, counts only) is also accepted; it must then contain
  exactly one column per bin in bin order, and constructs are identified by their row number
  starting at 0.
- Counts must be finite and non-negative. Fractional counts are truncated to whole numbers.

### `cells_bins.csv` (infer input, simulate output)

One row with the number of cells sorted into each bin, in the same order as the bins in
`sequencing.csv` (`bin_1` first). The number of values sets the number of bins. A header row is
optional (a leading `#` is tolerated, as written by numpy); the header names are not used.

```csv
bin_1,bin_2,bin_3
100000,250000,180000
```

Or equivalently, without a header:

```csv
100000,250000,180000
```

Values must be finite and non-negative, with exactly one value per bin (no extra columns).

## Output files and how to read them

### `results.csv` (infer)

One row per construct. Columns use `a`/`b` (gamma shape/scale) or `mu`/`sigma` (lognormal) as
appropriate; values are rounded to 3 decimals.

| Column | Meaning |
| --- | --- |
| `ID` | Construct ID from `sequencing.csv` |
| `a_mle`, `b_mle` | **Maximum-likelihood estimates** of the two distribution parameters. These are the main result. `0` when no fit was possible |
| `a_se`, `b_se` | Standard errors of the estimates. `0` when unavailable |
| `mean`, `st_dev` | Quick method-of-moments mean and standard deviation, computed from bin midpoints. Coarse (limited by bin width), but always available for sequenced constructs. For `lognormal` these are in log space |
| `inference_grade` | Quality flag, see below |
| `score` | Fraction of the construct's estimated cells sitting in the first and last bins. High values mean the distribution is clipped by the bin range, so estimates may be biased |
| `nll` | Negative log-likelihood at the optimum (lower is a better fit; useful for comparing models on the same construct) |
| `a_mom`, `b_mom` | The method-of-moments `mean`/`st_dev` converted to the model's parameters (the optimiser's starting point) |

**Inference grades**

| Grade | Meaning | What to do |
| --- | --- | --- |
| `1` | Maximum-likelihood fit succeeded and standard errors are valid | Use it |
| `2` | Fit found, but the likelihood curvature was unusable, so no standard errors | Treat with caution: usually a poor model fit or too little data |
| `3` | Reads in a single bin only, so no likelihood fit is possible; only `mean`/`st_dev` are given | Resolution is limited by the bin width; consider narrower bins |
| `4` | No usable reads for this construct | Nothing can be inferred |

**Turning parameters into expression and variation**

- Gamma: mean = `a_mle · b_mle`, standard deviation = `√a_mle · b_mle`.
- Lognormal: mean = `exp(mu + sigma²/2)`; `sigma` is the spread on the log scale.

### Other files

- `metadata_simulation.csv` / `metadata_inference.csv`: `key,value` records of the run settings
  (distribution, seed, bin ranges, paths and so on).
- `sequencing.csv`, `cells_bins.csv`: see above. `infer` copies the inputs next to its results.

## Tips and troubleshooting

- **"N bin ranges given but the data has M bins"**: the bin option passed to `infer` doesn't
  match the number of bins (values in `cells_bins.csv`). With `--upper-bounds`, N values mean N+1 bins.
- **"bin ranges ... overlaps or precedes ..."**: bins must be listed from lowest to highest and
  must not overlap.
- **Many grade 2 or 3 results**: the data is too coarse or clipped. More bins, wider coverage of
  the fluorescence range (check `score`), more cells or more reads usually helps. Use `simulate`
  to test this before running an experiment.
- **Reproducibility**: use `--seed` with `simulate`; `infer` is deterministic.
- **Large libraries**: `infer --count all` runs in parallel across all cores; use `--workers` to
  limit it on a shared machine.

## Differences from the Python version

- Flat options replace the `auto_bin` / `custom_bin` sub-commands, and `--count N|all` replaces
  `--last_index sample|all`.
- `--bin-ranges` is new, and `--upper-bounds` with `n` values now gives `n+1` bins (the Python
  custom-bin mode gave an inconsistent bin count).
- Standard errors use a finite-difference Hessian and the optimiser is `argmin`'s Nelder-Mead,
  so estimates agree closely with the Python version but not bit-for-bit.
- Only the core simulate and infer functionality is included; the Python plotting and
  investigation scripts are not ported.

Dependencies: `rand`/`rand_distr` (sampling), `statrs` (gamma and normal distribution functions),
`argmin` (Nelder-Mead), `csv`, `indicatif` (progress bars). Licensed under MIT; see
[LICENSE](LICENSE).
