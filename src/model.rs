//! Shared model definitions: fluorescence distributions, binning and the Flow-seq
//! `Experiment` (data plus derived quantities used by the likelihood).

use crate::special::{gamma_bin_probability, normal_bin_probability};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Distribution {
    Gamma,
    LogNormal,
}

impl Distribution {
    pub fn parse(s: &str) -> Result<Self, String> {
        match s {
            "gamma" => Ok(Distribution::Gamma),
            "lognormal" => Ok(Distribution::LogNormal),
            _ => Err(format!("unknown distribution '{s}' (expected gamma or lognormal)")),
        }
    }

    pub fn name(self) -> &'static str {
        match self {
            Distribution::Gamma => "gamma",
            Distribution::LogNormal => "lognormal",
        }
    }

    /// Probability that a cell with parameters (p1, p2) falls in (lo, hi]. For the lognormal
    /// model the edges are in log space and (p1, p2) are the mean and sd of log fluorescence.
    pub fn bin_probability(self, lo: f64, hi: f64, p1: f64, p2: f64) -> f64 {
        match self {
            Distribution::Gamma => gamma_bin_probability(lo, hi, p1, p2),
            Distribution::LogNormal => normal_bin_probability(lo, hi, p1, p2),
        }
    }
}

/// Fluorescence boundaries `[0, b1, .., b_{bins-1}]` with `bins` log-spaced points from 1 to `f_max`.
pub fn auto_boundaries(f_max: f64, bins: usize) -> Result<Vec<f64>, String> {
    if bins == 0 {
        return Err("--bins must be at least 1".into());
    }
    if !(f_max > 1.0) {
        return Err("--f-max must be greater than 1".into());
    }
    let top = f_max.log10();
    let mut b: Vec<f64> = (0..bins)
        .map(|k| if bins == 1 { 1.0 } else { 10f64.powf(top * k as f64 / (bins - 1) as f64) })
        .collect();
    b[0] = 0.0;
    Ok(b)
}

/// Boundaries from user-supplied upper bounds: `[0, ub1, .., ubn]` (n + 1 bins, last is open).
pub fn custom_boundaries(upper: &[f64]) -> Result<Vec<f64>, String> {
    if upper.is_empty() {
        return Err("--upper-bounds needs at least one value".into());
    }
    let mut b = vec![0.0];
    b.extend_from_slice(upper);
    if upper.iter().any(|x| !x.is_finite()) || b.windows(2).any(|w| !(w[1] > w[0])) || upper[0] <= 1.0 {
        return Err("upper bounds must be increasing and the first must be greater than 1".into());
    }
    Ok(b)
}

/// A fluorescence range `(lower, upper)` for one bin; `upper` may be infinite.
pub type Range = (f64, f64);

/// Contiguous ranges from boundaries `[b0, b1, .., b_{n-1}]`: bin j is `(b_j, b_{j+1})` and
/// the last bin is open-ended.
pub fn ranges_from_boundaries(b: &[f64]) -> Vec<Range> {
    (0..b.len()).map(|j| (b[j], b.get(j + 1).copied().unwrap_or(f64::INFINITY))).collect()
}

/// Parse `lo:hi,lo:hi,..` (e.g. `0:10,10:100,100:inf`).
pub fn parse_ranges(text: &str) -> Result<Vec<Range>, String> {
    let num = |t: &str| {
        t.trim().parse::<f64>().map_err(|_| format!("bad bin range value '{}'", t.trim()))
    };
    text.split(',')
        .filter(|r| !r.trim().is_empty())
        .map(|r| {
            let (lo, hi) = r
                .split_once(':')
                .ok_or_else(|| format!("bin range '{}' must look like LOWER:UPPER", r.trim()))?;
            Ok((num(lo)?, num(hi)?))
        })
        .collect()
}

/// Validate bin ranges: finite non-negative lower bound below the upper bound, bins ordered
/// and non-overlapping (gaps between bins are allowed: cells there are simply not sorted).
pub fn validate_ranges(ranges: &[Range]) -> Result<(), String> {
    if ranges.is_empty() {
        return Err("at least one bin range is required".into());
    }
    for (j, &(lo, hi)) in ranges.iter().enumerate() {
        if !lo.is_finite() || lo < 0.0 || hi.is_nan() || !(hi > lo) {
            return Err(format!(
                "bin ranges: bin {j} ({lo}:{hi}) must have a finite lower bound >= 0 below its upper bound"
            ));
        }
    }
    if ranges[..ranges.len() - 1].iter().any(|r| r.1.is_infinite()) {
        return Err("bin ranges: only the last bin can be unbounded above".into());
    }
    if let Some(j) = (1..ranges.len()).find(|&j| ranges[j].0 < ranges[j - 1].1) {
        return Err(format!(
            "bin ranges: bin {j} ({}:{}) overlaps or precedes bin {} ({}:{}); bins must be in increasing order",
            ranges[j].0, ranges[j].1, j - 1, ranges[j - 1].0, ranges[j - 1].1
        ));
    }
    Ok(())
}

/// Bin edges in the space of the model (log space for lognormal, where 0 stays 0).
pub struct Bins {
    pub lo: Vec<f64>,
    pub hi: Vec<f64>,
}

impl Bins {
    pub fn new(dist: Distribution, ranges: &[Range]) -> Result<Bins, String> {
        validate_ranges(ranges)?;
        let tf = |x: f64| match dist {
            Distribution::Gamma => x,
            Distribution::LogNormal if x == 0.0 => 0.0,
            Distribution::LogNormal => x.ln(),
        };
        let bins = Bins {
            lo: ranges.iter().map(|r| tf(r.0)).collect(),
            hi: ranges.iter().map(|r| tf(r.1)).collect(),
        };
        if bins.lo.iter().zip(&bins.hi).any(|(l, h)| !(h > l)) {
            return Err("bin ranges are too narrow after the log transform".into());
        }
        Ok(bins)
    }

    pub fn len(&self) -> usize {
        self.lo.len()
    }

    pub fn is_empty(&self) -> bool {
        self.lo.is_empty()
    }
}

/// Flow-seq experiment with the derived quantities used by the likelihood.
pub struct Experiment {
    pub bins: usize,
    pub diversity: usize,
    /// Cells sorted into each bin.
    pub nj: Vec<f64>,
    /// Reads per bin after normalisation.
    pub reads: Vec<f64>,
    /// Row-major diversity x bins read counts after normalisation.
    pub sequencing: Vec<f64>,
    pub distribution: Distribution,
    pub bins_model: Bins,
    pub mean_assigned: Vec<f64>,
    /// Estimated number of cells per construct and bin (row-major, integer valued).
    pub nijhat: Vec<f64>,
    /// Estimated number of cells per construct.
    pub nihat: Vec<f64>,
}

impl Experiment {
    pub fn new(
        nj: Vec<f64>,
        sequencing: Vec<f64>,
        distribution: Distribution,
        ranges: &[Range],
    ) -> Result<Self, String> {
        let bins = nj.len();
        if bins == 0 || sequencing.len() % bins != 0 {
            return Err("sequencing matrix does not match the number of bins".into());
        }
        if ranges.len() != bins {
            return Err(format!("{} bin ranges given but the data has {bins} bins", ranges.len()));
        }
        let diversity = sequencing.len() / bins;
        let bins_model = Bins::new(distribution, ranges)?;
        let mean_assigned: Vec<f64> = (0..bins)
            .map(|j| {
                let (lo, hi) = (bins_model.lo[j], bins_model.hi[j]);
                if hi.is_finite() { (lo + hi) / 2.0 } else { lo }
            })
            .collect();

        let mut reads = vec![0.0; bins];
        for row in sequencing.chunks(bins) {
            for (r, s) in reads.iter_mut().zip(row) {
                *r += s;
            }
        }
        let enrich: Vec<f64> =
            nj.iter().zip(&reads).map(|(&n, &r)| if r != 0.0 { n / r } else { 0.0 }).collect();
        let nijhat: Vec<f64> = sequencing
            .chunks(bins)
            .flat_map(|row| row.iter().zip(&enrich).map(|(s, e)| (s * e).trunc()))
            .collect();
        let nihat: Vec<f64> = nijhat.chunks(bins).map(|r| r.iter().sum()).collect();

        // Normalise read counts so confidence interval widths reflect estimated cell numbers.
        let mut sequencing = sequencing;
        for (row, &nh) in sequencing.chunks_mut(bins).zip(&nihat) {
            let total: f64 = row.iter().sum();
            if total > nh {
                row.iter_mut().for_each(|s| *s = *s * nh / total);
            }
        }
        let mut reads = vec![0.0; bins];
        for row in sequencing.chunks(bins) {
            for (r, s) in reads.iter_mut().zip(row) {
                *r += s;
            }
        }
        Ok(Experiment {
            bins,
            diversity,
            nj,
            reads,
            sequencing,
            distribution,
            bins_model,
            mean_assigned,
            nijhat,
            nihat,
        })
    }

    pub fn nijhat_row(&self, i: usize) -> &[f64] {
        &self.nijhat[i * self.bins..(i + 1) * self.bins]
    }

    pub fn sequencing_row(&self, i: usize) -> &[f64] {
        &self.sequencing[i * self.bins..(i + 1) * self.bins]
    }
}

/// Gamma (mean, sd) -> (shape, scale); zero where undefined.
pub fn ms_to_ab(m: f64, s: f64) -> (f64, f64) {
    let shape = if s != 0.0 { m * m / (s * s) } else { 0.0 };
    let scale = if m != 0.0 { s * s / m } else { 0.0 };
    (shape, scale)
}
