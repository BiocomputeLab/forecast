//! Flow-seq simulation: sampling constructs, sorting cells into bins, PCR and sequencing.

use crate::model::{Bins, Distribution, Range};
use crate::rng::Rng;

/// Receives progress updates; `begin` starts a stage of `total` work units, `advance` is
/// called (from any thread) as units complete.
pub trait Progress: Sync {
    fn begin(&self, label: &str, total: u64);
    fn advance(&self, units: u64);
}

pub struct NoProgress;

impl Progress for NoProgress {
    fn begin(&self, _: &str, _: u64) {}
    fn advance(&self, _: u64) {}
}

pub struct SimParams {
    pub distribution: Distribution,
    pub size: u64,
    pub reads: f64,
    pub ratio_amplification: f64,
    pub bias_library: bool,
    pub theta1: Vec<f64>,
    pub theta2: Vec<f64>,
    pub ranges: Vec<Range>,
    pub seed: u64,
    pub workers: usize,
}

pub struct SimOutput {
    pub bins: usize,
    /// Row-major diversity x bins read counts.
    pub sequencing: Vec<u64>,
    /// Cells sorted into each bin.
    pub cells_per_bin: Vec<u64>,
}

/// Run `work(index, &mut chunk_item)` over `items` split across `workers` threads.
pub fn parallel_chunks<T: Send, F>(items: &mut [T], chunk_len: usize, workers: usize, work: F)
where
    F: Fn(usize, &mut [T]) + Sync,
{
    let n_chunks = items.len().div_ceil(chunk_len.max(1));
    let workers = workers.clamp(1, n_chunks.max(1));
    if workers == 1 {
        for (i, c) in items.chunks_mut(chunk_len.max(1)).enumerate() {
            work(i, c);
        }
        return;
    }
    let per = n_chunks.div_ceil(workers);
    std::thread::scope(|s| {
        for (w, group) in items.chunks_mut(per * chunk_len).enumerate() {
            let work = &work;
            s.spawn(move || {
                for (k, c) in group.chunks_mut(chunk_len).enumerate() {
                    work(w * per + k, c);
                }
            });
        }
    });
}

pub fn simulate(p: &SimParams, progress: &dyn Progress) -> Result<SimOutput, String> {
    let diversity = p.theta1.len();
    let edges = Bins::new(p.distribution, &p.ranges)?;
    let bins = edges.len();
    if diversity == 0 || p.theta2.len() != diversity {
        return Err("library must contain at least one construct with two parameters".into());
    }

    // Step 1: relative abundance of each construct.
    let mut master = Rng::derive(p.seed, u64::MAX);
    let weights: Vec<f64> = if p.bias_library {
        (0..diversity).map(|_| master.gamma(0.5, 1.0)).collect()
    } else {
        vec![1.0; diversity]
    };

    // Step 2: number of cells sampled per construct.
    let mut ni = vec![0u64; diversity];
    master.multinomial(p.size, &weights, &mut ni);

    progress.begin("Sorting cells ", diversity as u64);
    // Step 3: sort cells into bins (cells outside the first edge are lost, as in a real sorter).
    let mut nij = vec![0u64; diversity * bins];
    parallel_chunks(&mut nij, 256 * bins, p.workers, |chunk_idx, chunk| {
        for (r, row) in chunk.chunks_mut(bins).enumerate() {
            let i = chunk_idx * 256 + r;
            let mut rng = Rng::derive(p.seed, i as u64);
            if r % 16 == 15 {
                progress.advance(16);
            }
            for _ in 0..ni[i] {
                let x = match p.distribution {
                    Distribution::LogNormal => rng.normal(p.theta1[i], p.theta2[i]),
                    Distribution::Gamma => rng.gamma(p.theta1[i], p.theta2[i]),
                };
                // cells outside every bin (below the first, or in a gap) are not collected
                let idx = edges.lo.partition_point(|&l| l <= x);
                if idx >= 1 && x < edges.hi[idx - 1] {
                    row[idx - 1] += 1;
                }
            }
        }
    });

    progress.advance(256); // rounding slack; bar is clamped at its length
    // Step 5: allocate reads proportionally to cells sorted per bin.
    let mut nj = vec![0u64; bins];
    for row in nij.chunks(bins) {
        for (t, &v) in nj.iter_mut().zip(row) {
            *t += v;
        }
    }
    let n: u64 = nj.iter().sum();
    let reads: Vec<u64> =
        nj.iter().map(|&c| (c as f64 * p.reads / (n as f64 + 0.001)).floor() as u64).collect();

    // Steps 4 and 6: PCR amplification (a uniform scaling) then multinomial read sampling per bin.
    progress.begin("Sequencing  ", (diversity * bins) as u64);
    let mut columns: Vec<Vec<u64>> = vec![vec![0; diversity]; bins];
    parallel_chunks(&mut columns, 1, p.workers, |j, col| {
        let weights: Vec<f64> =
            (0..diversity).map(|i| nij[i * bins + j] as f64 * p.ratio_amplification).collect();
        let mut rng = Rng::derive(p.seed, (1u64 << 40) + j as u64);
        rng.multinomial_with(reads[j], &weights, &mut col[0], |k| progress.advance(k));
    });
    let mut sequencing = vec![0u64; diversity * bins];
    for (j, col) in columns.iter().enumerate() {
        for (i, &v) in col.iter().enumerate() {
            sequencing[i * bins + j] = v;
        }
    }
    Ok(SimOutput { bins, sequencing, cells_per_bin: nj })
}
