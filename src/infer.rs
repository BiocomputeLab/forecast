//! Maximum-likelihood inference of per-construct fluorescence distributions.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

use crate::model::{ms_to_ab, Distribution, Experiment};
use crate::simulate::Progress;
use crate::optim::{hessian, nelder_mead};

/// Result for one construct. For gamma, `p1/p2` are shape/scale; for lognormal, mu/sigma.
#[derive(Clone, Debug, Default)]
pub struct Inference {
    pub mle: [f64; 2],
    pub se: [f64; 2],
    pub mean: f64,
    pub st_dev: f64,
    /// 1: ML with CI, 2: ML without usable CI, 3: naive (single bin), 4: no data.
    pub grade: f64,
    pub score: f64,
    pub nll: f64,
}

impl Inference {
    pub fn mom(&self, dist: Distribution) -> (f64, f64) {
        match dist {
            Distribution::Gamma => ms_to_ab(self.mean, self.st_dev),
            Distribution::LogNormal => (self.mean, self.st_dev),
        }
    }
}

/// Empirical mean and variance of bin-centre assigned values (the method-of-moments start).
fn starting_point(i: usize, e: &Experiment) -> (f64, f64) {
    let w: Vec<f64> = e.nijhat_row(i).iter().map(|x| x.ceil()).collect();
    let (mut tmax, mut tmin) = (f64::NEG_INFINITY, f64::INFINITY);
    for (&wj, &m) in w.iter().zip(&e.mean_assigned) {
        if wj > 0.0 {
            tmax = tmax.max(m);
            tmin = tmin.min(m);
        }
    }
    if tmax == tmin {
        let j = e.mean_assigned.iter().position(|&m| m == tmax).unwrap_or(0);
        return (tmax, ((e.bins_model.hi[j] - e.bins_model.lo[j]) / 4.0).powi(2));
    }
    if tmax == 0.0 && tmin == 0.0 {
        return (0.0, 0.0);
    }
    let total: f64 = w.iter().sum();
    let mu = w.iter().zip(&e.mean_assigned).map(|(a, m)| a * m).sum::<f64>() / total;
    let var = w.iter().zip(&e.mean_assigned).map(|(a, m)| a * (m - mu).powi(2)).sum::<f64>()
        / (total - 1.0);
    (mu, var)
}

/// Negative log likelihood in log parameters theta = (ln p1, ln p2).
fn neg_ll_rep(theta: &[f64], i: usize, e: &Experiment) -> f64 {
    let (p1, p2) = (theta[0].exp(), theta[1].exp());
    let seq = e.sequencing_row(i);
    let mut nl = 0.0;
    for j in 0..e.bins {
        if e.nj[j] == 0.0 {
            continue;
        }
        let pb = e.distribution.bin_probability(e.bins_model.lo[j], e.bins_model.hi[j], p1, p2);
        let intensity = e.nihat[i] * pb * e.reads[j] / e.nj[j];
        if seq[j] != 0.0 {
            if intensity > 0.0 {
                nl += intensity - seq[j] * intensity.ln();
            }
        } else {
            nl += intensity;
        }
    }
    nl
}

fn confidence_intervals(i: usize, c: f64, d: f64, e: &Experiment) -> ([f64; 2], f64) {
    let h = hessian(|x| neg_ll_rep(x, i, e), &[c, d]);
    let finite = h.iter().flatten().all(|v| v.is_finite());
    if finite {
        let det = h[0][0] * h[1][1] - h[0][1] * h[1][0];
        // positive definite <=> leading minor and determinant positive
        if h[0][0] > 0.0 && det > 0.0 {
            let inv00 = h[1][1] / det;
            let inv11 = h[0][0] / det;
            return ([c.exp() * inv00.sqrt(), d.exp() * inv11.sqrt()], 1.0);
        }
    }
    ([0.0, 0.0], 2.0)
}

/// Infer the distribution of construct `i`.
pub fn infer_construct(i: usize, e: &Experiment) -> Inference {
    let t = e.nijhat_row(i);
    let total: f64 = t.iter().sum();
    let mut r = Inference::default();
    if total == 0.0 {
        r.grade = 4.0;
        return r;
    }
    r.score = (t[0] + t[e.bins - 1]) / total;
    let (mu, var) = starting_point(i, e);
    r.mean = mu;
    r.st_dev = var.sqrt();
    if t.iter().filter(|&&x| x != 0.0).count() == 1 {
        r.grade = 3.0;
        return r;
    }
    let start = match e.distribution {
        Distribution::LogNormal => [mu.ln(), var.sqrt().ln()],
        Distribution::Gamma => {
            let (a, b) = ms_to_ab(mu, var.sqrt());
            [a.ln(), b.ln()]
        }
    };
    if !start.iter().all(|v| v.is_finite()) {
        r.grade = 3.0;
        return r;
    }
    let x = nelder_mead(|x| neg_ll_rep(x, i, e), &start);
    let (c, d) = (x[0], x[1]);
    r.mle = [c.exp(), d.exp()];
    let (se, grade) = confidence_intervals(i, c, d, e);
    r.se = se;
    r.grade = grade;
    r.nll = neg_ll_rep(&[c, d], i, e);
    r
}

/// Infer constructs `first..last` using `workers` threads. Returns results in index order.
pub fn parallel_inference(
    first: usize,
    last: usize,
    e: &Experiment,
    workers: usize,
    progress: &dyn Progress,
) -> Vec<Inference> {
    const CHUNK: usize = 8;
    let n = last.saturating_sub(first);
    let next = AtomicUsize::new(0);
    progress.begin("Inferring ", n as u64);
    let out: Mutex<Vec<(usize, Inference)>> = Mutex::new(Vec::with_capacity(n));
    let workers = workers.clamp(1, n.max(1));
    std::thread::scope(|s| {
        for _ in 0..workers {
            s.spawn(|| loop {
                let start = next.fetch_add(CHUNK, Ordering::Relaxed);
                if start >= n {
                    break;
                }
                let end = (start + CHUNK).min(n);
                let local: Vec<(usize, Inference)> = (start..end)
                    .map(|k| {
                        let r = (k, infer_construct(first + k, e));
                        progress.advance(1);
                        r
                    })
                    .collect();
                out.lock().unwrap().extend(local);
            });
        }
    });
    let mut v = out.into_inner().unwrap();
    v.sort_by_key(|(k, _)| *k);
    v.into_iter().map(|(_, r)| r).collect()
}
