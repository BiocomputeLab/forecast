//! Random sampling helpers built on `rand` / `rand_distr`, including a multinomial sampler
//! (sequential conditional binomials) which neither crate provides.

use rand::rngs::StdRng;
use rand::{Rng as _, SeedableRng};
use rand_distr::{Binomial, Distribution, Gamma, Normal};

pub struct Rng(StdRng);

impl Rng {
    pub fn seed_from(seed: u64) -> Self {
        Rng(StdRng::seed_from_u64(seed))
    }

    /// Independent stream derived from a master seed, so results do not depend on thread count.
    pub fn derive(seed: u64, stream: u64) -> Self {
        let mixed = seed ^ stream.wrapping_add(1).wrapping_mul(0x9E37_79B9_7F4A_7C15).rotate_left(23);
        Self::seed_from(mixed)
    }

    pub fn normal(&mut self, mean: f64, sd: f64) -> f64 {
        match Normal::new(mean, sd) {
            Ok(d) => d.sample(&mut self.0),
            Err(_) => mean,
        }
    }

    pub fn gamma(&mut self, shape: f64, scale: f64) -> f64 {
        match Gamma::new(shape, scale) {
            Ok(d) => d.sample(&mut self.0),
            Err(_) => 0.0,
        }
    }

    pub fn binomial(&mut self, n: u64, p: f64) -> u64 {
        if n == 0 || !(p > 0.0) {
            return 0;
        }
        if p >= 1.0 {
            return n;
        }
        Binomial::new(n, p).map(|d| d.sample(&mut self.0)).unwrap_or(0)
    }

    pub fn uniform(&mut self) -> f64 {
        self.0.gen::<f64>()
    }

    /// Multinomial draw of `n` items over categories with non-negative `weights`
    /// (need not be normalised); writes the counts to `out`.
    pub fn multinomial(&mut self, n: u64, weights: &[f64], out: &mut [u64]) {
        self.multinomial_with(n, weights, out, |_| {});
    }

    /// As [`Rng::multinomial`], calling `tick(k)` as categories are processed (`k` at a time).
    pub fn multinomial_with(
        &mut self,
        n: u64,
        weights: &[f64],
        out: &mut [u64],
        mut tick: impl FnMut(u64),
    ) {
        const STEP: usize = 4096;
        debug_assert_eq!(weights.len(), out.len());
        out.iter_mut().for_each(|x| *x = 0);
        let mut remaining_w: f64 = weights.iter().sum();
        if !(remaining_w > 0.0) {
            return;
        }
        let mut remaining_n = n;
        let last = weights.iter().rposition(|&w| w > 0.0).unwrap_or(0);
        for (k, &w) in weights.iter().enumerate() {
            if k % STEP == STEP - 1 {
                tick(STEP as u64);
            }
            if remaining_n == 0 {
                continue;
            }
            if k == last {
                out[k] = remaining_n;
                remaining_n = 0;
                continue;
            }
            if w <= 0.0 {
                continue;
            }
            let x = self.binomial(remaining_n, (w / remaining_w).min(1.0));
            out[k] = x;
            remaining_n -= x;
            remaining_w -= w;
        }
        tick((weights.len() % STEP) as u64);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn multinomial_sums() {
        let mut r = Rng::seed_from(4);
        let w = [1.0, 0.0, 3.0, 6.0];
        let mut out = [0u64; 4];
        r.multinomial(100_000, &w, &mut out);
        assert_eq!(out.iter().sum::<u64>(), 100_000);
        assert_eq!(out[1], 0);
        assert!((out[3] as f64 / 100_000.0 - 0.6).abs() < 0.01);
    }

    #[test]
    fn derived_streams_are_reproducible() {
        let a = Rng::derive(7, 3).gamma(2.0, 1.0);
        let b = Rng::derive(7, 3).gamma(2.0, 1.0);
        assert_eq!(a, b);
    }
}
