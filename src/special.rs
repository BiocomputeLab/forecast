//! Bin-probability helpers for the gamma and normal fluorescence models, built on `statrs`.

use statrs::function::erf::erfc;
use statrs::function::gamma::{gamma_lr, gamma_ur};

const SQRT_2: f64 = std::f64::consts::SQRT_2;

fn gamma_p(shape: f64, x: f64) -> f64 {
    if x <= 0.0 {
        0.0
    } else if x.is_infinite() {
        1.0
    } else {
        gamma_lr(shape, x)
    }
}

fn norm_cdf(z: f64) -> f64 {
    if z == f64::INFINITY { 1.0 } else if z == f64::NEG_INFINITY { 0.0 } else { 0.5 * erfc(-z / SQRT_2) }
}

fn norm_sf(z: f64) -> f64 {
    if z == f64::INFINITY { 0.0 } else if z == f64::NEG_INFINITY { 1.0 } else { 0.5 * erfc(z / SQRT_2) }
}

/// P(lo < X <= hi) for X ~ Gamma(shape, scale); `hi` may be infinite.
pub fn gamma_bin_probability(lo: f64, hi: f64, shape: f64, scale: f64) -> f64 {
    let (xl, xh) = (lo / scale, hi / scale);
    if xl > shape {
        // upper tail: differencing survival functions keeps precision
        let q = |x: f64| if x.is_infinite() { 0.0 } else { gamma_ur(shape, x) };
        q(xl) - q(xh)
    } else {
        gamma_p(shape, xh) - gamma_p(shape, xl)
    }
}

/// P(lo < X <= hi) for X ~ Normal(mean, sd); `hi` may be infinite.
pub fn normal_bin_probability(lo: f64, hi: f64, mean: f64, sd: f64) -> f64 {
    let (zl, zh) = ((lo - mean) / sd, (hi - mean) / sd);
    if zl > 0.0 {
        norm_sf(zl) - norm_sf(zh)
    } else {
        norm_cdf(zh) - norm_cdf(zl)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, tol: f64) {
        assert!((a - b).abs() <= tol * b.abs(), "{a} vs {b}");
    }

    #[test]
    fn bin_probabilities() {
        close(gamma_bin_probability(0.0, f64::INFINITY, 2.5, 3.0), 1.0, 1e-12);
        close(normal_bin_probability(f64::NEG_INFINITY, f64::INFINITY, 1.0, 2.0), 1.0, 1e-12);
        close(gamma_bin_probability(1.0, 2.0, 1.0, 1.0), (-1.0f64).exp() - (-2.0f64).exp(), 1e-10);
        close(gamma_bin_probability(10.0, 20.0, 1.0, 1.0), (-10.0f64).exp() - (-20.0f64).exp(), 1e-9);
        close(normal_bin_probability(0.0, f64::INFINITY, 0.0, 1.0), 0.5, 1e-12);
        close(normal_bin_probability(1.96, f64::INFINITY, 0.0, 1.0), 1.0 - 0.9750021048517795, 1e-9);
    }
}
