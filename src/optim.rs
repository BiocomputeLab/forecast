//! Derivative-free minimisation (argmin Nelder-Mead) and a finite difference Hessian.

use argmin::core::{CostFunction, Executor, State};
use argmin::solver::neldermead::NelderMead;

/// Minimise `f` from `x0` with argmin's Nelder-Mead, starting from the same initial simplex
/// SciPy uses (5% perturbation of each coordinate). Non-finite values are treated as +infinity.
pub fn nelder_mead<F: Fn(&[f64]) -> f64>(f: F, x0: &[f64]) -> Vec<f64> {
    struct Problem<G>(G);
    impl<G: Fn(&[f64]) -> f64> CostFunction for Problem<G> {
        type Param = Vec<f64>;
        type Output = f64;
        fn cost(&self, p: &Vec<f64>) -> Result<f64, argmin::core::Error> {
            let v = (self.0)(p);
            Ok(if v.is_nan() { f64::INFINITY } else { v })
        }
    }

    let mut simplex = vec![x0.to_vec()];
    for k in 0..x0.len() {
        let mut y = x0.to_vec();
        y[k] = if y[k] != 0.0 { 1.05 * y[k] } else { 0.00025 };
        simplex.push(y);
    }
    let run = || -> Result<Vec<f64>, argmin::core::Error> {
        let solver = NelderMead::new(simplex).with_sd_tolerance(1e-9)?;
        let res = Executor::new(Problem(f), solver)
            .configure(|state| state.max_iters(400))
            .run()?;
        res.state().get_best_param().cloned().ok_or_else(|| argmin::core::Error::msg("no solution"))
    };
    run().unwrap_or_else(|_| x0.to_vec())
}

/// Central finite-difference Hessian with one step of Richardson extrapolation.
pub fn hessian<F: FnMut(&[f64]) -> f64>(mut f: F, x: &[f64]) -> Vec<Vec<f64>> {
    let n = x.len();
    let base: Vec<f64> = x.iter().map(|v| 1e-3 * v.abs().max(1.0)).collect();
    let f0 = f(x);

    let mut estimate = |scale: f64| -> Vec<Vec<f64>> {
        let h: Vec<f64> = base.iter().map(|b| b * scale).collect();
        let mut at = |di: &[(usize, f64)]| {
            let mut y = x.to_vec();
            for &(k, d) in di {
                y[k] += d;
            }
            f(&y)
        };
        let mut hess = vec![vec![0.0; n]; n];
        for i in 0..n {
            let (fp, fm) = (at(&[(i, h[i])]), at(&[(i, -h[i])]));
            hess[i][i] = (fp - 2.0 * f0 + fm) / (h[i] * h[i]);
            for j in (i + 1)..n {
                let fpp = at(&[(i, h[i]), (j, h[j])]);
                let fpm = at(&[(i, h[i]), (j, -h[j])]);
                let fmp = at(&[(i, -h[i]), (j, h[j])]);
                let fmm = at(&[(i, -h[i]), (j, -h[j])]);
                let v = (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j]);
                hess[i][j] = v;
                hess[j][i] = v;
            }
        }
        hess
    };
    let coarse = estimate(1.0);
    let fine = estimate(0.5);
    (0..n)
        .map(|i| (0..n).map(|j| (4.0 * fine[i][j] - coarse[i][j]) / 3.0).collect())
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimises_rosenbrock_like_bowl() {
        let x = nelder_mead(|p| (p[0] - 1.5).powi(2) + 3.0 * (p[1] + 0.5).powi(2), &[0.0, 0.0]);
        assert!((x[0] - 1.5).abs() < 1e-2 && (x[1] + 0.5).abs() < 1e-2, "{x:?}");
    }

    #[test]
    fn hessian_of_quadratic() {
        let h = hessian(|p| 2.0 * p[0] * p[0] + p[0] * p[1] + 3.0 * p[1] * p[1], &[0.7, -0.2]);
        assert!((h[0][0] - 4.0).abs() < 1e-4);
        assert!((h[1][1] - 6.0).abs() < 1e-4);
        assert!((h[0][1] - 1.0).abs() < 1e-4 && (h[1][0] - 1.0).abs() < 1e-4);
    }
}
