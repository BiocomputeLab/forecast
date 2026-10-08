use std::process::Command;

fn bin() -> Command {
    Command::new(env!("CARGO_BIN_EXE_forecast"))
}

#[test]
fn simulate_then_infer_recovers_parameters() {
    let dir = std::env::temp_dir().join(format!("forecast_rt_{}", std::process::id()));
    let (sim, inf) = (dir.join("sim"), dir.join("inf"));
    let ok = |c: &mut Command| assert!(c.status().unwrap().success());
    ok(bin().args(["simulate", "--reads", "1e7", "--size", "2e6", "--seed", "11", "--quiet", "--out-path"]).arg(&sim));
    ok(bin().args(["infer", "--count", "50", "--quiet", "--data"]).arg(&sim).arg("--out-path").arg(&inf));

    let lib = std::fs::read_to_string("data/library_gamma.csv").unwrap();
    let truth: Vec<(f64, f64)> = lib
        .lines()
        .skip(1)
        .map(|l| {
            let mut f = l.split(',').map(|x| x.parse::<f64>().unwrap());
            (f.next().unwrap(), f.next().unwrap())
        })
        .collect();
    let results = std::fs::read_to_string(inf.join("results.csv")).unwrap();
    let (mut good, mut n) = (0, 0);
    for (i, line) in results.lines().skip(1).enumerate() {
        let f: Vec<f64> = line.split(',').skip(1).map(|x| x.parse().unwrap()).collect();
        if f[6] == 1.0 {
            n += 1;
            let (mean_true, mean_est) = (truth[i].0 * truth[i].1, f[0] * f[1]);
            if (mean_est / mean_true - 1.0).abs() < 0.15 {
                good += 1;
            }
        }
    }
    assert!(n >= 40 && good as f64 / n as f64 > 0.8, "{good}/{n} means recovered");
    std::fs::remove_dir_all(dir).ok();
}
