//! Bad input must produce a clear error message and exit code 2, never a panic.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_forecast")).args(args).arg("--quiet").output().unwrap()
}

fn scratch(name: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("forecast_bad_{}_{name}", std::process::id()));
    let _ = fs::remove_dir_all(&d);
    fs::create_dir_all(&d).unwrap();
    d
}

fn p(path: &Path) -> &str {
    path.to_str().unwrap()
}

/// Assert a clean failure whose message contains `needle`.
fn fails(args: &[&str], needle: &str) {
    let o = run(args);
    let err = String::from_utf8_lossy(&o.stderr);
    assert_eq!(o.status.code(), Some(2), "args {args:?}: stderr {err}");
    assert!(!err.contains("panicked"), "panic for {args:?}: {err}");
    assert!(err.contains(needle), "args {args:?}: expected '{needle}' in '{err}'");
}

fn data_dir(name: &str, cells: &str, seq: &str) -> PathBuf {
    let d = scratch(name);
    fs::write(d.join("cells_bins.csv"), cells).unwrap();
    fs::write(d.join("sequencing.csv"), seq).unwrap();
    d
}

const CELLS: &str = "bin_0,bin_1,bin_2\n100,200,300\n";
const SEQ: &str = "ID,bin_0,bin_1,bin_2\nc0,5,10,15\nc1,1,2,3\n";

fn infer_args<'a>(d: &'a Path, out: &'a Path, extra: &[&'a str]) -> Vec<&'a str> {
    let mut v = vec!["infer", "--data", p(d), "--out-path", p(out), "--f-max", "1000"];
    v.extend_from_slice(extra);
    v
}

#[test]
fn good_baseline_works() {
    let d = data_dir("good", CELLS, SEQ);
    let out = scratch("good_out");
    let o = run(&infer_args(&d, &out, &[]));
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
    assert!(out.join("results.csv").exists());
}

#[test]
fn overlapping_or_unsorted_custom_bins() {
    let o = scratch("o");
    for bad in ["10,5,100", "10,10,100", "10,100,50", "0.5,10", "1,10", "-5,10"] {
        fails(&["simulate", "--upper-bounds", bad, "--out-path", p(&o)], "upper bounds");
    }
    let d = data_dir("ov", CELLS, SEQ);
    fails(&infer_args(&d, &o, &["--upper-bounds", "100,50"]), "upper bounds");
}

#[test]
fn bad_bound_values() {
    let o = scratch("bv");
    fails(&["simulate", "--upper-bounds", "10,abc", "--out-path", p(&o)], "bad upper bound");
    fails(&["simulate", "--upper-bounds", "10,NaN", "--out-path", p(&o)], "upper bounds");
    fails(&["simulate", "--upper-bounds", "10,inf", "--out-path", p(&o)], "upper bounds");
}

#[test]
fn bin_count_mismatch_with_data() {
    let d = data_dir("mm", CELLS, SEQ);
    let o = scratch("mm_o");
    // 3 bins in the data need 2 upper bounds
    fails(&infer_args(&d, &o, &["--upper-bounds", "10,100,1000"]), "bins");
}

#[test]
fn bad_simulation_options() {
    let o = scratch("so");
    for (args, needle) in [
        (vec!["--bins", "0"], "--bins"),
        (vec!["--bins", "2.5"], "--bins"),
        (vec!["--f-max", "1"], "--f-max"),
        (vec!["--size", "0"], "--size"),
        (vec!["--reads", "-3"], "--reads"),
        (vec!["--size", "lots"], "not a number"),
        (vec!["--distribution", "weibull"], "unknown distribution"),
        (vec!["--distribution", "lognormal"], "--library is required"),
        (vec!["--seed", "-1"], "--seed"),
        (vec!["--bogus", "1"], "unknown option"),
        (vec!["--bins"], "needs a value"),
        (vec!["stray"], "unexpected argument"),
    ] {
        let mut a = vec!["simulate"];
        a.extend(args);
        a.extend(["--out-path", p(&o)]);
        fails(&a, needle);
    }
}

#[test]
fn bad_libraries() {
    let d = scratch("lib");
    let o = scratch("lib_o");
    for (name, body, needle) in [
        ("missing", None, "cannot read"),
        ("empty", Some(""), "no constructs"),
        ("header_only", Some("a,b\n"), "no constructs"),
        ("one_col", Some("a,b\n1\n"), "two parameters"),
        ("text", Some("a,b\n1,x\n"), "not a number"),
        ("negative", Some("a,b\n-1,2\n"), "invalid parameters"),
        ("zero", Some("a,b\n1,0\n"), "invalid parameters"),
        ("nan", Some("a,b\nNaN,2\n"), "invalid parameters"),
        ("inf", Some("a,b\n1,inf\n"), "invalid parameters"),
        ("unterminated_quote", Some("a,b\n\"1,2\n"), "lib"),
    ] {
        let path = d.join(format!("{name}.csv"));
        if let Some(b) = body {
            fs::write(&path, b).unwrap();
        }
        fails(&["simulate", "--library", p(&path), "--out-path", p(&o)], needle);
    }
}

#[test]
fn mangled_sequencing_files() {
    let o = scratch("ms_o");
    for (name, seq, needle) in [
        ("ragged", "ID,bin_0,bin_1,bin_2\nc0,5,10\n", "columns"),
        ("extra_col", "ID,bin_0,bin_1,bin_2\nc0,5,10,15,20\n", "columns"),
        ("text_count", "ID,bin_0,bin_1,bin_2\nc0,5,ten,15\n", "not a number"),
        ("negative", "ID,bin_0,bin_1,bin_2\nc0,5,-10,15\n", "non-negative"),
        ("nan", "ID,bin_0,bin_1,bin_2\nc0,5,NaN,15\n", "non-negative"),
        ("empty", "", "no sequencing data"),
        ("header_only", "ID,bin_0,bin_1,bin_2\n", "no sequencing data"),
        ("binary_junk", "\u{0}\u{1}\u{2}\n\u{ff}\n", "columns"),
    ] {
        let d = data_dir(&format!("seq_{name}"), CELLS, seq);
        fails(&infer_args(&d, &o, &[]), needle);
    }
}

#[test]
fn mangled_cells_files() {
    let o = scratch("mc_o");
    for (name, cells, needle) in [
        ("empty", "", "empty"),
        ("header_only", "bin_0,bin_1,bin_2\n", "empty"),
        ("text", "bin_0,bin_1,bin_2\n1,two,3\n", "not a number"),
        ("negative", "bin_0,bin_1,bin_2\n1,-2,3\n", "non-negative"),
        ("wrong_len", "bin_0,bin_1\n1,2\n", "columns"),
    ] {
        let d = data_dir(&format!("cells_{name}"), cells, SEQ);
        fails(&infer_args(&d, &o, &[]), needle);
    }
}

#[test]
fn missing_files_and_bad_infer_options() {
    let o = scratch("mf_o");
    let empty = scratch("mf");
    fails(&infer_args(&empty, &o, &[]), "cannot read");
    fails(&infer_args(&scratch("nonexistent_never_created").join("nope"), &o, &[]), "cannot read");
    let d = data_dir("opts", CELLS, SEQ);
    fails(&infer_args(&d, &o, &["--count", "many"]), "--count");
    fails(&infer_args(&d, &o, &["--first-index", "99"]), "beyond");
    fails(&infer_args(&d, &o, &["--distribution", "poisson"]), "unknown distribution");
    fails(&infer_args(&d, &o, &["--bogus"]), "unknown option");
}

#[test]
fn unknown_command_and_help() {
    fails(&["frobnicate"], "unknown command");
    for c in ["simulate", "infer", "help"] {
        let o = Command::new(env!("CARGO_BIN_EXE_forecast")).args([c, "--help"]).output().unwrap();
        assert!(o.status.success() && !o.stdout.is_empty());
    }
}

#[test]
fn degenerate_but_valid_data_does_not_panic() {
    // constructs with zero reads, a single occupied bin, and all-zero cell bins
    let seq = "ID,bin_0,bin_1,bin_2\nzero,0,0,0\nsingle,0,40,0\nlast,0,0,9\nfirst,7,0,0\n";
    for (name, cells) in [("deg", CELLS), ("deg0", "bin_0,bin_1,bin_2\n0,0,0\n")] {
        let d = data_dir(name, cells, seq);
        let o = scratch(&format!("{name}_o"));
        let r = run(&infer_args(&d, &o, &[]));
        assert!(r.status.success(), "{}", String::from_utf8_lossy(&r.stderr));
        assert_eq!(fs::read_to_string(o.join("results.csv")).unwrap().lines().count(), 5);
    }
}

#[test]
fn bad_bin_ranges() {
    let o = scratch("br");
    for (bad, needle) in [
        ("0:10,5:100", "overlaps"),          // overlap
        ("0:10,10:100,50:200", "overlaps"),  // later overlap
        ("10:100,0:10", "overlaps"),         // out of order
        ("10:10", "upper bound"),            // empty range
        ("20:10", "upper bound"),            // reversed
        ("-1:10", "lower bound"),            // negative
        ("0:inf,10:20", "last bin"),         // unbounded in the middle
        ("0:10,abc:20", "bad bin range"),    // not a number
        ("0-10,10-20", "LOWER:UPPER"),       // wrong syntax
        ("nan:10", "lower bound"),           // NaN
        ("", "at least one"),                // nothing
    ] {
        fails(&["simulate", "--bin-ranges", bad, "--out-path", p(&o)], needle);
    }
    fails(&["simulate", "--bin-ranges", "0:10", "--upper-bounds", "10", "--out-path", p(&o)], "not both");
    let d = data_dir("br_data", CELLS, SEQ);
    fails(&infer_args(&d, &o, &["--bin-ranges", "0:10,10:100"]), "bin ranges given");
}

#[test]
fn explicit_ranges_with_gap_roundtrip() {
    let sim = scratch("gap_sim");
    let inf = scratch("gap_inf");
    let ranges = "0:5,5:20,30:100,100:1000,1000:inf"; // gap between 20 and 30
    let ok = |c: &mut Command| assert!(c.status().unwrap().success());
    ok(Command::new(env!("CARGO_BIN_EXE_forecast"))
        .args(["simulate", "--quiet", "--seed", "5", "--size", "1e6", "--reads", "1e7", "--bin-ranges", ranges, "--out-path", p(&sim)]));
    let header = fs::read_to_string(sim.join("sequencing.csv")).unwrap();
    assert_eq!(header.lines().next().unwrap().split(',').count(), 6); // ID + 5 bins
    ok(Command::new(env!("CARGO_BIN_EXE_forecast"))
        .args(["infer", "--quiet", "--count", "20", "--bin-ranges", ranges, "--data", p(&sim), "--out-path", p(&inf)]));
    assert_eq!(fs::read_to_string(inf.join("results.csv")).unwrap().lines().count(), 21);
}
