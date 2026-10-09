//! Command line interface: argument parsing and the `simulate` / `infer` commands.

use std::collections::{HashMap, HashSet};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use crate::csvio::{parse_f64, parse_table, read_table, write_rows};
use crate::infer::parallel_inference;
use crate::model::{
    auto_boundaries, custom_boundaries, parse_ranges, ranges_from_boundaries, Distribution, Experiment, Range,
};
use crate::simulate::{simulate, Progress, SimParams};

/// Terminal progress bars (hidden when `quiet`, and automatically when stderr is not a TTY).
struct Bars {
    bar: std::sync::Mutex<Option<indicatif::ProgressBar>>,
    quiet: bool,
}

impl Bars {
    fn new(quiet: bool) -> Self {
        Bars { bar: std::sync::Mutex::new(None), quiet }
    }

    fn finish(&self) {
        if let Some(b) = self.bar.lock().unwrap().take() {
            b.finish();
        }
    }
}

impl Progress for Bars {
    fn begin(&self, label: &str, total: u64) {
        self.finish();
        if self.quiet {
            return;
        }
        let bar = indicatif::ProgressBar::new(total);
        bar.set_style(
            indicatif::ProgressStyle::with_template(
                "{prefix} [{bar:40.cyan/blue}] {pos}/{len} ({percent}%) {elapsed_precise} ETA {eta}",
            )
            .unwrap()
            .progress_chars("=> "),
        );
        bar.set_prefix(label.to_string());
        *self.bar.lock().unwrap() = Some(bar);
    }

    fn advance(&self, units: u64) {
        if let Some(b) = self.bar.lock().unwrap().as_ref() {
            b.inc(units);
        }
    }
}

const DEFAULT_GAMMA_LIBRARY: &str = include_str!("../data/library_gamma.csv");

const USAGE: &str = "\
forecast - design and analysis of massively parallel reporter assays (Flow-seq)

USAGE:
    forecast <COMMAND> [OPTIONS]

COMMANDS:
    simulate    Generate a synthetic Flow-seq data set
    infer       Infer construct expression levels and variation from Flow-seq data
    help        Show help (also: forecast <COMMAND> --help)
";

const SIMULATE_USAGE: &str = "\
forecast simulate - generate a synthetic Flow-seq data set

OPTIONS:
    --distribution <gamma|lognormal>  Fluorescence distribution            [default: gamma]
    --library <FILE>                  CSV with two parameters per construct (gamma: shape,
                                      scale; lognormal: mu,sigma). Built-in gamma library
                                      is used when omitted for gamma.
    --size <N>                        Cells sorted through the FACS         [default: 1e6]
    --reads <N>                       Total sequencing reads                [default: 1e5]
    --ratio-amplification <X>         PCR amplification ratio               [default: 100]
    --bias-library                    Draw unequal construct abundances (Dirichlet 0.5)
    --f-amp <X>                       Fluorescence/protein ratio            [default: 1]
    --bins <N>                        Number of log-spaced bins (automatic) [default: 12]
    --f-max <X>                       Max fluorescence of the FACS (auto)   [default: 1e5]
    --upper-bounds <A,B,..>           Custom bin upper bounds (instead of --bins/--f-max);
                                      gives one more bin than values (last bin is open)
    --quiet                           Suppress progress bars
    --bin-ranges <LO:HI,..>           Explicit range per bin, e.g. 0:10,10:100,100:inf. Bins
                                      must be ordered and non-overlapping; gaps are allowed
                                      (cells in a gap are discarded). Only the last bin may
                                      be unbounded.
    --seed <N>                        RNG seed (results are reproducible)   [default: time]
    --workers <N>                     Threads, -1 = all cores               [default: -1]
    --out-path <DIR>                  Output folder   [default: out/simulation_<timestamp>]

Writes sequencing.csv, cells_bins.csv and metadata_simulation.csv.
";

const INFER_USAGE: &str = "\
forecast infer - infer construct expression levels and variation

OPTIONS:
    --data <DIR>                      Folder with sequencing.csv and cells_bins.csv
                                      [default: newest out/simulation_*]
    --distribution <gamma|lognormal>  Fluorescence distribution            [default: gamma]
    --f-max <X>                       Max fluorescence of the FACS (automatic bins)
                                                                           [default: 1e5]
    --upper-bounds <A,B,..>           Custom bin upper bounds; must give bins-1 values
    --bin-ranges <LO:HI,..>           Explicit range per bin (one per data column), as in
                                      `simulate`
    --first-index <N>                 Index of the first construct          [default: 0]
    --count <N|all>                   Number of constructs to infer         [default: 100]
    --workers <N>                     Threads, -1 = all cores               [default: -1]
    --quiet                           Suppress progress bars
    --out-path <DIR>                  Output folder   [default: out/inference_<timestamp>]

Writes results.csv, metadata_inference.csv and copies of the input data.
Inference grades: 1 ML with CI, 2 ML without usable CI, 3 naive (one bin), 4 no data.
";

struct Args {
    opts: HashMap<String, String>,
    flags: HashSet<String>,
}

impl Args {
    fn parse(argv: &[String], valued: &[&str], flags: &[&str]) -> Result<Args, String> {
        let mut a = Args { opts: HashMap::new(), flags: HashSet::new() };
        let mut i = 0;
        while i < argv.len() {
            let raw = argv[i].strip_prefix("--").ok_or_else(|| format!("unexpected argument '{}'", argv[i]))?;
            let (name, inline) = match raw.split_once('=') {
                Some((n, v)) => (n.replace('_', "-"), Some(v.to_string())),
                None => (raw.replace('_', "-"), None),
            };
            if flags.contains(&name.as_str()) {
                a.flags.insert(name);
            } else if valued.contains(&name.as_str()) {
                let mut value = match inline {
                    Some(v) => v,
                    None => {
                        i += 1;
                        argv.get(i)
                            .filter(|v| !v.starts_with("--"))
                            .cloned()
                            .ok_or_else(|| format!("--{name} needs a value"))?
                    }
                };
                if name == "upper-bounds" {
                    // allow space separated values
                    while argv.get(i + 1).map_or(false, |n| !n.starts_with("--")) {
                        i += 1;
                        value = format!("{value},{}", argv[i]);
                    }
                }
                a.opts.insert(name, value);
            } else {
                return Err(format!("unknown option '--{name}'"));
            }
            i += 1;
        }
        Ok(a)
    }

    fn num(&self, name: &str, default: f64) -> Result<f64, String> {
        match self.opts.get(name) {
            None => Ok(default),
            Some(v) => v.parse::<f64>().map_err(|_| format!("--{name}: '{v}' is not a number")),
        }
    }

    fn get(&self, name: &str) -> Option<&String> {
        self.opts.get(name)
    }

    fn workers(&self) -> Result<usize, String> {
        let w = self.num("workers", -1.0)?;
        Ok(if w < 1.0 {
            std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
        } else {
            w as usize
        })
    }

    fn distribution(&self) -> Result<Distribution, String> {
        Distribution::parse(self.get("distribution").map_or("gamma", String::as_str))
    }

    /// Bin ranges from `--bin-ranges` or `--upper-bounds`; `None` means automatic bins.
    fn custom_ranges(&self) -> Result<Option<Vec<Range>>, String> {
        match (self.get("bin-ranges"), self.upper_bounds()?) {
            (Some(_), Some(_)) => Err("use either --bin-ranges or --upper-bounds, not both".into()),
            (Some(r), None) => Ok(Some(parse_ranges(r)?)),
            (None, Some(ub)) => Ok(Some(ranges_from_boundaries(&custom_boundaries(&ub)?))),
            (None, None) => Ok(None),
        }
    }

    fn upper_bounds(&self) -> Result<Option<Vec<f64>>, String> {
        self.get("upper-bounds")
            .map(|s| {
                s.split(',')
                    .filter(|x| !x.trim().is_empty())
                    .map(|x| x.trim().parse::<f64>().map_err(|_| format!("bad upper bound '{x}'")))
                    .collect()
            })
            .transpose()
    }
}

fn timestamp() -> String {
    // seconds since epoch rendered as UTC YYYYMMDD-HHMMSS (no date crate needed)
    let secs = SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0) as i64;
    let (days, rem) = (secs.div_euclid(86_400), secs.rem_euclid(86_400));
    let z = days + 719_468;
    let era = z.div_euclid(146_097);
    let doe = z.rem_euclid(146_097);
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = yoe + era * 400 + i64::from(m <= 2);
    format!("{y:04}{m:02}{d:02}-{:02}{:02}{:02}", rem / 3600, rem % 3600 / 60, rem % 60)
}

fn out_path(a: &Args, prefix: &str) -> Result<PathBuf, String> {
    let p = a
        .get("out-path")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(format!("out/{prefix}_{}", timestamp())));
    fs::create_dir_all(&p).map_err(|e| format!("cannot create {}: {e}", p.display()))?;
    Ok(p)
}

fn write_metadata(path: &Path, rows: Vec<(String, String)>) -> Result<(), String> {
    write_rows(path, rows.into_iter().map(|(k, v)| [k, v]))
}

fn fmt_ranges(r: &[Range]) -> String {
    r.iter().map(|(lo, hi)| format!("{lo}:{hi}")).collect::<Vec<_>>().join(",")
}

fn run_simulate(argv: &[String]) -> Result<(), String> {
    let a = Args::parse(
        argv,
        &[
            "distribution", "library", "size", "reads", "ratio-amplification", "f-amp", "bins",
            "f-max", "upper-bounds", "bin-ranges", "seed", "workers", "out-path",
        ],
        &["bias-library", "quiet"],
    )?;
    let dist = a.distribution()?;
    let (size, reads) = (a.num("size", 1e6)?, a.num("reads", 1e5)?);
    let ratio = a.num("ratio-amplification", 1e2)?;
    let f_amp = a.num("f-amp", 1.0)?;
    if !(size >= 1.0 && reads >= 1.0 && ratio > 0.0 && f_amp >= 1.0) {
        return Err("--size and --reads must be >= 1, --ratio-amplification > 0, --f-amp >= 1".into());
    }

    let library = match a.get("library") {
        Some(p) => {
            let t = read_table(Path::new(p))?;
            (t, PathBuf::from(p))
        }
        None if dist == Distribution::Gamma => {
            let name = PathBuf::from("<built-in gamma library>");
            (parse_table(DEFAULT_GAMMA_LIBRARY.as_bytes(), &name)?, name)
        }
        None => return Err("--library is required for the lognormal distribution".into()),
    };
    let (mut theta1, mut theta2) = (Vec::new(), Vec::new());
    for row in &library.0.rows {
        if row.len() < 2 {
            return Err(format!("{}: each row needs two parameters", library.1.display()));
        }
        theta1.push(parse_f64(&row[0], &library.1)?);
        theta2.push(parse_f64(&row[1], &library.1)?);
    }
    if theta1.is_empty() {
        return Err(format!("{}: library contains no constructs", library.1.display()));
    }
    let valid = |t1: f64, t2: f64| {
        t1.is_finite() && t2.is_finite() && t2 > 0.0 && (dist == Distribution::LogNormal || t1 > 0.0)
    };
    if let Some(i) = (0..theta1.len()).find(|&i| !valid(theta1[i], theta2[i])) {
        return Err(format!(
            "{}: row {} has invalid parameters ({}, {}); values must be finite and positive",
            library.1.display(), i + 1, theta1[i], theta2[i]
        ));
    }
    match dist {
        Distribution::Gamma => theta2.iter_mut().for_each(|t| *t *= f_amp),
        Distribution::LogNormal => theta1.iter_mut().for_each(|t| *t += f_amp.ln()),
    }

    let ranges = match a.custom_ranges()? {
        Some(r) => r,
        None => {
            let f_max = a.num("f-max", 1e5)?;
            let bins = a.num("bins", 12.0)?;
            if bins < 1.0 || bins.fract() != 0.0 {
                return Err("--bins must be a positive integer".into());
            }
            ranges_from_boundaries(&auto_boundaries(f_max, bins as usize)?)
        }
    };
    let seed = match a.get("seed") {
        Some(s) => s.parse::<u64>().map_err(|_| "--seed must be a non-negative integer".to_string())?,
        None => SystemTime::now().duration_since(UNIX_EPOCH).map(|d| d.as_nanos() as u64).unwrap_or(0),
    };
    let out = out_path(&a, "simulation")?;

    let bars = Bars::new(a.flags.contains("quiet"));
    let result = simulate(&SimParams {
        distribution: dist,
        size: size as u64,
        reads,
        ratio_amplification: ratio,
        bias_library: a.flags.contains("bias-library"),
        theta1: theta1.clone(),
        theta2,
        ranges: ranges.clone(),
        seed,
        workers: a.workers()?,
    }, &bars);
    bars.finish();
    let result = result?;

    let bins = result.bins;
    let bin_names: Vec<String> = (0..bins).map(|j| format!("bin_{}", j + 1)).collect();
    let header = std::iter::once("ID".to_string()).chain(bin_names.iter().cloned());
    let body = result.sequencing.chunks(bins).enumerate().map(|(i, row)| {
        std::iter::once(i.to_string()).chain(row.iter().map(u64::to_string)).collect::<Vec<_>>()
    });
    write_rows(
        &out.join("sequencing.csv"),
        std::iter::once(header.collect::<Vec<_>>()).chain(body),
    )?;
    write_rows(
        &out.join("cells_bins.csv"),
        [bin_names.clone(), result.cells_per_bin.iter().map(u64::to_string).collect()],
    )?;
    write_metadata(
        &out.join("metadata_simulation.csv"),
        vec![
            ("distribution".into(), dist.name().into()),
            ("size".into(), size.to_string()),
            ("reads".into(), reads.to_string()),
            ("ratio_amplification".into(), ratio.to_string()),
            ("bias_library".into(), a.flags.contains("bias-library").to_string()),
            ("library".into(), library.1.display().to_string()),
            ("out_path".into(), out.display().to_string()),
            ("f_amp".into(), f_amp.to_string()),
            ("bins".into(), bins.to_string()),
            ("seed".into(), seed.to_string()),
            ("bin_ranges".into(), fmt_ranges(&ranges)),
        ],
    )?;
    println!("simulated {} constructs over {bins} bins -> {}", theta1.len(), out.display());
    Ok(())
}

fn latest_simulation() -> Result<PathBuf, String> {
    let entries = fs::read_dir("out").map_err(|_| "no --data given and no out/ folder found".to_string())?;
    entries
        .flatten()
        .filter(|e| e.file_name().to_string_lossy().starts_with("simulation_"))
        .filter_map(|e| Some((e.metadata().ok()?.modified().ok()?, e.path())))
        .max()
        .map(|(_, p)| p)
        .ok_or_else(|| "no --data given and no out/simulation_* folder found".to_string())
}

fn round3(x: f64) -> f64 {
    if x.is_finite() { (x * 1000.0).round() / 1000.0 } else { x }
}

fn run_infer(argv: &[String]) -> Result<(), String> {
    let a = Args::parse(
        argv,
        &["data", "distribution", "f-max", "upper-bounds", "bin-ranges", "first-index", "count", "workers", "out-path"],
        &["quiet"],
    )?;
    let dist = a.distribution()?;
    let data = match a.get("data") {
        Some(p) => PathBuf::from(p),
        None => latest_simulation()?,
    };
    println!("flow seq data from: {}", data.display());

    let cb_path = data.join("cells_bins.csv");
    let cb = read_table(&cb_path)?;
    let nj: Vec<f64> = cb
        .rows
        .first()
        .ok_or("cells_bins.csv is empty")?
        .iter()
        .map(|x| parse_f64(x, &cb_path))
        .collect::<Result<_, _>>()?;
    let bins = nj.len();
    if nj.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(format!("{}: cell counts must be finite and non-negative", cb_path.display()));
    }

    let sq_path = data.join("sequencing.csv");
    let sq = read_table(&sq_path)?;
    // With a header, columns are found by name, in any order: the ID from the column named ID
    // (falling back to the first column) and the counts from bin_1..bin_N (N = number of bins in
    // cells_bins.csv); any other columns are ignored.
    // Without a header the file must be counts only, one column per bin.
    let has_ids = sq.header.is_some();
    let cols: Vec<usize> = match &sq.header {
        Some(h) if h.iter().any(|c| *c == format!("bin_{}", bins + 1)) => {
            return Err(format!(
                "{}: has more bin columns than the {bins} in cells_bins.csv (found bin_{})",
                sq_path.display(),
                bins + 1
            ))
        }
        Some(h) => (1..=bins)
            .map(|k| {
                let name = format!("bin_{k}");
                h.iter().position(|c| *c == name).ok_or_else(|| {
                    format!(
                        "{}: no column '{name}' in the header; expected an ID column plus columns bin_1..bin_{bins}",
                        sq_path.display()
                    )
                })
            })
            .collect::<Result<_, _>>()?,
        None => (0..bins).collect(),
    };
    let id_col = sq.header.as_ref().and_then(|h| h.iter().position(|c| c == "ID")).unwrap_or(0);
    let needed = cols.iter().copied().chain([id_col]).max().map_or(0, |m| m + 1);
    let mut ids = Vec::with_capacity(sq.rows.len());
    let mut seq = Vec::with_capacity(sq.rows.len() * bins);
    for (r, row) in sq.rows.iter().enumerate() {
        if row.len() < needed || (!has_ids && row.len() != bins) {
            let expected = if has_ids { format!("at least {needed}") } else { bins.to_string() };
            return Err(format!("{}: row {} has {} columns, expected {expected}", sq_path.display(), r + 1, row.len()));
        }
        ids.push(if has_ids { row[id_col].clone() } else { r.to_string() });
        for &c in &cols {
            let v = parse_f64(&row[c], &sq_path)?;
            if !v.is_finite() || v < 0.0 {
                return Err(format!("{}: row {}: read counts must be finite and non-negative", sq_path.display(), r + 1));
            }
            seq.push(v.trunc());
        }
    }

    let ranges = match a.custom_ranges()? {
        Some(r) => r,
        None => ranges_from_boundaries(&auto_boundaries(a.num("f-max", 1e5)?, bins)?),
    };
    if seq.is_empty() {
        return Err(format!("{}: no sequencing data", sq_path.display()));
    }
    let experiment = Experiment::new(nj, seq, dist, &ranges)?;

    let first = a.num("first-index", 0.0)? as usize;
    let last = match a.get("count").map(String::as_str) {
        Some("all") => experiment.diversity,
        other => {
            let c = match other {
                Some(s) => s.parse::<usize>().map_err(|_| "--count must be an integer or 'all'".to_string())?,
                None => 100,
            };
            (first + c).min(experiment.diversity)
        }
    };
    if first >= experiment.diversity {
        return Err(format!("--first-index {first} is beyond the {} constructs", experiment.diversity));
    }

    let bars = Bars::new(a.flags.contains("quiet"));
    let results = parallel_inference(first, last, &experiment, a.workers()?, &bars);
    bars.finish();

    let out = out_path(&a, "inference")?;
    let (n1, n2) = match dist {
        Distribution::Gamma => ("a", "b"),
        Distribution::LogNormal => ("mu", "sigma"),
    };
    let header: Vec<String> =
        format!("ID,{n1}_mle,{n2}_mle,{n1}_se,{n2}_se,mean,st_dev,inference_grade,score,nll,{n1}_mom,{n2}_mom")
            .split(',')
            .map(str::to_string)
            .collect();
    let body = results.iter().enumerate().map(|(k, r)| {
        let (m1, m2) = r.mom(dist);
        let vals = [
            r.mle[0], r.mle[1], r.se[0], r.se[1], r.mean, r.st_dev, r.grade, r.score, r.nll, m1, m2,
        ];
        std::iter::once(ids[first + k].clone())
            .chain(vals.iter().map(|&v| format!("{:?}", round3(v))))
            .collect::<Vec<_>>()
    });
    write_rows(&out.join("results.csv"), std::iter::once(header).chain(body))?;
    for f in ["cells_bins.csv", "sequencing.csv"] {
        fs::copy(data.join(f), out.join(f)).map_err(|e| e.to_string())?;
    }
    write_metadata(
        &out.join("metadata_inference.csv"),
        vec![
            ("distribution".into(), dist.name().into()),
            ("data".into(), data.display().to_string()),
            ("out_path".into(), out.display().to_string()),
            ("first_index".into(), first.to_string()),
            ("last_index".into(), last.to_string()),
            ("bin_ranges".into(), fmt_ranges(&ranges)),
        ],
    )?;
    let _ = std::io::stdout().flush();
    println!("inferred {} constructs -> {}", results.len(), out.display());
    Ok(())
}

/// Entry point; returns the process exit code.
pub fn run(argv: Vec<String>) -> i32 {
    let rest = if argv.len() > 2 { &argv[2..] } else { &[][..] };
    let wants_help = rest.iter().any(|a| a == "--help" || a == "-h");
    let result = match argv.get(1).map(String::as_str) {
        Some("simulate") if wants_help => Ok(print!("{SIMULATE_USAGE}")),
        Some("infer") if wants_help => Ok(print!("{INFER_USAGE}")),
        Some("simulate") => run_simulate(rest),
        Some("infer") => run_infer(rest),
        None | Some("help") | Some("--help") | Some("-h") => Ok(print!("{USAGE}")),
        Some(c) => Err(format!("unknown command '{c}'\n\n{USAGE}")),
    };
    match result {
        Ok(()) => 0,
        Err(e) => {
            eprintln!("error: {e}");
            2
        }
    }
}
