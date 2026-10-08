//! CSV reading/writing helpers on top of the `csv` crate.

use std::path::Path;

pub struct Table {
    pub header: Option<Vec<String>>,
    pub rows: Vec<Vec<String>>,
}

/// Read a CSV file. A first line containing any non-numeric field is a header (a leading `#`
/// is tolerated, as written by `numpy.savetxt`). Whitespace around fields is trimmed.
pub fn read_table(path: &Path) -> Result<Table, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("cannot read {}: {e}", path.display()))?;
    parse_table(file, path)
}

/// Parse CSV text (e.g. an embedded file); `name` is used in error messages.
pub fn parse_table<R: std::io::Read>(reader: R, name: &Path) -> Result<Table, String> {
    let path = name;
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .flexible(true)
        .trim(csv::Trim::All)
        .from_reader(reader);
    let mut header = None;
    let mut rows = Vec::new();
    for (n, rec) in rdr.records().enumerate() {
        let rec = rec.map_err(|e| format!("{}: {e}", path.display()))?;
        let mut fields: Vec<String> = rec.iter().map(str::to_string).collect();
        if n == 0 {
            if let Some(f) = fields.first_mut() {
                *f = f.trim_start_matches('#').trim().to_string();
            }
        }
        if fields.iter().all(String::is_empty) {
            continue;
        }
        if n == 0 && fields.iter().any(|f| f.parse::<f64>().is_err()) {
            header = Some(fields);
        } else {
            rows.push(fields);
        }
    }
    Ok(Table { header, rows })
}

pub fn parse_f64(s: &str, path: &Path) -> Result<f64, String> {
    s.parse::<f64>().map_err(|_| format!("{}: '{s}' is not a number", path.display()))
}

/// Write rows to a CSV file.
pub fn write_rows<I, R>(path: &Path, rows: I) -> Result<(), String>
where
    I: IntoIterator<Item = R>,
    R: IntoIterator<Item = String>,
{
    let mut w = csv::WriterBuilder::new()
        .flexible(true)
        .from_path(path)
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    for r in rows {
        w.write_record(r).map_err(|e| format!("{}: {e}", path.display()))?;
    }
    w.flush().map_err(|e| format!("{}: {e}", path.display()))
}
