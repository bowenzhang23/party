use std::fs::File;
use std::io::prelude::*;
use std::path::Path;

pub fn get_name<F>(_: &F) -> &'static str {
    std::any::type_name::<F>()
}

pub fn asym(f: f64, b: f64) -> f64 {
    (f - b) / (f + b)
}

pub fn eta(t: f64) -> f64 {
    -(t / 2.0).tan().ln()
}

pub fn write_to_file<P: AsRef<Path>>(data: Vec<f64>, path: P) -> std::io::Result<()> {
    let mut f = File::create(path)?;
    for e in data {
        f.write_all(format!("{}\n", e).as_bytes())?;
    }
    Ok(())
}
