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

/**
 * This is a fake parton distribution function set!
 * see page12 of https://pdg.lbl.gov/2019/reviews/rpp2019-rev-structure-functions.pdf
 */
pub fn pdf(flavour: &str, x: f64, _q2: f64) -> f64 {
    match flavour {
        "u" => 2. * x * (1. - x), // max at 1/3, makes sense..
        "d" => x * (1. - x), // max at 1/3, makes sense..
        "u~" => 2. - 2. * x * (1. - x).exp(), // some function..
        "d~" => 2. - 2. * x * (1. - x).exp(), // same as u~
        _ => panic!("Flavour {} not supported", flavour),
    }
}