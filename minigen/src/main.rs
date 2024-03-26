
use minigen::cross_section::*;
use minigen::monte_carlo::*;
use minigen::utils::*;

fn main() -> std::io::Result<()> {
    let mut ecm = Vec::new();
    let mut xs = Vec::new();
    let mut xs_err = Vec::new();
    let q = "u";
    let range = (-1.0, 1.0);
    
    for e in 50..=130 {
        ecm.push(e as f64);
        let f = |cost: f64| qq_zy_mumu(cost, q, e as f64);
        let result = integral(&f, range, Some(1000));
        xs.push(result.int);
        xs_err.push(result.err);
    }
    write_to_file(ecm, "outputs/uu_zy_mumu_ecm.txt")?;
    write_to_file(xs, "outputs/uu_zy_mumu_xs.txt")?;
    write_to_file(xs_err, "outputs/uu_zy_mumu_xs_err.txt")?;
    
    Ok(())
}