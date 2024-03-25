use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use std::fmt::Debug;

pub struct IntegralResult {
    pub int: f64,
    pub err: f64,
    pub wmax: f64,
}

impl Debug for IntegralResult {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(
            f,
            "{:.4} +- {:.4}, wmax = {:.4}",
            self.int, self.err, self.wmax
        )
    }
}

pub fn integral<F>(f: &F, range: (f64, f64), npoints: Option<usize>) -> IntegralResult
where
    F: Fn(f64) -> f64,
{
    let x1 = range.0;
    let x2 = range.1;
    let rho_gen = Uniform::new(0.0, 1.0);
    let mut rng = thread_rng();

    let n = npoints.unwrap_or(10000);
    let mut sum = 0.0;
    let mut var = 0.0;
    let mut wmax = 0.0;
    for _ in 0..n {
        let xi = (x2 - x1) * rng.sample(rho_gen) + x1;
        let wi = f(xi) * (x2 - x1);
        sum += wi;
        var += wi * wi;
        wmax = f64::max(wmax, wi);
    }
    sum /= n as f64;
    var /= n as f64;
    var -= sum * sum;
    // division by n comes from the Central Limit Theorem
    IntegralResult {
        int: sum,
        err: (var / n as f64).sqrt(),
        wmax,
    }
}

pub fn generate<F>(
    f: &F,
    result: IntegralResult,
    range: (f64, f64),
    npoints: Option<usize>,
) -> Vec<f64>
where
    F: Fn(f64) -> f64,
{
    let x1 = range.0;
    let x2 = range.1;
    let rho_gen = Uniform::new(0.0, 1.0);
    let r_gen = Uniform::new(0.0, 1.0);
    let mut rng = thread_rng();

    let n = npoints.unwrap_or(10000);
    let wmax = result.wmax;
    let mut events: Vec<f64> = Vec::new();
    let mut count = 0_usize;

    while count < n {
        let xi = (x2 - x1) * rng.sample(rho_gen) + x1;
        let wi = f(xi) * (x2 - x1);
        let wi_rel = wi / wmax;
        if wi_rel > rng.sample(r_gen) {
            count += 1;
            events.push(xi);
        }
    }
    events
}
