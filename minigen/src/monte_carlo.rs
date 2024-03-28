use rand::distributions::Uniform;
use rand::{thread_rng, Rng};
use std::fmt::Debug;

/**
 * Maintain results of the integration
 */
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

// ----------------------------------------------------------------------------
// 1D
// ----------------------------------------------------------------------------

/**
 * Integrate a 1D function, given the range and number of MC points
 */
pub fn integral<F>(f: &F, range_min: f64, range_max: f64, npoints: Option<usize>) -> IntegralResult
where
    F: Fn(f64) -> f64,
{
    let x1 = range_min;
    let x2 = range_max;
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

/**
 * Generate events with one var, given integration result, range and number
 * of events.
 */
pub fn generate<F>(
    f: &F,
    result: IntegralResult,
    range_min: f64,
    range_max: f64,
    npoints: Option<usize>,
) -> Vec<f64>
where
    F: Fn(f64) -> f64,
{
    let x1 = range_min;
    let x2 = range_max;
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

// ----------------------------------------------------------------------------
// Multi-D
// ----------------------------------------------------------------------------

/**
 * Integrate a N-D function, given the range and number of MC points.
 * The range must be given by N->1 function, the first m vars are used to
 * determine the integration limits of the m+1 th var.
 * The integration order of the variables is indicated by the array index.
 */
pub fn integral_ndim<F, Fr, const N: usize>(
    f: &F,
    range_min: [Fr; N],
    range_max: [Fr; N],
    npoints: Option<usize>,
) -> IntegralResult
where
    F: Fn([f64; N]) -> f64,
    Fr: Fn([Option<f64>; N]) -> f64,
{
    let rho_gen = Uniform::new(0.0, 1.0);
    let mut rng = thread_rng();

    let n = npoints.unwrap_or(10000);
    let mut sum = 0.0;
    let mut var = 0.0;
    let mut wmax = 0.0;
    for _ in 0..n {
        let mut wi = 1.;
        let mut xi: [Option<f64>; N] = [None; N];
        for i in 0..N {
            let x2 = range_max.get(i).unwrap()(xi);
            let x1 = range_min.get(i).unwrap()(xi);
            xi[i] = Some((x2 - x1) * rng.sample(rho_gen) + x1);
            wi *= x2 - x1;
        }
        wi *= f(xi.map(|x| x.unwrap()));
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

/**
 * Generate events with N vars, given integration result, range and number
 * of events. See `integrate_ndim` for the definition of `range`.
 */
pub fn generate_ndim<F, Fr, const N: usize>(
    f: &F,
    result: IntegralResult,
    range_min: [Fr; N],
    range_max: [Fr; N],
    npoints: Option<usize>,
) -> Vec<[f64; N]>
where
    F: Fn([f64; N]) -> f64,
    Fr: Fn([Option<f64>; N]) -> f64,
{
    let rho_gen = Uniform::new(0.0, 1.0);
    let r_gen = Uniform::new(0.0, 1.0);
    let mut rng = thread_rng();

    let n = npoints.unwrap_or(10000);
    let wmax = result.wmax;
    let mut events: Vec<[f64; N]> = Vec::new();
    let mut count = 0_usize;

    while count < n {
        let mut wi = 1.;
        let mut xi: [Option<f64>; N] = [None; N];
        for i in 0..N {
            let x2 = range_max.get(i).unwrap()(xi);
            let x1 = range_min.get(i).unwrap()(xi);
            xi[i] = Some((x2 - x1) * rng.sample(rho_gen) + x1);
            wi *= x2 - x1;
        }
        let xi_v = xi.map(|x| x.unwrap());
        wi *= f(xi_v);
        let wi_rel = wi / wmax;
        if wi_rel > rng.sample(r_gen) {
            count += 1;
            events.push(xi_v);
        }
    }
    events
}
