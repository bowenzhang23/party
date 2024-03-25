pub trait Square {
    fn sq(self) -> Self;
}

impl Square for f64 {
    fn sq(self) -> Self {
        self * self
    }
}

// ----------------------------------------------------------------------------
// Functions
// ----------------------------------------------------------------------------

pub fn constant(_x: f64, c: Option<f64>) -> f64 {
    c.unwrap_or(1.0)
}

pub fn breit_wigner(x2: f64, m: f64, gamma: f64) -> f64 {
    ((x2 - m.sq()).sq() + m.sq() * gamma.sq()).recip()
}

// ----------------------------------------------------------------------------
// Important Samplings
// ----------------------------------------------------------------------------

pub fn breit_wigner_is(
    m: f64,
    gamma: f64,
) -> (
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
    impl Fn(f64) -> f64,
) {
    let trans = move |rho: f64| m * gamma * rho.tan() + m.sq();
    let trans_inv = move |x2: f64| ((x2 - m.sq()) / (m * gamma)).atan();
    let f_new = move |_rho: f64| 1.0 / (m * gamma);
    return (trans, trans_inv, f_new);
}
