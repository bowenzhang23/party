pub trait Square {
    fn sq(self) -> Self;
}

impl Square for f64 {
    fn sq(self) -> Self {
        self * self
    }
}

pub struct C;

#[allow(dead_code)]
impl C {
    pub const PI: f64 = std::f64::consts::PI;
    pub const CONVERSION_FACTOR: f64 = 3.894e8; // pb per GeV^-2
    pub const Z_MASS: f64 = 91.188; // GeV
    pub const Z_WIDTH: f64 = 2.4414; // GeV
    pub const QED_RUNNING_COUPLING: f64 = 1.0 / 132.507;
    pub const FERMI_CONSTANT: f64 = 1.16639e-5; // GeV^-2
    pub const SIN_THETAW2: f64 = 0.222246; // Weinberg angle
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

// result in pb
pub fn ee_y_mumu(cost: f64, ecm: f64) -> f64 {
    let int_phi = 2.0 * C::PI;
    let coeff = C::CONVERSION_FACTOR * C::QED_RUNNING_COUPLING.sq() / (4.0 * ecm.sq());
    int_phi * coeff * (1.0 + cost.sq())
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
