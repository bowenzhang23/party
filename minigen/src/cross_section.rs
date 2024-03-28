use crate::utils::pdf;

/**
 * Basically do x.sq() = x * x
 */
pub trait Square {
    fn sq(self) -> Self;
}

/**
 * Implement square for f64, used as a shortcut
 */
impl Square for f64 {
    fn sq(self) -> Self {
        self * self
    }
}

/**
 * Particle Physics constants
 */
pub struct C;

/**
 * Values of Particle Physics constants
 */
impl C {
    pub const PI: f64 = std::f64::consts::PI;
    pub const SQRT_2: f64 = std::f64::consts::SQRT_2;
    pub const CONVERSION_FACTOR: f64 = 3.894e8; // pb per GeV^-2
    pub const Z_MASS: f64 = 91.188; // GeV
    pub const Z_WIDTH: f64 = 2.4414; // GeV
    pub const QED_RUNNING_COUPLING: f64 = 1.0 / 132.507;
    pub const FERMI_CONSTANT: f64 = 1.16639e-5; // GeV^-2
    pub const SIN_THETAW2: f64 = 0.222246; // Weinberg angle
}

/**
 * EW parameters for fermion
 */
pub struct Fermion {
    pub charge: f64,
    pub v: f64,
    pub a: f64,
}

/**
 * EW parameters for fermion: charge and V,A in neutral current interaction
 */
impl Fermion {
    pub fn new(label: &str) -> Self {
        let s2 = C::SIN_THETAW2;
        let lower_label = label.to_lowercase();
        let lower_label = lower_label.as_str();
        match lower_label {
            "u" | "c" | "t" => Fermion {
                charge: 2. / 3.,
                v: 0.5 - 4. / 3. * s2,
                a: 0.5,
            },
            "d" | "s" | "b" => Fermion {
                charge: -1. / 3.,
                v: -0.5 - 2. / 3. * s2,
                a: -0.5,
            },
            "ve" | "vm" | "vt" => Fermion {
                charge: 0.,
                v: 0.5,
                a: 0.5,
            },
            "e" | "mu" | "tau" => Fermion {
                charge: -1.,
                v: -0.5 + 2. * s2,
                a: -0.5,
            },
            _ => panic!("No! You are defining BSM fermions!"),
        }
    }
}

/**
 * Fermion into tuple, use as a shortcut
 */
impl Into<(f64, f64, f64)> for Fermion {
    fn into(self) -> (f64, f64, f64) {
        (self.charge, self.v, self.a)
    }
}

/**
 * Calculate the a0 and a1 in Zy matrix element
 */
pub fn drell_yan_matrix_element_constants(q: &str, l: &str, ecm: f64) -> (f64, f64) {
    let s = ecm.sq();
    let mz2 = C::Z_MASS.sq();
    let gz2 = C::Z_WIDTH.sq();
    let kappa = C::SQRT_2 * C::FERMI_CONSTANT * mz2 / (4. * C::PI * C::QED_RUNNING_COUPLING);
    let denom = ((s - mz2).sq() + gz2 * mz2).recip();
    let chi1 = kappa * s * (s - mz2) * denom;
    let chi2 = kappa.sq() * s.sq() * denom;
    let lepton = Fermion::new(l);
    let quark = Fermion::new(q);
    let (_ql, vl, al) = lepton.into();
    let (qq, vq, aq) = quark.into();
    let a0 = qq.sq() - 2. * qq * vl * vq * chi1 + (al.sq() + vl.sq()) * (aq.sq() + vq.sq()) * chi2;
    let a1 = -4. * qq * al * aq * chi1 + 8. * al * vl * aq * vq * chi2;
    (a0, a1)
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

pub fn qq_zy_mumu(cost: f64, q: &str, ecm: f64) -> f64 {
    let (a0, a1) = drell_yan_matrix_element_constants(q, "mu", ecm);
    let int_phi = 2.0 * C::PI;
    let coeff = C::CONVERSION_FACTOR * C::QED_RUNNING_COUPLING.sq() / (4.0 * ecm.sq());
    int_phi * coeff * (a0 * (1.0 + cost.sq()) + a1 * cost)
}

pub fn pp_zy_mumu(cost: f64, tau: f64, eta: f64, q: &str, ecm: f64) -> f64 {
    // only consider one parton of the proton, i.e. u
    let q_prime: String = q.to_string() + "~";
    let int_phi = 2.0 * C::PI;
    let x1 = tau.sqrt() * (eta).exp();
    let x2 = tau.sqrt() * (-eta).exp();
    let s_hat = tau * ecm.sq();
    let fq = pdf(q, x1, s_hat);
    let fq_prime = pdf(q_prime.as_str(), x2, s_hat);
    let coeff = C::CONVERSION_FACTOR * C::QED_RUNNING_COUPLING.sq() / (4.0 * s_hat);
    let (a0, a1) = drell_yan_matrix_element_constants(q, "mu", s_hat.sqrt());
    int_phi * coeff * (a0 * (1.0 + cost.sq()) + a1 * cost) * fq * fq_prime
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
