pub mod cross_section;
pub mod monte_carlo;
pub mod utils;

#[cfg(test)]
mod tests {
    use crate::cross_section::*;
    use crate::monte_carlo::*;
    use crate::utils::*;

    #[test]
    fn test_constant() {
        let f = |_x: f64| constant(_x, Some(1.0));
        let result = integral(&f, (-1.0, 1.0), None);
        println!("{} => {:?}", get_name(&f), result);
        assert!((result.int - 2.0).abs() < 1e-6);
        assert!((result.err - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_breit_wigner() {
        let (m, gamma) = (10.0, 1.0);
        let range = (0.0, 20.0 * 20.0);

        let f = |x2: f64| breit_wigner(x2, m, gamma);
        let result = integral(&f, range, None);
        println!("{} => {:?}", get_name(&f), result);

        let (trans, trans_inv, f_new) = breit_wigner_is(m, gamma);
        let range_new = (trans_inv(range.0), trans_inv(range.1));
        let result_is = integral(&f_new, range_new, None);
        println!("{} => {:?}", get_name(&f), result_is);

        assert!((result.int - result_is.int).abs() < result.err);
        assert!((result_is.err - 0.0).abs() < 1e-6);

        let events_rho = generate(&f_new, result_is, range_new, Some(10));
        let mut events_x = Vec::new();
        for event in events_rho {
            events_x.push(trans(event).sqrt());
        }
        println!("{:?}", events_x);
    }

    #[test]
    fn test_ee_y_mumu() {
        let ecm: f64 = 9.0; // GeV
        let range = (-1.0, 1.0);
        let range_f = (0.0, 1.0);
        let range_b = (-1.0, 0.0);

        // total cross section
        let sigma = 4. * C::PI * C::QED_RUNNING_COUPLING.powf(2.) / 3. / ecm.powf(2.)
            * C::CONVERSION_FACTOR;
        let f = |cost: f64| ee_y_mumu(cost, ecm);
        let result = integral(&f, range, Some(100_000));
        println!("{} => {:?}, exact = {:.4}", get_name(&f), result, sigma);
        assert!((result.int - sigma).abs() < 3. * result.err);

        // distribution of pseudorapidity
        let events_cost = generate(&f, result, range, None);
        let mut events_eta = Vec::new();
        for event in events_cost {
            events_eta.push(eta(event.acos()));
        }
        write_to_file(events_eta, "outputs/ee_y_mumu_eta_dist.txt").expect("failed to write");

        // asymmetry
        let result_f = integral(&f, range_f, Some(100_000));
        let result_b = integral(&f, range_b, Some(100_000));
        println!("{} => {:?}", get_name(&f), result_f);
        println!("{} => {:?}", get_name(&f), result_b);
        let asym = asym(result_f.int, result_b.int);
        println!("asymmetry = {}", asym);
        assert!(asym.abs() < 1e-3);
    }

    #[test]
    fn test_qq_zy_mumu() {
        let ecm: f64 = 90.0; // GeV
        let s = ecm.sq();
        let range = (-1.0, 1.0);
        let range_f = (0.0, 1.0);
        let range_b = (-1.0, 0.0);
        let q = "u";

        let mz2 = C::Z_MASS.sq();
        let gz2 = C::Z_WIDTH.sq();
        let kappa = C::SQRT_2 * C::FERMI_CONSTANT * mz2 / (4. * C::PI * C::QED_RUNNING_COUPLING);
        let denom = ((s - mz2).sq() + gz2 * mz2).recip();
        let chi1 = kappa * s * (s - mz2) * denom;
        let chi2 = kappa.sq() * s.sq() * denom;
        let muon = Fermion::new("mu");
        let quark = Fermion::new(q);
        let (_qm, vm, am) = muon.into();
        let (qf, vf, af) = quark.into();
        let a0 =
            qf.sq() - 2. * qf * vm * vf * chi1 + (am.sq() + vm.sq()) * (af.sq() + vf.sq()) * chi2;

        // total cross section
        let sigma =
            4. * C::PI * C::QED_RUNNING_COUPLING.powf(2.) / 3. / s * a0 * C::CONVERSION_FACTOR;
        let f = |cost: f64| qq_zy_mumu(cost, q, ecm);
        let result = integral(&f, range, Some(100_000));
        println!("{} => {:?}, exact = {:.4}", get_name(&f), result, sigma);
        assert!((result.int - sigma).abs() < 3. * result.err);

        // asymmetry
        let result_f = integral(&f, range_f, Some(100_000));
        let result_b = integral(&f, range_b, Some(100_000));
        println!("{} => {:?}", get_name(&f), result_f);
        println!("{} => {:?}", get_name(&f), result_b);
        let asym = asym(result_f.int, result_b.int);
        println!("asymmetry = {}", asym);
        assert!(asym.abs() > 5e-2);
    }
}
