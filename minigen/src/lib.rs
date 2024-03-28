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
        let result = integral(&f, -1.0, 1.0, None);
        println!("{} => {:?}", get_name(&f), result);
        assert!((result.int - 2.0).abs() < 1e-6);
        assert!((result.err - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_breit_wigner() {
        let (m, gamma) = (10.0, 1.0);
        let range = (0.0, 20.0 * 20.0);

        let f = |x2: f64| breit_wigner(x2, m, gamma);
        let result = integral(&f, range.0, range.1, None);
        println!("{} => {:?}", get_name(&f), result);

        let (trans, trans_inv, f_new) = breit_wigner_is(m, gamma);
        let range_new = (trans_inv(range.0), trans_inv(range.1));
        let result_is = integral(&f_new, range_new.0, range_new.1, None);
        println!("{} => {:?}", get_name(&f), result_is);

        assert!((result.int - result_is.int).abs() < 3. * result.err);
        assert!((result_is.err - 0.0).abs() < 1e-6);

        let events_rho = generate(&f_new, result_is, range_new.0, range_new.1, Some(10));
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
        let result = integral(&f, range.0, range.1, Some(100_000));
        println!("{} => {:?}, exact = {:.4}", get_name(&f), result, sigma);
        assert!((result.int - sigma).abs() < 3. * result.err);

        // distribution of pseudorapidity
        let events_cost = generate(&f, result, range.0, range.1, None);
        let mut events_eta = Vec::new();
        for event in events_cost {
            events_eta.push(eta(event.acos()));
        }
        write_to_file(events_eta, "outputs/ee_y_mumu_eta_dist.txt").expect("failed to write");

        // asymmetry
        let result_f = integral(&f, range_f.0, range_f.1, Some(100_000));
        let result_b = integral(&f, range_b.0, range_b.1, Some(100_000));
        println!("{} => {:?}", get_name(&f), result_f);
        println!("{} => {:?}", get_name(&f), result_b);
        let asym = asym(result_f.int, result_b.int);
        println!("asymmetry = {}", asym);
        assert!(asym.abs() < 1e-3);
    }

    #[test]
    fn test_qq_zy_mumu() {
        const ECM: f64 = 13600.0; // GeV
        const S: f64 = ECM * ECM;
        const Q: &str = "u";
        let range = (-1.0, 1.0);
        let range_f = (0.0, 1.0);
        let range_b = (-1.0, 0.0);
        let (a0, _a1) = drell_yan_matrix_element_constants(Q, "mu", ECM);

        // total cross section
        let sigma =
            4. * C::PI * C::QED_RUNNING_COUPLING.powf(2.) / 3. / S * a0 * C::CONVERSION_FACTOR;
        let f = |cost: f64| qq_zy_mumu(cost, Q, ECM);
        let result = integral(&f, range.0, range.1, Some(100_000));
        println!("{} => {:?}, exact = {:.4}", get_name(&f), result, sigma);
        assert!((result.int - sigma).abs() < 3. * result.err);

        // asymmetry
        let result_f = integral(&f, range_f.0, range_f.1, Some(100_000));
        let result_b = integral(&f, range_b.0, range_b.1, Some(100_000));
        println!("{} => {:?}", get_name(&f), result_f);
        println!("{} => {:?}", get_name(&f), result_b);
        let asym = asym(result_f.int, result_b.int);
        println!("asymmetry = {}", asym);
        assert!(asym.abs() > 5e-2);
    }

    #[test]
    fn test_integral_ndim() {
        let f = |x: [f64; 3]| x[0] + x[1].powf(2.) + x[2].powf(3.);
        let range_min = [|_| 2., |_| 1., |_| 0.];
        let range_max = [|_| 3., |_| 2., |_| 1.];
        let result = integral_ndim(&f, range_min, range_max, None);
        println!("{:?}", result);
        let analytic = 2.5 + 7. / 3. + 0.25;
        assert!((result.int - analytic).abs() < 3. * result.err);
    }

    #[test]
    fn test_generate_ndim() {
        let f = |x: [f64; 3]| x[0] + x[1].powf(2.) + x[2].powf(3.);
        let range_min = [|_| 2., |_| 1., |_| 0.];
        let range_max = [|_| 3., |_| 2., |_| 1.];
        let result = integral_ndim(&f, range_min, range_max, None);
        let events = generate_ndim(&f, result, range_min, range_max, Some(5));
        assert_eq!(events.len(), 5);
        assert_eq!(events[0].len(), 3);
    }

    #[test]
    fn test_pp_zy_mumu() -> std::io::Result<()> {
        const ECM: f64 = 1000.; // GeV
        const CUTOFF: f64 = 60.; // GeV
        const S: f64 = ECM * ECM;
        const Q: &str = "u";
        let range_min = [
            |_| -1.,
            |_| CUTOFF * CUTOFF / S,
            |x: [Option<f64>; 3]| 0.5 * x[1].unwrap().ln(),
        ]; // variables are 0:cost, 1:tau, 2: eta
        let range_max = [
            |_| 1.,
            |_| 1.,
            |x: [Option<f64>; 3]| -0.5 * x[1].unwrap().ln(),
        ]; // variables are 0:cost, 1:tau, 2: eta

        // total cross section
        let f = |x: [f64; 3]| pp_zy_mumu(x[0], x[1], x[2], Q, ECM);
        let result = integral_ndim(&f, range_min, range_max, Some(100_000));
        println!("{} => {:?}", get_name(&f), result);

        let events = generate_ndim(&f, result, range_min, range_max, Some(10_000));
        let mut cost: Vec<f64> = Vec::new();
        let mut beta: Vec<f64> = Vec::new();
        let mut ecm_hat: Vec<f64> = Vec::new();
        for event in events {
            cost.push(event[0]);
            let x1 = event[1].sqrt() * event[2].exp();
            let x2 = event[1].sqrt() * (-event[2]).exp();
            beta.push((x2 - x1) / (x2 + x1));
            ecm_hat.push((event[1] * S).sqrt());
        }
        write_to_file(cost, "outputs/pp_zy_mumu_cost.txt")?;
        write_to_file(beta, "outputs/pp_zy_mumu_beta.txt")?;
        write_to_file(ecm_hat, "outputs/pp_zy_mumu_ecm_hat.txt")?;
        Ok(())
    }
}
