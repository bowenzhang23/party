pub mod cross_section;
pub mod monte_carlo;
pub mod utils;

#[cfg(test)]
mod tests {
    use crate::*;

    #[test]
    fn test_constant() {
        let f = |_x: f64| cross_section::constant(_x, Some(1.0));
        let result = monte_carlo::integral(&f, (-1.0, 1.0), None);
        println!("{} => {:?}", utils::get_name(&f), result);
        assert!((result.int - 2.0).abs() < 1e-6);
        assert!((result.err - 0.0).abs() < 1e-6);
    }

    #[test]
    fn test_breit_wigner() {
        let (m, gamma) = (10.0, 1.0);
        let range = (0.0, 20.0 * 20.0);

        let f = |x2: f64| cross_section::breit_wigner(x2, m, gamma);
        let result = monte_carlo::integral(&f, range, None);
        println!("{} => {:?}", utils::get_name(&f), result);

        let (trans, trans_inv, f_new) = cross_section::breit_wigner_is(m, gamma);
        let range_new = (trans_inv(range.0), trans_inv(range.1));
        let result_is = monte_carlo::integral(&f_new, range_new, None);
        println!("{} => {:?}", utils::get_name(&f), result_is);

        assert!((result.int - result_is.int).abs() < result.err);
        assert!((result_is.err - 0.0).abs() < 1e-6);

        let events_rho = monte_carlo::generate(&f_new, result_is, range_new, Some(10));
        let mut events_x = Vec::new();
        for event in events_rho {
            events_x.push(trans(event).sqrt());
        }
        println!("{:?}", events_x);
    }
}
