pub fn get_name<F>(_: &F) -> &'static str {
    std::any::type_name::<F>()
}

pub fn asym(f: f64, b: f64) -> f64 {
    (f - b) / (f + b)
}