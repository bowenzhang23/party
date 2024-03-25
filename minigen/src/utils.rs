pub fn get_name<F>(_: &F) -> &'static str {
    std::any::type_name::<F>()
}
