fn main() {
    std::process::exit(forecast::cli::run(std::env::args().collect()));
}
