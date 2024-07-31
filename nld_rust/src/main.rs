use nld::with_threads;
use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Equilibrium density
    #[arg(short='n', long="eq_den", default_value_t = 0.5)]
    eq_den: f64,

    /// Total number of samples
    #[arg(short='s', long="samples", default_value_t = 100)]
    samples: usize,

    /// Total number of time steps
    #[arg(short='t', long="t_max", default_value_t = 100)]
    t_max: usize,

    /// Total number of lattice sites 
    #[arg(short='l', long="chain_length", default_value_t = 999)]
    chain_length: usize,

    /// wave number mode of initial perturbation. It is constrained to be 0<q<l/2
    #[arg(short='q', long="wave_number", default_value_t = 99)]
    wave_number: usize,

    /// Amplitude of the intital perturbation from equilibrium density at n. n is constrained to be 0<n<1
    #[arg(short='a', long="amplitude", default_value_t = 0.1)]
    amplitude: f64,

}


fn main() {

    let args = Args::parse();

    let result = with_threads(args.samples, args.t_max, args.chain_length, args.wave_number, args.amplitude, args.eq_den);

    println!("{:?}", result.m_rho);
    println!("{:?}", result.v_rho);
    println!("{:?}", result.s_rho);
    println!("{:?}", result.k_rho);
}