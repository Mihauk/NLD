use nld::with_threads;
use clap::Parser;
use hdf5::File;
//use ndarray::Array1;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Equilibrium density
    #[arg(short='n', long="eq_den", default_value_t = 0.5)]
    eq_den: f64,

    /// Total number of samples
    #[arg(short='s', long="samples", default_value_t = 100)]
    samples: usize,

    #[arg(short='S', long="outer_samples", default_value_t = 100)]
    outer_samples: usize,

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

    println!("Start!");
    println!("Running 1D non-linear simulation with the following parameters:");
    println!("{:?}", args);

    let n0 = args.eq_den;
    let a = args.amplitude;
    let l = args.chain_length;
    let q = args.wave_number;
    let tmax = args.t_max;
    let samples = args.samples;
    let t_samples = args.outer_samples;

    let filename = format!("./target/debug/data/rho_dotp-n_{}-A_{}-q_{}-t_samples_{}-samples_each_run_{}-tmax_{}-l-{}.h5",
                           n0, a, q, t_samples, samples, tmax, l);
    let dataset_names = ["m1", "m2", "m3", "m4"];

    let data = with_threads(samples, tmax, l, q, a, n0);


    let file = File::create(filename).unwrap();

    for (i, &name) in dataset_names.iter().enumerate() {
        // Create a new dataset
        let dataset = file.new_dataset::<f64>().shape(data[i].len()).create(name).unwrap();
        // Write the combined data to the dataset
        dataset.write(&data[i]).unwrap();
    }


 /*   
    for outer in 0..=t_samples {
        let data= with_threads(samples, tmax, l, q, a, n0);

        let file = File::open_rw(&filename).unwrap_or_else(|_| File::create(&filename).unwrap());

        for (i, &dataset_name) in dataset_names.iter().enumerate() {
            let dataset = file.new_dataset::<f64>().shape(data[i].len()).create(dataset_name).unwrap();
            if outer == 0 {
                dataset.write(&data[i]).unwrap();
            } else {
                let mut previous_data = dataset.read().unwrap();
                for (j, &value) in data[i].iter().enumerate() {
                    previous_data[j] += value;
                }
                dataset.write(&previous_data).unwrap();
            }
        }
    }

    for outer in 0..=t_samples {
        // Initialize or read previous data
        let previous_data: [Vec<f64>] = if outer == 0 {
            vec![Vec<f64>::zeros(tmax); 4]
        } else {
            let file = File::open(filename).unwrap();
            dataset_names.iter().map(|&name| {
                let dataset = file.dataset(name).unwrap();
                dataset.read_1d::<f64>().unwrap()
            }).collect()
        };

        let data: [Vec<f64>; 4] = with_threads(samples, tmax, l, q, a, n0); // your computations here

        // Save the data at the end of the inner loop
        let file = File::create(filename).unwrap();
        for (i, &name) in dataset_names.iter().enumerate() {
            let dataset = file.new_dataset::<f64>().create(name).unwrap();
            let combined_data = &previous_data[i] + &data[i];
            dataset.write(&combined_data).unwrap();
        }
    }
 */

    println!("Done!");
}
