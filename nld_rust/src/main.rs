use nld::with_threads;
use clap::Parser;
use hdf5::File;
use std::thread;
use std::time::Instant;
//use threadpool::ThreadPool;
//use std::sync::mpsc::channel;
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

    let args: Args = Args::parse();

    println!("Start!");
    println!("Running 1D non-linear simulation with the following parameters:");
    println!("{:?}", args);

    let n0: f64 = args.eq_den;
    let a: f64 = args.amplitude;
    let l: usize = args.chain_length;
    let q: usize = args.wave_number;
    let tmax: usize = args.t_max;
    let samples: usize = args.samples;
    let t_samples: usize = args.outer_samples;

    let filename: String = format!("./target/debug/data/rho_dotp-n_{}-A_{}-q_{}-t_samples_{}-samples_each_run_{}-tmax_{}-l-{}.h5",
                           n0, a, q, t_samples, samples, tmax, l);
    let dataset_names: [&str; 4] = ["m1", "m2", "m3", "m4"];

    let now: Instant = Instant::now();

    let mut data: [Vec<f64>; 4] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];


//############################################################################################################################################################

    // Concurrency using ThreadPool

    /*
    let n_workers: usize = 56;
    //let n_jobs: i32 = 8;
    let pool: ThreadPool = ThreadPool::new(n_workers);

    let (tx, rx) = channel();

    for _ in 0..t_samples {
        let tx = tx.clone();
        pool.execute(move|| {
            tx.send(with_threads(samples, tmax, l, q, a, n0)).expect("channel will be there waiting for the pool");
        });
    }

    for received_data in rx.iter().take(t_samples) {
        for j in 0..received_data.len() {
            if data[j].is_empty() {
                data[j].clone_from(&received_data[j])
            } else {
                data[j] = data[j].iter()
                    .zip(received_data[j].iter())
                    .map(|(&x, &y)| x + y)
                    .collect();
            }
        }
    }

    */

//############################################################################################################################################################
    
    // Concurrency using Std::Threads

    
    // here t_samples is the number of threads

    let handles: Vec<_> = (0..t_samples).map(|_| {
        thread::spawn(move || {
            with_threads(samples, tmax, l, q, a, n0)
        })
    }).collect();
        
    for handle in handles {
        let thread_data:[Vec<f64>; 4] = handle.join().unwrap();
        for j in 0..thread_data.len() {
            if data[j].is_empty() {
                data[j].clone_from(&thread_data[j])
            } else {
                data[j] = data[j].iter()
                    .zip(thread_data[j].iter())
                    .map(|(&x, &y)| x + y)
                    .collect();
            }
        }
    }
    


    //############################################################################################################################################################



    let elapsed_time: std::time::Duration = now.elapsed();

    println!("Running time: {} seconds.", elapsed_time.as_secs());

    // Write the data to a hdf5 file
    println!("Writing data to file: {}", filename);

    let file = File::create(filename).unwrap();

    for (i, &name) in dataset_names.iter().enumerate() {
        // Create a new dataset
        let dataset: hdf5::Dataset = file.new_dataset::<f64>().shape(data[i].len()).create(name).unwrap();
        // Write the combined data to the dataset
        dataset.write(&data[i]).unwrap();
    }

    println!("Done!");
}
