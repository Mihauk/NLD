use rand::Rng;
//use rustfft::num_complex::Complex;
//use rustfft::FftPlanner;
use std::f64::consts::PI;
use rayon::prelude::*;
//use ndarray::arr2;

#[derive(Debug)]
struct Conf {
    st: Vec<bool>,
}

#[derive(Debug)]
struct TotalData {
    m1: Vec<f64>,
    m2: Vec<f64>,
    m3: Vec<f64>,
    m4: Vec<f64>,
}

fn rfftfreq(n: usize) -> Vec<f64> {
    (0..n / 2 + 1).map(|i: usize| i as f64 / n as f64).collect()
    }

//##################################################################################################################################################
    
/// Update the configuration of the system sequentially
/* 
fn update(config: &mut Conf, p: &[bool]) {
    for (j, &pj) in p.iter().enumerate() {
        let idx: usize = 3 * j;
        if config.st[idx] && !config.st[idx + 1] && !config.st[idx + 2] {
            config.st[idx] = false;
            config.st[idx + 1] = pj;
            config.st[idx + 2] = !pj;
        } else if !config.st[idx] && config.st[idx + 1] && !config.st[idx + 2] {
            config.st[idx] = pj;
            config.st[idx + 1] = false;
            config.st[idx + 2] = !pj;
        } else if !config.st[idx] && !config.st[idx + 1] && config.st[idx + 2] {
            config.st[idx] = pj;
            config.st[idx + 1] = !pj;
            config.st[idx + 2] = false;
        }
    }
}
*/

//##################################################################################################################################################

/// Update the configuration of the system in parallel
fn update(config: &mut Conf, p: &[bool]) {
    config.st.par_chunks_mut(3).zip(p.par_iter()).for_each(|(chunk, &pj)| {
        match chunk {
            [true, false, false] => {
                chunk[0] = false;
                chunk[1] = pj;
                chunk[2] = !pj;
            }
            [false, true, false] => {
                chunk[0] = pj;
                chunk[1] = false;
                chunk[2] = !pj;
            }
            [false, false, true] => {
                chunk[0] = pj;
                chunk[1] = !pj;
                chunk[2] = false;
            }
            _ => {}
        }
    });
}

//##################################################################################################################################################


pub fn with_threads(samples: usize, tmax: usize, l: usize, q: usize, a: f64, n0: f64) -> [Vec<f64>; 4] {
    let freq: Vec<f64> = rfftfreq(l);
    let q_freq: f64 = freq[q+1];
    let ax: Vec<f64> = (0..l).into_par_iter().map(|i: usize| (q_freq * i as f64 * 2.0 * PI).cos()).collect();
    let rho_in: Vec<f64> = ax.par_iter().map(|&x| n0 + a * x).collect();

    let mut my_data: TotalData = TotalData {
        m1: vec![0.0; tmax],
        m2: vec![0.0; tmax],
        m3: vec![0.0; tmax],
        m4: vec![0.0; tmax],
    };

    for _ in 0..samples {
        let mut rng: rand::prelude::ThreadRng = rand::thread_rng();
        let mut config: Conf = Conf {
            st: (0..l).map(|k: usize| rng.gen::<f64>() <= rho_in[k]).collect(),
        };

        //let rho_0: f64 = config.st.iter().zip(&ax.iter()).map(|(&s, &ax)| if s { a } else { 0.0 }).sum::<f64>() / l as f64;
        let dot_p: f64 = config.st.par_iter()
            .zip(ax.par_iter())
            .map(|(&st, &ax)| if st { ax } else { 0.0 })
            .sum();

        let rho_0: f64 = dot_p / l as f64;

        my_data.m1[0] += rho_0;
        my_data.m2[0] += rho_0.powi(2);
        my_data.m3[0] += rho_0.powi(3);
        my_data.m4[0] += rho_0.powi(4);

        for t in 1..tmax {
            let off: usize = rand::thread_rng().gen_range(0..3);
            config.st.rotate_right(off);
            let p: Vec<bool> = (0..(l / 3)).map(|_| rand::random()).collect();
            //println!("{:?}", p);
            update(&mut config, &p);
            config.st.rotate_left(off);
            //let rho_t: f64 = config.st.iter().zip(&ax.iter()).map(|(&s, &ax)| if s { a } else { 0.0 }).sum::<f64>() / l as f64;
            let dot_p_t: f64 = config.st.par_iter()
                .zip(ax.par_iter())
                .map(|(&st, &ax)| if st { ax } else { 0.0 })
                .sum();
            let rho_t: f64 = dot_p_t / l as f64;
            my_data.m1[t] += rho_t;
            my_data.m2[t] += rho_t.powi(2);
            my_data.m3[t] += rho_t.powi(3);
            my_data.m4[t] += rho_t.powi(4);
        }
    }
    
    return [ my_data.m1.iter().map(|&x| x / samples as f64).collect::<Vec<f64>>(), my_data.m2.iter().map(|&x| x / samples as f64).collect::<Vec<f64>>(), my_data.m3.iter().map(|&x| x / samples as f64).collect::<Vec<f64>>(), my_data.m4.iter().map(|&x| x / samples as f64).collect::<Vec<f64>>() ];
}
