use rand::Rng;
//use rustfft::num_complex::Complex;
//use rustfft::FftPlanner;
use std::f64::consts::PI;
//use ndarray::arr2;

#[derive(Debug)]
pub struct Conf {
    st: Vec<bool>,
}

#[derive(Debug)]
pub struct TotalData {
    pub m_rho: Vec<f64>,
    pub v_rho: Vec<f64>,
    pub s_rho: Vec<f64>,
    pub k_rho: Vec<f64>,
}

pub fn rfftfreq(n: usize) -> Vec<f64> {
    (0..n / 2 + 1).map(|i| i as f64 / n as f64).collect()
}

pub fn update(config: &mut Conf, p: &[bool]) {
    for (j, &pj) in p.iter().enumerate() {
        let idx = 3 * j;
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

pub fn with_threads(samples: usize, tmax: usize, l: usize, q: usize, a: f64, n0: f64) -> TotalData {
    let freq = rfftfreq(l);
    let q_freq = freq[q+1];
    let ax: Vec<f64> = (0..l).map(|i| (q_freq * i as f64 * 2.0 * PI).cos()).collect();
    let rho_in: Vec<f64> = ax.iter().map(|&x| n0 + a * x).collect();

    let mut my_data = TotalData {
        m_rho: vec![0.0; tmax + 1],
        v_rho: vec![0.0; tmax + 1],
        s_rho: vec![0.0; tmax + 1],
        k_rho: vec![0.0; tmax + 1],
    };

    for _ in 0..samples {
        let mut rng = rand::thread_rng();
        let mut config = Conf {
            st: (0..l).map(|k| rng.gen::<f64>() <= rho_in[k]).collect(),
        };

        //let rho_0: f64 = config.st.iter().zip(&ax.iter()).map(|(&s, &ax)| if s { a } else { 0.0 }).sum::<f64>() / l as f64;
        let dot_p: f64 = config.st.iter()
            .zip(ax.iter())
            .map(|(&st, &ax)| if st { ax } else { 0.0 })
            .sum();

        let rho_0 = dot_p / l as f64;

        my_data.m_rho[0] += rho_0;
        my_data.v_rho[0] += rho_0.powi(2);
        my_data.s_rho[0] += rho_0.powi(3);
        my_data.k_rho[0] += rho_0.powi(4);

        for t in 0..tmax {
            let off = rand::thread_rng().gen_range(0..3);
            config.st.rotate_right(off);
            let p: Vec<bool> = (0..(l / 3)).map(|_| rand::random()).collect();
            //println!("{:?}", p);
            update(&mut config, &p);
            config.st.rotate_left(off);
            //let rho_t: f64 = config.st.iter().zip(&ax.iter()).map(|(&s, &ax)| if s { a } else { 0.0 }).sum::<f64>() / l as f64;
            let dot_p_t: f64 = config.st.iter()
                .zip(ax.iter())
                .map(|(&st, &ax)| if st { ax } else { 0.0 })
                .sum();
            let rho_t = dot_p_t / l as f64;
            my_data.m_rho[t + 1] += rho_t;
            my_data.v_rho[t + 1] += rho_t.powi(2);
            my_data.s_rho[t + 1] += rho_t.powi(3);
            my_data.k_rho[t + 1] += rho_t.powi(4);
        }
    }
    my_data
}