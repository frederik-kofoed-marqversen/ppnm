extern crate scientific;

use scientific::rand::Rng;

fn monte_carlo_sampling(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, samples: u32, rng: &mut Rng) -> (f64, f64) {
    assert!(samples > 0);
    
    let dim = a.len();
    let volume = std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));
    let (mut sum, mut sum2) = (0.0, 0.0);
    let mut x = vec![0.0; dim];
    for n in 0..samples {
        if n >= 100 && n % 100 == 0 {
            let mean = sum / n as f64;
            let variance = sum2/n as f64 - mean*mean;
            let sigma = volume * f64::sqrt(variance / n as f64);
            println!("{n} {sigma}");
        }

        for i in 0..dim {x[i] = a[i] + rng.f64()*(b[i] - a[i])}
        let fx = f(&x);
        sum += fx;
        sum2 += fx * fx;
    }
    let mean = sum / samples as f64;
    let variance = sum2/samples as f64 - mean*mean;
    
    return (mean * volume, volume * f64::sqrt(variance / samples as f64))
}

fn main() {
    let mut rng = scientific::rand::Rng::new(1234);
    let samples = 1e4 as u32;
    
    let f = |x: &Vec<f64>| {
        if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
            return x[0] + x[1]
        } else {
            return 0.0
        }
    };
    let (a, b) = (vec![0.0, 0.0], vec![2.0, 2.0]);

    // plain monte carlo
    monte_carlo_sampling(&f, a.clone(), b.clone(), samples, &mut rng);
}