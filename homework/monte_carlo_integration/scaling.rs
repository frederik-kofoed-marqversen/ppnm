extern crate scientific;

use scientific::rand::Rng;

fn monte_carlo_sampling(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, samples: u32, rng: &mut Rng) -> (f64, f64) {
    assert!(samples > 0);
    
    let dim = a.len();
    let volume = std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));
    let (mut sum, mut sum2) = (0.0, 0.0);
    let mut x = vec![0.0; dim];
    for n in 1..=samples {
        if n % 200 == 0 {
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

fn corput_number(mut n: u32, base: u32) -> f64 {
    let mut result = 0.0;
    let mut denom = 1;
    while n > 0 {
        denom *= base;
        result += (n % base) as f64 / denom as f64;
        n /= base;
    }
    return result
}

fn halton_number(n: u32, dimension: usize) -> Vec<f64> {
    // prime numbers makes sure that all pairs are coprime
    const BASE: [u32; 19] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67];
    assert!(dimension < BASE.len());
    return BASE[..dimension].iter().map(|b| corput_number(n, *b)).collect()
}

fn low_descrepancy_sampling(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, samples: u32, seed: u32) -> (f64, f64) {
    assert!(samples > 0);
    
    let dimension = a.len();
    let volume = std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));
    let (mut sum1, mut sum2) = (0.0, 0.0);
    
    let mut x = vec![0.0; dimension];
    for (i, n) in (seed..seed+samples/2).enumerate() {
        if i*2 % 200 == 0 {
            let int1 = volume * sum1 / i as f64;
            let int2 = volume * sum2 / i as f64;
            let sigma = (int1 - int2).abs();
            println!("{} {sigma}", i*2);
        }
        
        for (i, random) in halton_number(n, dimension).iter().enumerate() {x[i] = a[i] + random*(b[i] - a[i])}
        let fx = f(&x);
        sum1 += fx;

        for (i, random) in halton_number(n+samples/2, dimension).iter().enumerate() {x[i] = a[i] + random*(b[i] - a[i])}
        let fx = f(&x);
        sum2 += fx;
    }
    let int1 = volume * sum1 / (samples/2) as f64;
    let int2 = volume * sum2 / (samples/2) as f64;
    
    return (int1, (int1 - int2).abs())
}

fn main() { 
    let f = |x: &Vec<f64>| {
        if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
            return x[0] + x[1]
        } else {
            return 0.0
        }
    };
    let (a, b) = (vec![0.0, 0.0], vec![2.0, 2.0]);


    let samples = 1e4 as u32;
    // plain monte carlo
    let mut rng = scientific::rand::Rng::new(1234);
    monte_carlo_sampling(&f, a.clone(), b.clone(), samples, &mut rng);
    println!("\n\n");
    // low discepancy
    low_descrepancy_sampling(&f, a.clone(), b.clone(), samples, 1234);
}