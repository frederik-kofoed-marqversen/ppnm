extern crate scientific;

use scientific::rand::Rng;
use scientific::integration::{recursive_stratified_sampling};

fn main() {
    let samples = 5000 as u32;
    let r1 = 1.0;
    let r2 = 0.5;
    let f = |x: &Vec<f64>| {
        println!("{} {}", x[0], x[1]);
        let dist = x.iter().fold(0.0, |sum, val| sum + val*val);
        if  dist <= r1*r1 && dist >= r2*r2 {
            return 1.0
        } else {
            return 0.0
        }
    };
    let a = vec![0.0, 0.0];
    let b = vec![1.0, 1.0];
    
    let dimension = 2;
    let min_calls = 8 * dimension;  // GSL uses 16
    let estimate_fraction = 0.1;  // GSL uses 0.1
    let min_calls_per_bisection = 16 * min_calls;  // GSL uses 32
    
    let mut rng = Rng::new(1234);    
    let options = Some((estimate_fraction, min_calls, min_calls_per_bisection));
    
    recursive_stratified_sampling(&f, a, b, samples, &mut rng, options);
}