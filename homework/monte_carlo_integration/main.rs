extern crate scientific;

use scientific::rand::Rng;
use scientific::integration::{monte_carlo_sampling, low_descrepancy_sampling};
use std::f64::consts::PI;

fn main() {
    let mut rng = Rng::new(1234);
    let seed: u32 = 1234;
    let samples = 1e5 as u32;

    println!("MONTE CARLO INTEGRATION (using {samples} samples)");

    let f = |x: &Vec<f64>| {
        if x.iter().fold(0.0, |sum, val| sum + val*val) <= 1.0 {
            return 1.0
        } else {
            return 0.0
        }
    };
    let expected = PI / 4.0;
    println!("\nIntegral of f(x, y)=1 in the unit disc in the first quadrant.");
    println!("Analytic result: pi/4 ≈ {expected}");
    
    let (int, sigma) = monte_carlo_sampling(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, &mut rng);
    println!("Plain Monte Carlo: {int}({sigma})");
    let (int, sigma) = low_descrepancy_sampling(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, seed);
    println!("Low descrepancy: {int}({sigma})");
    

    let f = |x: &Vec<f64>| {
        if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
            return x[0] + x[1]
        } else {
            return 0.0
        }
    };
    let expected = 13.0 / 20.0;
    println!("\nIntegral of f(x, y)=x+y for x>=0 above the curve y=x^2 and below y=1.");
    println!("Analytic result: 13/20 ≈ {expected}");
    let (int, sigma) = monte_carlo_sampling(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, &mut rng);
    println!("Plain Monte Carlo: {int}({sigma})");
    let (int, sigma) = low_descrepancy_sampling(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, seed);
    println!("Low descrepancy: {int}({sigma})");


    let f = |x: &Vec<f64>| {
        1.0 / (1.0 - x[0].cos()*x[1].cos()*x[2].cos())
    };
    let expected = 1.3932039296856768591842462603255;
    println!("\n∫_0^π  dx/π ∫_0^π  dy/π ∫_0^π  dz/π [1-cos(x)cos(y)cos(z)]^(-1)");
    println!("Analytic result: Γ(1/4)4/(4π3) ≈ {expected}");
    let (mut int, sigma) = monte_carlo_sampling(&f, vec![0.0; 3], vec![PI; 3], samples, &mut rng);
    int /= PI.powi(3);
    println!("Plain Monte Carlo: {int}({sigma})");
    let (mut int, sigma) = low_descrepancy_sampling(&f, vec![0.0; 3], vec![PI; 3], samples, seed);
    int /= PI.powi(3);
    println!("Low descrepancy: {int}({sigma})");
}