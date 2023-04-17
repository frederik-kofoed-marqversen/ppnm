extern crate scientific;

use scientific::rand::Rng;
use scientific::integration::{monte_carlo_sampling, low_descrepancy_sampling, recursive_stratified_sampling};
use std::f64::consts::PI;

const SAMPLES: u32 = 10000;

fn test_integrators(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>) {
    let mut rng = Rng::new(1234);

    let (int, sigma) = monte_carlo_sampling(&f, a.clone(), b.clone(), SAMPLES, &mut rng);
    println!("Plain Monte Carlo: {int:.4}({sigma:.4})");
    let (int, sigma) = low_descrepancy_sampling(&f, a.clone(), b.clone(), SAMPLES, 1234);
    println!("Low descrepancy:   {int:.4}({sigma:.4})");
    let (int, sigma) = recursive_stratified_sampling(&f, a.clone(), b.clone(), SAMPLES, &mut rng, None);
    println!("Stratified:        {int:.4}({sigma:.4})");
}

fn main() {
    println!("MONTE CARLO INTEGRATION (using {SAMPLES} samples)");

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
    test_integrators(&f, vec![0.0, 0.0], vec![2.0, 2.0]);
    

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
    test_integrators(&f, vec![0.0, 0.0], vec![2.0, 2.0]);


    let f = |x: &Vec<f64>| {
        1.0 / (1.0 - x[0].cos()*x[1].cos()*x[2].cos())
    };
    let expected = 1.3932039296856768591842462603255 * PI.powi(3);
    println!("\n∫_0^π  dx/π ∫_0^π  dy/π ∫_0^π  dz/π [1-cos(x)cos(y)cos(z)]^(-1)");
    println!("Analytic result: Γ(1/4)4/(4π3) ≈ {expected}");
    test_integrators(&f, vec![0.0; 3], vec![PI; 3]);
}