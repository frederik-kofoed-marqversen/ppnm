extern crate scientific;
use scientific::integration::{plain_monte_carlo};
use std::f64::consts::PI;

fn main() {
    let mut rng = scientific::rand::Rng::new(1234);
    let samples = 1e6 as u64;

    println!("PLAIN MONTE CARLO INTEGRATION (using {samples} samples)");

    let f = |x: &Vec<f64>| {
        if x.iter().fold(0.0, |sum, val| sum + val*val).sqrt() <= 1.0 {
            return 1.0
        } else {
            return 0.0
        }
    };
    let expected = PI / 4.0;
    let (int, sigma) = plain_monte_carlo(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, &mut rng);
    println!("\nIntegral of f(x, y)=1 in the unit disc in the first quadrant.");
    println!("Expected result: pi/4 ≈ {expected}");
    println!("Calculated result: {int} \nEstimated error: {sigma}");
    println!("Matches the expected: {}", (int-expected).abs() < sigma);

    let f = |x: &Vec<f64>| {
        if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
            return x[0] + x[1]
        } else {
            return 0.0
        }
    };
    let (int, sigma) = plain_monte_carlo(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, &mut rng);
    let expected = 13.0 / 20.0;
    println!("\nIntegral of f(x, y)=x+y for x>=0 above the curve y=x^2 and below y=1.");
    println!("Expected result: 13/20 ≈ {expected}");
    println!("Calculated result: {int} \nEstimated error: {sigma}");
    println!("Matches the expected: {}", (int-expected).abs() < sigma);

    let f = |x: &Vec<f64>| {
        if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
            return 1.0 / (1.0 - x[0].cos()*x[1].cos()*x[2].cos())
        } else {
            return 0.0
        }
    };
    let (mut int, sigma) = plain_monte_carlo(&f, vec![0.0; 3], vec![PI; 3], samples, &mut rng);
    int /= PI.powi(3);
    let expected = 1.3932039296856768591842462603255;
    println!("\n∫_0^π  dx/π ∫_0^π  dy/π ∫_0^π  dz/π [1-cos(x)cos(y)cos(z)]^(-1)");
    println!("Expected result: Γ(1/4)4/(4π3) ≈ {expected}");
    println!("Calculated result: {int} \nEstimated error: {sigma}");
    println!("Matches the expected: {}", (int-expected).abs() < sigma);

}