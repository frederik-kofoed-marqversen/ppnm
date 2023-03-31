extern crate sfuns;
extern crate scientific;

use sfuns::map_range;
use scientific::integration::*;
use std::f64::consts::PI;
use std::io::Write;

fn erf(z: f64) -> f64 {
    integrate(&|x: f64| f64::exp(-x*x) * 2.0 / f64::sqrt(PI), 0.0, z, Some((0.0001, 0.0))).0
}

fn print_test(file: &mut std::fs::File, f: impl Fn(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>, expected: f64, do_cc: bool) -> std::io::Result<()> {
    let (int, _err, evals) = integrate(&f, a, b, precision);
    
    file.write(format!(
        "Integrate -> {int}\n\
        \tAbsolute accuracy < 0.01: {}\n\
        \tIntegrand evaluations: {evals}\n",
        (int - expected).abs() < 0.01,
    ).as_bytes())?;

    if do_cc {
        let (int, _err, evals) = clenshaw_curtis(&f, a, b, precision);
        file.write(format!(
            "Clenshaw Curtis -> {int}\n\
            \tAbsolute accuracy < 0.01: {}\n\
            \tIntegrand evaluations: {evals}\n",
            (int - expected).abs() < 0.01,
        ).as_bytes())?;
    }

    Ok(())
}

fn main() -> std::io::Result<()> {
    // write erf.data
    let mut file = std::fs::File::create("erf.data")?;
    let num_points = 100;
    for i in 0..100 {
        let z = map_range(i as f64, &(0.0, num_points as f64), &(-3.0, 3.0));
        file.write(format!("{z}\t{}\n", erf(z)).as_bytes())?;
    }
    
    // Output file
    let mut file = std::fs::File::create("Out.txt")?;
    file.write("RESULTS FOR INTEGRATION HOMEWORK\n".as_bytes())?;

    let (a, b) = (0.0, 1.0);
    
    // Function 1
    let f = |x: f64| f64::sqrt(x);
    let expected = 2.0/3.0;
    file.write(format!("\n∫_0^1 dx √(x) = 2/3 = {expected}\n").as_bytes())?;
    print_test(&mut file, f, a, b, None, expected, true)?;

    // Function 2
    let f = |x: f64| 1.0/f64::sqrt(x);
    let expected = 2.0;
    file.write(format!("\n∫_0^1 dx 1/√(x) = {expected}\n").as_bytes())?;
    print_test(&mut file, f, a, b, None, expected, true)?;

    // Function 3
    let f = |x: f64| 4.0 * f64::sqrt(1.0 - x*x);
    let expected = PI;
    file.write(format!("\n∫_0^1 dx 4√(1-x²) = π = {expected}\n").as_bytes())?;
    print_test(&mut file, f, a, b, None, expected, true)?;

    // Function 4
    let f = |x: f64| f64::ln(x)/f64::sqrt(x);
    let expected = -4.0;
    file.write(format!("\n∫_0^1 dx 1n(x)/√(x) = {expected}\n").as_bytes())?;
    print_test(&mut file, f, a, b, Some((0.01, 0.0)), expected, true)?;

    // Integral representation of arctan
    let f = |x: f64| 1.0 / (1.0 + x * x);
    let expected = PI/2.0;
    file.write(format!("\n∫_0^∞ dx 1/(1+x^2) = π/2 = {expected}\n").as_bytes())?;
    print_test(&mut file, f, 0.0, f64::INFINITY, None, expected, false)?;

    // Polynomial improper integral
    let f = |x: f64| x.powi(-2);
    let expected = 1.0;
    file.write(format!("\n∫_1^∞ dx 1/x^2 = {expected}\n").as_bytes())?;
    print_test(&mut file, f, 1.0, f64::INFINITY, None, expected, false)?;

    // The Gaussian integral
    let f = |x: f64| f64::exp(-x*x);
    let expected = PI.sqrt();
    file.write(format!("\n∫_(-∞)^∞ dx e^(-x^2) = √π = {expected}\n").as_bytes())?;
    print_test(&mut file, f, f64::NEG_INFINITY, f64::INFINITY, None, expected, false)?;

    // Integral representation of the Gamma function
    let f = |x: f64| x.powi(4) * f64::exp(-x);
    let expected = 24.0;
    file.write(format!("\n∫_0^∞ dx x^4 * e^(-x) = Γ(4 + 1) = 4! = {expected}\n").as_bytes())?;
    print_test(&mut file, f, 0.0, f64::INFINITY, None, expected, false)?;

    return Ok(());
}