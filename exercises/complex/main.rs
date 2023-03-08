extern crate num_traits;
extern crate num_complex;
use num_complex::Complex;
const I: Complex<f64> = Complex::<f64>::new(0.0, 1.0);
const ONE: Complex<f64> = Complex::<f64>::new(1.0, 0.0);
const PI: f64 = std::f64::consts::PI;

fn main() {
    println!("sqrt(-1) = {}", (-ONE).sqrt());
    println!("sqrt(i) = {}", I.sqrt());
    println!("e^i = {}", I.exp());
    println!("e^(i * pi) = {}", (I*PI).exp());
    println!("i^i = {}", I.powc(I));
    println!("ln(i) = {}", I.ln());
    println!("sin(i * pi) = {}", (I*PI).sin());
    println!("sinh(i * pi) = {}", (I*PI).sinh());
    println!("cosh(i * pi) = {}", (I*PI).cosh());
}