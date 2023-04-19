extern crate scientific;
extern crate matrix;
use std::iter::zip;
use matrix::linalg::optimisation::newton_root;
use scientific::ode::*;

fn main() {
    println!("Finding minimum of f(x,y) = (1-x)^2+100(y-x^2)^2");
    println!("Expected result = [1.0, 1.0]");
    let f = |x: &Vec<f64>| -> Vec<f64> {
        vec![
            -2.0*(1.0-x[0]) - 400.0*x[0]*(x[1]-x[0]*x[0]), 
            200.0*(x[1]-x[0]*x[0])
        ]
    };
    let x = vec![0.0, 0.0];
    let result = newton_root(&f, x, None).unwrap();
    println!("Calculated result = {result:?}");



    println!("\nNow solving SE for Hydrogen atom");
    println!("True ground state energy = -0.5");
    let r_max = 8.0;
    let r_min = 0.001;

    let stepper = RungeKuttaStepper::rk45();
    let driver = AdaptiveStepSizeDriver::new(stepper, Some((1e-3, 1e-3)));
    let solve_schodinger = |energy: f64| -> (Vec<f64>, Vec<f64>) {
        let f = |r: f64, mut y: Vec<f64>| ->  Vec<f64> {
            let f = y[0];
            let fp = y[1];
            let fpp = -2.0 * (1.0/r + energy) * f;
            y[0] = fp;
            y[1] = fpp;
            return y;
        };
        let y0 = vec![r_min-r_min*r_min, 1.0-2.0*r_min];
        let (rs, ys) = driver.run(f, r_min, y0, r_max, 0.01);
        return (rs, ys.iter().map(|item| item[0]).collect())
    };

    let f_at_r_max = |x: &Vec<f64>| -> Vec<f64> {
        let energy = x[0];
        let (_, fs) = solve_schodinger(energy);
        return vec![*fs.last().unwrap()]
    };

    let e0 = newton_root(&f_at_r_max, vec![-1.0], Some((1000, 1e-3))).unwrap()[0];
    println!("Ground state energy = {e0}");

    println!("\n\nGround state f-function");
    let (rs, fs) = solve_schodinger(e0);
    for (r, f) in zip(rs, fs) {
        println!("{r} {f}");
    }
}