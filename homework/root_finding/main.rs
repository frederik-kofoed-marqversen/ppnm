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



    println!("\nSolving SE for Hydrogen atom");
    println!("Analytical bound-state energies: -0.5/n^2 = {:?}", [-0.5, -0.5/4.0, -0.5/9.0]);
    let r_min = 0.001;
    let r_max = 14.0;

    let stepper = RungeKuttaStepper::rk45();
    let driver = AdaptiveStepSizeDriver::new(stepper, Some((1e-6, 0.0)));
    let solve_schrodinger = |energy: f64| -> (Vec<f64>, Vec<Vec<f64>>) {
        let f = |r: f64, mut y: Vec<f64>| ->  Vec<f64> {
            let f = y[0];
            let fp = y[1];
            let fpp = -2.0 * (energy + 1.0/r) * f;
            y[0] = fp;
            y[1] = fpp;
            return y;
        };
        let y0 = vec![r_min-r_min*r_min, 1.0-2.0*r_min];
        return driver.run(f, r_min, y0, r_max, 0.01);
    };

    let objective = |x: &Vec<f64>| -> Vec<f64> {
        let energy = x[0];
        let k = f64::sqrt(2.0 * energy.abs());

        let (rs, ys) = solve_schrodinger(energy);
        let (r_end, y_end) = (rs.last().unwrap(), ys.last().unwrap());
        let (f, df) = (y_end[0], y_end[1]);

        return vec![df * r_end - f * (1.0 - k*r_end)];
    };

    let e0 = newton_root(&objective, vec![-1.0], Some((1000, 1e-3))).unwrap()[0];
    let e1 = newton_root(&objective, vec![-0.1], Some((1000, 1e-3))).unwrap()[0];
    let e2 = newton_root(&objective, vec![-0.01], Some((1000, 1e-3))).unwrap()[0];
    println!("Ground state energy = {e0} => n = {}", (-2.0*e0).powf(-0.5));
    println!("First excited state energy = {e1} => n = {}", (-2.0*e1).powf(-0.5));
    println!("Second excited state energy = {e2} => n = {}", (-2.0*e2).powf(-0.5));

    println!("\n--------Bound state reduced wavefunctions--------");

    for e in [e0, e1, e2] {
        println!("\n");
        let (rs, ys) = solve_schrodinger(e);
        for (r, y) in zip(rs, ys) {
            println!("{r} {}", y[0]);
        }
    }
}