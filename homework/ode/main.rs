extern crate sfuns;
extern crate scientific;
use sfuns::linspace;
use scientific::spline::CubicSpline;
use scientific::ode::*;

fn main() {
    // setup ODE solver
    let stepper = RungeKuttaStepper::rk45();
    let driver = AdaptiveStepSizeDriver::new(stepper);

    // Damped harmonic occilator
    let pendulum = |_t: f64, mut y: Vec<f64>| ->  Vec<f64> {
        let (b, c) = (0.25, 5.0);
        let (theta, omega) = (y[0], y[1]);
        y[0] = omega;
        y[1] = -b * omega - c * f64::sin(theta);
        return y;
    };

    let x0 = 0.0;
    let y0 = vec![std::f64::consts::PI - 0.1, 0.0];
    let endpoint = 10.0;

    let (xs, ys) = driver.run(pendulum, x0, y0, endpoint, 0.1);
    
    let theta = ys.iter().map(|item| item[0]).collect();
    let omega = ys.iter().map(|item| item[1]).collect();
    let solution = CubicSpline::new(xs.clone(), theta);
    let derivative = CubicSpline::new(xs.clone(), omega);

    for x in linspace(x0, endpoint, 101) {
        println!("{x}\t{}\t{}", solution.evaluate(x), derivative.evaluate(x));
    }
    println!("\n");

    // Lotka Volterra equation
    let lotka_volterra = |_t: f64, mut z: Vec<f64>| -> Vec<f64> {
        let (a, b, c, d) = (1.5, 1.0, 3.0, 1.0);
        let (x, y) = (z[0], z[1]);
        z[0] = a*x - b*x*y;
        z[1] = -c*y + d*x*y;
        return z;
    };

    let x0 = 0.0;
    let y0 = vec![10.0, 15.0];
    let endpoint = 15.0;

    let (xs, ys) = driver.run(lotka_volterra, x0, y0, endpoint, 0.1);
    
    let x = ys.iter().map(|item| item[0]).collect();
    let y = ys.iter().map(|item| item[1]).collect();
    let fx = CubicSpline::new(xs.clone(), x);
    let fy = CubicSpline::new(xs.clone(), y);

    for x in linspace(x0, endpoint, 300) {
        println!("{x}\t{}\t{}", fx.evaluate(x), fy.evaluate(x));
    }
}