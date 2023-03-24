extern crate scientific;
use scientific::ode::*;
use std::iter::zip;

fn main() {
    // setup ODE solver
    let stepper = RungeKuttaStepper::rk45();
    let driver = AdaptiveStepSizeDriver::new(stepper, None);
    // Passing 'None' to the driver the tolerances take default values abs=rel=0.01

    // Damped harmonic occilator
    let pendulum = |_t: f64, mut y: Vec<f64>| ->  Vec<f64> {
        let (b, c) = (0.25, 5.0);
        let (theta, omega) = (y[0], y[1]);
        y[0] = omega;
        y[1] = -b * omega - c * f64::sin(theta);
        return y;
    };

    let t0 = 0.0;
    let y0 = vec![std::f64::consts::PI - 0.1, 0.0];
    let endpoint = 10.0;

    let (ts, ys) = driver.run(pendulum, t0, y0, endpoint, 0.1);
    for (t, y) in zip(ts, ys) {
        println!("{t}\t{}\t{}", y[0], y[1]);
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

    let t0 = 0.0;
    let z0 = vec![10.0, 15.0];
    let endpoint = 15.0;

    let (ts, zs) = driver.run(lotka_volterra, t0, z0, endpoint, 0.1);
    for (t, z) in zip(ts, zs) {
        println!("{t}\t{}\t{}", z[0], z[1]);
    }
    println!("\n");

    // Three body system
    let three_body = |_t: f64, mut y: Vec<f64>| -> Vec<f64> {
        /* ordering of y: [
            x1, y1, x2, y2, x3, y3, 
            vx1, vy1, vx2, vy2, vx3, vy3
        ] */
        let m = [1.0, 1.0, 1.0];
        let g = 1.0;
        let mut a = [0.0; 6];
        for i in 0..3 {
            for j in 0..3 {
                if j==i {continue;}
                let dx = y[2*j] - y[2*i];
                let dy = y[2*j+1] - y[2*i+1];
                let factor = g * m[j] / (dx*dx + dy*dy).powf(1.5);
                a[2*i] += dx * factor;
                a[2*i+1] += dy * factor;
            }
        }
        for i in 0..6 {
            y[i] = y[i+6];
            y[i+6] = a[i];
        }
        return y;
    };

    let v3 = [-0.93240737, -0.86473146];
    let x0 = 0.0;
    let y0 = vec![
        0.97000436, -0.24308753, -0.97000436, 0.24308753, 0.0, 0.0,  // r0s
        -v3[0]/2.0, -v3[1]/2.0, -v3[0]/2.0, -v3[1]/2.0, v3[0], v3[1],  // v0s
    ];
    let period = 6.32591398;
    let (_xs, ys) = driver.run(three_body, x0, y0, period, 0.05);

    for y in ys {
        println!("{}\t{}\t{}\t{}\t{}\t{}", y[0], y[1], y[2], y[3], y[4], y[5]);
    }
    println!("\n");

}