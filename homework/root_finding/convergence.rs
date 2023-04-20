extern crate scientific;
extern crate matrix;
extern crate sfuns;
use matrix::linalg::optimisation::newton_root;
use scientific::ode::*;
use sfuns::linspace;

fn diff_eq(r: f64, mut y: Vec<f64>, energy: f64) ->  Vec<f64> {
    let f = y[0];
    let fp = y[1];
    let fpp = -2.0 * (1.0/r + energy) * f;
    y[0] = fp;
    y[1] = fpp;
    return y;
}

fn solve_schodinger(energy: f64, r_max: f64, r_min: f64, driver: &AdaptiveStepSizeDriver) -> (Vec<f64>, Vec<Vec<f64>>) {
    let y0 = vec![r_min-r_min*r_min, 1.0-2.0*r_min]; 
    return driver.run(|r: f64, y: Vec<f64>| diff_eq(r, y, energy), r_min, y0, r_max, 0.01)
}

fn run_convergence_calculations(r_maxs: &Vec<f64>, r_mins: &Vec<f64>, abss: &Vec<f64>, rels: &Vec<f64>) {
    let energy_start_guess = -1.0;
    let max_iter = 1e6 as u32;
    let acc = 1e-6;

    for abs in abss { for rel in rels {
        let driver = AdaptiveStepSizeDriver::new(RungeKuttaStepper::rk45(), Some((*abs, *rel)));
        for r_max in r_maxs {for r_min in r_mins{
            let objective = |x: &Vec<f64>| -> Vec<f64> {
                let energy = x[0];
                let (_, ys) = solve_schodinger(energy, *r_max, *r_min, &driver);
                return vec![ys.iter().map(|y| y[0]).last().unwrap()]
            };
            let eigen_energy = newton_root(&objective, vec![energy_start_guess], Some((max_iter, acc))).unwrap()[0];
            println!("{eigen_energy:+.8}\t{r_max:.1}\t{r_min:.4}\t{abs:.2e}\t{rel:e}");
        }}
    }}
}

fn main() {
    println!("# energy\tr_max\tr_min\tabs\trel");
    let r_max = vec![10.0];
    let r_min = vec![0.001];
    let (abs, rel) = (vec![1e-4], vec![0.0]);

    let mut r_maxs = linspace(7.0, 10.0, 10);
    r_maxs.push(10.0);
    let mut r_mins = linspace(0.1, 0.001, 10);
    r_mins.push(0.001);
    let mut exponents = linspace(-2.0, -4.0, 20);
    exponents.push(-4.0);
    let abss = exponents.into_iter().map(|x| 10.0_f64.powf(x)).collect();

    dbg!("Running r_min convergence calculations");
    run_convergence_calculations(&r_max, &r_mins, &abs, &rel);
    println!("\n");
    
    dbg!("Running r_max convergence calculations");
    run_convergence_calculations(&r_maxs, &r_min, &abs, &rel);
    println!("\n");

    dbg!("Running abs convergence calculations");
    run_convergence_calculations(&r_max, &r_min, &abss, &rel);
    println!("\n");
}