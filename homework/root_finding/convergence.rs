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

fn run_convergence_calculations_part_c(energy_start_guess: f64, r_maxs: &Vec<f64>, r_mins: &Vec<f64>, abss: &Vec<f64>, rels: &Vec<f64>) {
    let max_iter = 1e3 as u32;
    let acc = 1e-3;

    for abs in abss { for rel in rels {
        let driver = AdaptiveStepSizeDriver::new(RungeKuttaStepper::rk45(), Some((*abs, *rel)));
        for r_max in r_maxs {for r_min in r_mins{
            let objective = |x: &Vec<f64>| -> Vec<f64> {
                let energy = x[0];
                let (rs, ys) = solve_schodinger(energy, *r_max, *r_min, &driver);

                let (r_end, y_end) = (rs.last().unwrap(), ys.last().unwrap());
                let (f, df) = (y_end[0], y_end[1]);
                let k = f64::sqrt(2.0 * energy.abs());
                
                return vec![df * r_end - f * (1.0 - k*r_end)];
            };
            let eigen_energy = match newton_root(&objective, vec![energy_start_guess], Some((max_iter, acc))) {
                Ok(energy) => energy[0],
                Err(energy) => energy[0],
            };
            println!("{eigen_energy:+.8}\t{r_max:.1}\t{r_min:.4}\t{abs:.2e}\t{rel:e}");
        }}
    }}
}

fn main() {
    let r_max = vec![12.0];
    let r_min = vec![0.001];
    let (abs, rel) = (vec![1e-4], vec![0.0]);

    let mut r_maxs = linspace(2.0, 12.0, 10);
    r_maxs.push(12.0);
    let mut r_mins = linspace(0.08, 0.001, 10);
    r_mins.push(0.001);
    let mut exponents = linspace(-2.5, -4.0, 20);
    exponents.push(-4.0);
    let abss = exponents.into_iter().map(|x| 10.0_f64.powf(x)).collect();




    // PART B CONVERGENCE CALCULATIONS

    println!("# PART B CONVERGENCE CALCULATIONS");
    println!("# energy\tr_max\tr_min\tabs\trel");

    dbg!("Running r_min convergence calculations");
    run_convergence_calculations(&r_max, &r_mins, &abs, &rel);
    println!("\n");
    
    dbg!("Running r_max convergence calculations");
    run_convergence_calculations(&r_maxs, &r_min, &abs, &rel);
    println!("\n");

    dbg!("Running abs convergence calculations");
    run_convergence_calculations(&r_max, &r_min, &abss, &rel);
    println!("\n");





    // PART C CONVERGENCE CALCULATIONS
    println!("# PART C CONVERGENCE CALCULATIONS");
    println!("# energy\tr_max\tr_min\tabs\trel");
    
    dbg!("Part C convergence calculations");
    let mut r_maxs = linspace(5.0, 25.0, 15);
    r_maxs.push(25.0);
    
    dbg!("Running for ground state");
    run_convergence_calculations_part_c(-0.55, &r_maxs, &r_min, &abs, &rel);
    println!("\n");

    let mut r_maxs = linspace(5.0, 30.0, 15);
    r_maxs.push(30.0);

    dbg!("Running for 1st excited state");
    run_convergence_calculations_part_c(-0.12, &r_maxs, &r_min, &abs, &rel);
    println!("\n");

    dbg!("Running for 2nd excited state");
    run_convergence_calculations_part_c(-0.05, &r_maxs, &r_min, &abs, &rel);
    println!("\n");
}