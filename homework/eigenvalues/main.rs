extern crate matrix;

use matrix::Matrix;
use matrix::linalg::eig::jacobi_cyclic;
use std::iter::zip;

fn main() {
    let mut num_states: usize = 0;
    let mut r_max = 0.0;
    let mut dr = 0.0;
    
    let mut args: Vec<String> = std::env::args().collect();
    args.remove(0); // ignore first argument since this will always be ./main.bin
    while let Some(arg) = args.pop() { 
        let words: Vec<&str> = arg.split(":").collect();
        match words[0] {
            "-r_max" => r_max = words[1].parse::<f64>().unwrap(),
            "-dr" => dr = words[1].parse::<f64>().unwrap(),
            "-N" => num_states = words[1].parse::<usize>().unwrap(),
            _ => panic!("Command '{}' not found.", words[0]),
        }
    }
    
    if r_max == 0.0 || dr == 0.0 {panic!("Both of r_max and dr must be given!");}
    solve_hydrogen(r_max, dr, num_states);
}

fn solve_hydrogen(r_max: f64, dr: f64, num_states: usize) {
    let num_points = (r_max / dr) as usize - 1;
    let r: Vec<f64> = (0..num_points).map(|i| dr*(i+1) as f64).collect();
    let mut h = Matrix::<f64>::zeros(num_points, num_points);
    let mut v = Matrix::<f64>::idty(num_points);
    for i in 0..num_points-1 {
        h[i][i] = -2.0;
        h[i+1][i] = 1.0;
        h[i][i+1] = 1.0;
    }
    h[num_points-1][num_points-1] = -2.0;
    let factor = -0.5/dr/dr;
    for elem in h.iter_mut() {
        *elem *= factor
    } 
    for i in 0..num_points {h[i][i] += -1.0/r[i];}

    jacobi_cyclic(&mut h, &mut v);

    if num_states == 0 {
        println!("{:.4} {r_max} {dr}", h[0][0]); // lowest energy
    } else {
        for i in 0..num_states {
            // println!("state {i} energy: {:.4}", h[i][i]);
            println!("0  0");
            for (rj, vij) in zip(r.iter(), v[i].iter()) {
                println!("{rj}  {vij}");
            }
            println!("{r_max}  0");
            println!("\n");
        }
    }
}