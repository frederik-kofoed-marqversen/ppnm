extern crate matrix;
use matrix::linalg::optimisation::{quasi_newton_min, dowhill_simplex};
use std::io::BufRead;

fn parse_string(string_to_parse: &str, split_delimiters: Vec<char>) -> Vec<f64> {
    let mut result = Vec::new();
    for entry in string_to_parse.split(&split_delimiters[..]).filter(|&x| !x.is_empty()) {
        let number = entry.parse::<f64>().expect(&format!("Could not interpret '{entry}' as an f64"));
        result.push(number);
    }
    return result;
}

fn main() -> std::io::Result<()> {
    let acc = 0.01;
    let max_iter = 10000;

    // Rosenbrock's valley function
    println!("Rosenbrock's valley function");
    println!("Starting point around (0, 2)");
    let f = |x: &Vec<f64>| -> f64 {(1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0]*x[0]).powi(2)};
    println!("Quasi Newton");
    let x0 = vec![0.0, 2.0];
    let (x, fx, iter) = quasi_newton_min(&f, x0, Some((max_iter, acc))).unwrap();
    println!("\tFound minimum x = [{}, {}]\n\twith f(x)={fx}\n\twithin {iter} iterations.", x[0], x[1]);
    println!("Downhill simplex");
    let x0 = vec![0.0, 2.0];
    let x1 = vec![0.0, 2.1];
    let x2 = vec![0.1, 2.0];
    let (x, fx, iter) = dowhill_simplex(&f, vec![x0, x1, x2], Some((max_iter, acc))).unwrap();
    println!("\tFound minimum x = [{}, {}]\n\twith f(x)={fx}\n\twithin {iter} iterations.", x[0], x[1]);

    // Himmelblau's function
    println!("\nHimmelblau's function");
    println!("Starting point around (0, 0)");
    let f = |x: &Vec<f64>| -> f64 {(x[0]*x[0] + x[1] - 11.0).powi(2) + (x[0] + x[1]*x[1] - 7.0).powi(2)};
    println!("Quasi Newton");
    let x0 = vec![0.0, 0.0];
    let (x, fx, iter) = quasi_newton_min(&f, x0, Some((max_iter, acc))).unwrap();
    println!("\tFound minimum x = [{}, {}]\n\twith f(x)={fx}\n\twithin {iter} iterations.", x[0], x[1]);
    println!("Downhill simplex");
    let x0 = vec![0.0, 0.0];
    let x1 = vec![0.0, 0.1];
    let x2 = vec![0.1, 0.0];
    let (x, fx, iter) = dowhill_simplex(&f, vec![x0, x1, x2], Some((max_iter, acc))).unwrap();
    println!("\tFound minimum x = [{}, {}]\n\twith f(x)={fx}\n\twithin {iter} iterations.", x[0], x[1]);


    // read data from standard input stream
    let mut data: Vec<Vec<f64>> = Vec::new();
    let mut lines = std::io::stdin().lock().lines();
    while let Some(line) = lines.next() {
        let line = line?;
        if line.len() == 0 {break}
        if line.chars().nth(0).unwrap() == '#' {continue}
        data.push(parse_string(&line, vec![' ', '\t']));
    }

    let breit_wigner = |e: f64, variables: &Vec<f64>| -> f64 {
        let (m, gamma, a) = (variables[0], variables[1], variables[2]);
        a * gamma / std::f64::consts::PI / ((e - m).powi(2) + gamma*gamma/4.0)
    };
    
    let cost = |variables: &Vec<f64>| -> f64 {
        let mut result = 0.0;
        for point in data.iter() {
            let (e, sigma, delta) = (point[0], point[1], point[2]);
            result += (breit_wigner(e, variables) - sigma).powi(2) / (delta * delta);
        }
        return result
    };

    let variables_0 = vec![125.0, 5.0, 5.0];
    let (x, _, _) = quasi_newton_min(&cost, variables_0.clone(), Some((max_iter, acc))).unwrap();
    
    println!("\nFitted Breit-Wigner function to Higgs data");
    println!("f(E) = A * (Γ / π) / [(E-m)^2 + Γ^2/4]");
    println!("with: m={:.2}, Γ={:.2}, and A={:.2}", x[0], x[1], x[2]);

    println!("\n--------------------------------------------");
    println!("Fit function data");
    println!("--------------------------------------------");
    
    println!("\n");
    let f = |e: f64| breit_wigner(e, &x);
    let range: [f64; 2] = [100.0, 160.0];
    let num = 500;
    for n in 0..=num {
        let e = n as f64/num as f64 * (range[1] - range[0]) + range[0];
        println!("{e} {}", f(e as f64));
    }

    Ok(())
}