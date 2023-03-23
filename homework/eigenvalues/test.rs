extern crate matrix;
extern crate rand;
extern crate sfuns;

use matrix::Matrix;
use matrix::linalg::eig::jacobi_cyclic;
use rand::Rng;
use sfuns::are_close;

use std::iter::zip;

fn random_symmetric_matrix(n: usize, rng: &Rng) -> Matrix<f64> {
    let mut mat = Matrix::zeros(n, n);
    for j in 0..n {
        mat[j][j] = rng.f64() * 10.0;
        for i in 0..j {
            let num = rng.f64() * 10.0;
            mat[j][i] = num;
            mat[i][j] = num;
        }
    }
    return mat
}

fn mat_are_close(mat1: &Matrix<f64>, mat2: &Matrix<f64>) -> bool {
    for (a, b) in zip(mat1.iter(), mat2.iter()) {
        if !are_close(*a, *b) {return false};
    }
    return true;
}

fn main() {
    let mut print: bool = false;
    let mut n: usize = 5;
    
    let mut args: Vec<String> = std::env::args().collect();
    args.remove(0); // ignore first argument since this will always be ./main.bin
    while let Some(arg) = args.pop() { 
        let words: Vec<&str> = arg.split(":").collect();
        match words[0] {
            "-print" => print = words[1].parse::<bool>().unwrap(),
            "-size" => n = words[1].parse::<usize>().unwrap(),
            _ => panic!("Command '{}' not found.", words[0]),
        }
    }

    let rng = Rng::new(1234);
    let mut a = random_symmetric_matrix(n, &rng);
    let mut v = Matrix::idty(n);

    if !print {
        jacobi_cyclic(&mut a, &mut v);
    } else {
        let a_copy = a.clone();
        jacobi_cyclic(&mut a, &mut v);

        println!("Generated random symmetric {n} by {n} matrix A");
        println!("A={a_copy}\n");
        println!("Perfomed cyclic Jacobi eigenvalue decomposition on A=V D V^T");
        println!("D={a}\n");
        println!("V={v}\n");
        println!("V^TV = I: {}", mat_are_close(&(v.transpose() * &v), &Matrix::idty(n)));
        println!("VV^T = I: {}", mat_are_close(&(&v * v.transpose()), &Matrix::idty(n)));
        println!("V D V^T = A: {}", mat_are_close(&(&v * &a * v.transpose()), &a_copy));
        println!("V^T A V = D: {}\n", mat_are_close(&(v.transpose() * &a_copy * &v), &a));
    }
}