extern crate matrix;
extern crate scientific;
extern crate sfuns;

use matrix::Matrix;
use matrix::linalg::eig::hessenberg as hessenberg_decomp;
use matrix::linalg::eig::determinant_upper_hessenberg as determinant;
use matrix::linalg::qr::decomp as qr_decomp;
use scientific::rand::Rng;
use sfuns::are_close;

use std::iter::zip;

fn random_matrix(rows: usize, cols: usize, rng: &Rng) -> Matrix<f64> {
    return Matrix::from_data((0..rows*cols).map(|_| rng.f64()).collect(), rows, cols)
}

fn is_upper_hessenberg(mat: &Matrix<f64>) -> bool {
    assert!(mat.num_cols == mat.num_rows, "Matrix is not square");
    for j in 0..mat.num_cols {
        for i in j+2..mat.num_rows {
            if !are_close(mat.get(i, j), 0.0) {return false};
        }
    }
    return true;
}

fn mat_are_close(mat1: &Matrix<f64>, mat2: &Matrix<f64>) -> bool {
    for (a, b) in zip(mat1.iter(), mat2.iter()) {
        if !are_close(*a, *b) {return false};
    }
    return true;
}

fn main() {
    let mut size: usize = 0;
    let mut test: bool = false;
    let mut method: &str = "Hessenberg";

    let args: Vec<String> = std::env::args().collect();
    for arg in &args[1..] {  // ignore first argument since this will always be ./main.bin
        let words: Vec<&str> = arg.split(":").collect();
        match words[0] {
            "-size" => size = words[1].parse::<usize>().unwrap(),
            "-method" => method = words[1],
            "-test" => test = true,
            _ => panic!("Command '{}' not found.", words[0]),
        }
    }
    if test {test_functions(); return}
    if size == 0 {panic!("No size given. Use: -size:<u32>")};
    match method {
        "hessenberg" => hessenberg(size),
        "qr" => qr(size),
        _ => panic!("Unknown method '{}'", method),
    }
}

fn hessenberg(size: usize) {
    let rng = Rng::new(2347);
    let mut a = random_matrix(size, size, &rng);
    // the hessenberg decomposition alg allocate memory for q internally
    let _q = hessenberg_decomp(&mut a);
}

fn qr(size: usize) {
    let rng = Rng::new(2347);
    let mut a = random_matrix(size, size, &rng);
    let mut r = Matrix::zeros(size, size);
    // the QR factorisation alg does not allocate memory for r
    qr_decomp(&mut a, &mut r);
}

fn test_functions() {
    let rng = Rng::new(2347);

    let n = 5;
    let a = random_matrix(n, n, &rng) * 10.0;
    let mut h = a.clone();
    let q = hessenberg_decomp(&mut h);

    println!("Generated random {n} by {n} matrix A");
    println!("A = {a}\n");
    println!("Computed Hessenberg decomposition: A = Q H Q^T");
    println!("H = {h}\n");
    println!("Q = {q}\n");
    println!("H is in Hessenberg form?: {}", is_upper_hessenberg(&h));
    println!("Q is orthonormal (Q^TQ = QQ^T = I)?: {}",
        mat_are_close(&(q.transpose() * &q), &Matrix::idty(n))
        &&
        mat_are_close(&(&q * q.transpose()), &Matrix::idty(n))
    );
    println!("Q H Q^T = A?: {}\n", mat_are_close(&(&q * &h * q.transpose()), &a));
    println!("det(A) = det(H) = {}", determinant(&h));
}