extern crate matrix;
extern crate scientific;
extern crate sfuns;

use matrix::Matrix;
use scientific::rand::Rng;
use sfuns::are_close;

use std::iter::zip;

fn qr_gs_decomp(mat: &mut Matrix<f64>, r: &mut Matrix<f64>) {
    let m = mat.num_cols;
    for i in 0..m {
        let ai = &mat[i];
        let norm = f64::sqrt(ai.iter().map(|x| x*x).sum());
        r[i][i] = norm;
        let qi: Vec<f64> = ai.iter().map(|x| x/norm).collect();
        mat[i].clone_from_slice(&qi[..]);
        for j in i+1..m {
            let aj = &mat[j];
            let inner_prod: f64 = qi.iter()
                                    .zip(aj.iter())
                                    .map(|x| x.0 * x.1)
                                    .sum();
            r[j][i] = inner_prod;
            let new_col: Vec<f64> = zip(aj.iter(), qi.iter()).map(|(a, q)| a - inner_prod * q).collect();
            mat[j].clone_from_slice(&new_col[..]);
        }
    }
}

fn back_substitution(a: &Matrix<f64>, b: &mut Matrix<f64>) {
    for i in (0..b.num_rows).rev() {
        let mut sum = 0.0;
        for j in i+1..b.num_rows {
            sum += a[j][i] * b[0][j];
        }
        b[0][i] = (b[0][i] - sum) / a[i][i];
    }
}

fn qr_gs_solve(q: &Matrix<f64>, r: &Matrix<f64>, b: &mut Matrix<f64>) {
    let temp = b.clone();
    b[0].clone_from_slice(&(q.transpose() * temp)[0]);
    back_substitution(r, b);
}

fn qr_gs_inverse(q: &Matrix<f64>, r: &Matrix<f64>) -> Matrix<f64> {
    let n = q.num_cols;
    let mut result = Matrix::zeros(n, n);
    for i in 0..n {
        let mut ei = Matrix::new(vec![vec![0.0; n]]);
        ei[0][i] = 1.0;
        qr_gs_solve(q, r, &mut ei);
        result[i].clone_from_slice(&ei[0]);
    }
    return result;
}

fn random_matrix(rows: usize, cols: usize, rng: &Rng) -> Matrix<f64> {
    return Matrix::from_data((0..rows*cols).map(|_| rng.f64()).collect(), rows, cols)
}

fn is_upper_tiangular(mat: &Matrix<f64>) -> bool {
    for i in 1..mat.num_rows {
        for j in 0..i {
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

    let args: Vec<String> = std::env::args().collect();
    for arg in &args[1..] {  // ignore first argument since this will always be ./main.bin
        let words: Vec<&str> = arg.split(":").collect();
        match words[0] {
            "-size" => size = words[1].parse::<usize>().unwrap(),
            "-test" => test = true,
            _ => panic!("Command '{}' not found.", words[0]),
        }
    }
    if test {test_functions(); return}
    if size == 0 {panic!("No size given. Use: -size:<u32>!")};
    solve(size);
}

fn solve(size: usize) {
    let rng = Rng::new(2347);
    let mut a = random_matrix(size, size, &rng);
    let mut r = Matrix::zeros(size, size);
    qr_gs_decomp(&mut a, &mut r);
}

fn test_functions() {
    let rng = Rng::new(2347);

    let (n, m) = (5, 3);
    let mut a = random_matrix(n, m, &rng) * 10.0;
    let mut r = Matrix::zeros(m, m);
    let a_copy = a.clone();
    qr_gs_decomp(&mut a, &mut r);

    println!("Generated random {n} by {m} matrix A");
    println!("{a_copy}\n");
    println!("Perfomed QR factorisation on A");
    println!("R is upper triangular?: {}", is_upper_tiangular(&r));
    println!("Q^TQ = I: {}", mat_are_close(&(&a.transpose() * &a), &Matrix::idty(m)));
    println!("QR = A: {}\n", mat_are_close(&(&a * &r), &a_copy));

    let mut a = random_matrix(n, n, &rng) * 10.0;
    let mut r = Matrix::zeros(n, n);
    let mut b = random_matrix(n, 1, &rng) * 10.0;
    let b_copy = b.clone();
    let a_copy = a.clone();
    qr_gs_decomp(&mut a, &mut r);
    qr_gs_solve(&a, &r, &mut b);

    println!("Generated random {n} by {n} matrix A");
    println!("{a_copy}\n");
    println!("and {n}-vector b\n{b_copy}\n");
    println!("Perfomed QR factorisation on A");
    println!("Solved Rx = Q^Tb for x by back-substitution");
    println!("QRx = b: {}\n", mat_are_close(&(&a_copy * &b), &b_copy));
    
    let a_inv = qr_gs_inverse(&a, &r);
    println!("Computed matrix inverse of A from above");
    println!("A^(-1)A = I: {}", mat_are_close(&(&a_copy * &a_inv), &Matrix::idty(n)));
}