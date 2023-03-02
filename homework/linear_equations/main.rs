extern crate matrix;
extern crate rand;
extern crate sfuns;

use matrix::Matrix;
use rand::Rng;
use sfuns::are_close;

fn qr_gs_decomp(mat: &mut Matrix, r: &mut Matrix) {
    let m = mat.num_cols;
    for i in 0..m {
        let ai = &mat.data[i];
        let norm = f64::sqrt(ai.iter().map(|x| x*x).sum());
        r.data[i][i] = norm;
        let qi: Vec<f64> = ai.iter().map(|x| x/norm).collect();
        mat.data[i] = qi.clone();
        for j in i+1..m {
            let aj = &mat.data[j];
            let inner_prod: f64 = qi.iter()
                                    .zip(aj.iter())
                                    .map(|x| x.0 * x.1)
                                    .sum();
            r.data[j][i] = inner_prod;
            mat.data[j] = aj.iter().zip(qi.iter()).map(|(a, q)| a - inner_prod * q).collect();
        }
    }
}

fn back_substitution(a: &Matrix, b: &mut Matrix) {
    for i in (0..b.num_rows).rev() {
        let mut sum = 0.0;
        for j in i+1..b.num_rows {
            sum += a.data[j][i] * b.data[0][j];
        }
        b.data[0][i] = (b.data[0][i] - sum) / a.data[i][i];
    }
}

fn qr_gs_solve(q: &Matrix, r: &Matrix, b: &mut Matrix) {
    b.data = (&q.transpose() * &b).data;
    back_substitution(r, b);
}

fn qr_gs_inverse(q: &Matrix, r: &Matrix) -> Matrix {
    let n = q.num_cols;
    let mut result = Matrix::zeros(n, n);
    for i in 0..n {
        let mut ei = Matrix::new(vec![vec![0.0; n]]);
        ei.data[0][i] = 1.0;
        qr_gs_solve(q, r, &mut ei);
        result.data[i] = ei.data[0].clone();
    }
    return result;
}

fn random_vec(len: usize, rng: &Rng) -> Vec<f64> {
    let mut result = vec![];
    result.resize_with(len, || rng.f64() * 10.0);
    return result;
}
fn random_matrix(rows: usize, cols: usize, rng: &Rng) -> Matrix {
    let mut data = vec![];
    data.resize_with(cols, || random_vec(rows, rng));
    return Matrix::new(data)
}

fn is_upper_tiangular(mat: &Matrix) -> bool {
    for (i, col) in mat.data.iter().enumerate() {
        for elem in &col[i+1..] {
            if !are_close(*elem, 0.0) {return false};
        }
    }
    return true;
}

fn mat_are_close(a: &Matrix, b: &Matrix) -> bool {
    for (a_col, b_col) in a.data.iter().zip(b.data.iter()) {
        for (a_elem, b_elem) in a_col.iter().zip(b_col.iter()) {
            if !are_close(*a_elem, *b_elem) {return false};
        }
    }
    return true;
}

fn main() {
    let rng = Rng::new(2347);

    let (n, m) = (5, 3);
    let mut a = random_matrix(n, m, &rng);
    let mut r = Matrix::zeros(m, m);
    let a_copy = a.clone();
    qr_gs_decomp(&mut a, &mut r);

    println!("Generated random {n} by {m} matrix A");
    println!("Perfomed QR factorisation on A");
    println!("R is upper triangular?: {}", is_upper_tiangular(&r));
    println!("Q^TQ = I: {}", mat_are_close(&(&a.transpose() * &a), &Matrix::idty(m)));
    println!("QR = A: {}\n", mat_are_close(&(&a * &r), &a_copy));

    let mut a = random_matrix(n, n, &rng);
    let mut r = Matrix::zeros(n, n);
    let mut b = random_matrix(n, 1, &rng);
    let b_copy = b.clone();
    let a_copy = a.clone();
    qr_gs_decomp(&mut a, &mut r);
    qr_gs_solve(&a, &r, &mut b);

    println!("Generated new random {n} by {n} matrix A and {n}-vector b");
    println!("Perfomed QR factorisation on A");
    println!("Solved Rx = Q^Tb for x by back-substitution");
    println!("QRx = b: {}\n", mat_are_close(&(&a_copy * &b), &b_copy));
    
    let a_inv = qr_gs_inverse(&a, &r);
    println!("Computed matrix inverse of A from above");
    println!("A^(-1)A = I: {}", mat_are_close(&(&a_copy * &a_inv), &Matrix::idty(n)));
}