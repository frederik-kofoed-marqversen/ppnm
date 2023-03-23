extern crate sfuns;
extern crate matrix;
extern crate rand;
use sfuns::are_close;
use matrix::Matrix;
use rand::Rng;

fn solve_tridiagonal_system(left_diag: &Vec<f64>, diag: &mut Vec<f64>, right_diag: &Vec<f64>, b: &mut Vec<f64>) {
    // solve (left_diag)_i x_{i-1} + (diag)_i x_i + (right_diag)_i x_{i+1} = b_i
    // all vectors must have same length n
    // a[0] and q[n-1] are simply ignored!
    
    // simplified Gauss elemination for tri-diagonal systems
    let n = left_diag.len();
    for i in 1..n {
        let w = left_diag[i] / diag[i-1];
        diag[i] -= w * right_diag[i-1];
        b[i] -= w * b[i-1];
    }
    // back substitution for bi-diagonal system
    b[n-1] = b[n-1] / diag[n-1];
    for i in (0..n-1).rev() {
        b[i] = (b[i] - right_diag[i] * b[i+1]) / diag[i];
    }
}

fn mat_are_close(mat1: &Matrix<f64>, mat2: &Matrix<f64>) -> bool {
    for (a, b) in std::iter::zip(mat1.iter(), mat2.iter()) {
        if !are_close(*a, *b) {return false};
    }
    return true;
}

fn main() {
    let rng = Rng::new(1234);
    
    let n = 4;

    let mut b = vec![];
    let mut left_diag = vec![];
    let mut diag = vec![];
    let mut right_diag = vec![];
    for _ in 0..n {
        b.push(rng.f64() * 10.0);
        left_diag.push(rng.f64() * 10.0);
        diag.push(rng.f64() * 10.0);
        right_diag.push(rng.f64() * 10.0);
    }
    
    left_diag[0] = 0.0;
    right_diag[n-1] = 0.0;

    let mut mat = Matrix::<f64>::zeros(n, n);
    let b_copy = Matrix::from_data(b.clone(), n, 1);

    mat[0][0] = diag[0];
    mat[1][0] = right_diag[0];
    for i in 1..n-1 {
        mat[i][i] = diag[i];
        mat[i+1][i] = right_diag[i];
        mat[i-1][i] = left_diag[i];
    }
    mat[n-1][n-1] = diag[n-1];
    mat[n-2][n-1] = left_diag[n-1];

    solve_tridiagonal_system(&left_diag, &mut diag, &right_diag, &mut b);
    let b = Matrix::from_data(b, n, 1);

    println!("Setting up random {n}-dimensional tridiagonal system Ax=b");
    println!("A = {}", &mat);
    println!("b = {}\n", &b_copy);
    println!("Solving the system in-place in O(n)");
    println!("x = {}", &b);
    println!("Ax == b? {}", mat_are_close(&(mat * b), &b_copy));
}