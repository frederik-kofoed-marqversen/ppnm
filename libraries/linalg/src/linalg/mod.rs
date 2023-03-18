use super::Matrix;

pub mod eig;
pub mod lstsq;
pub mod qr;

pub fn back_substitution(r: &Matrix<f64>, b: &mut Matrix<f64>) {
    for i in (0..b.num_rows).rev() {
        let mut sum = 0.0;
        for j in i+1..b.num_rows {
            sum += r[j][i] * b[0][j];
        }
        b[0][i] = (b[0][i] - sum) / r[i][i];
    }
}

pub fn solve_tridiagonal_system(left_diag: &Vec<f64>, diag: &mut Vec<f64>, right_diag: &Vec<f64>, b: &mut Vec<f64>) {
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

/* pub fn gauss_elemination(a: &mut Matrix<f64>) {
    pseudo code at: https://en.wikipedia.org/wiki/Gaussian_elimination
} */