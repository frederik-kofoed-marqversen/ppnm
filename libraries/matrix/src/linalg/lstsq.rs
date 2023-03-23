use super::{qr, back_substitution};
use super::Matrix;

pub fn fit(fns: &Vec<&dyn Fn(f64) -> f64>, x: &Matrix<f64>, y: &Matrix<f64>, dy: &Matrix<f64>) -> (Matrix<f64>, Matrix<f64>) {
    let mut a = Matrix::zeros(x.num_rows, fns.len());
    let mut r = Matrix::idty(a.num_cols);
    let mut b = Matrix::zeros(x.num_rows, 1);
    for k in 0..fns.len() { // build a matrix A
        let fk = fns[k];
        let a_k = &mut a[k];
        for i in 0..x.num_rows {
            a_k[i] = fk(x[0][i]) / dy[0][i]
        }
    }
    for i in 0..x.num_rows { // build b vector
        b[0][i] = y[0][i] / dy[0][i];
    }
    
    qr::decomp(&mut a, &mut r); // QR decompose A
    b = a.transpose() * b;
    back_substitution(&r, &mut b); // solve linear equation R c = Q^T b
    r = qr::inverse(&Matrix::idty(r.num_rows), &r); // inverse of R
    let sigma = &r * r.transpose(); // compute covariance matrix
    return (b, sigma)
}

pub fn correlations(covariance: &mut Matrix<f64>) { // only updates upper triangular part
    covariance[0][0] = f64::sqrt(covariance[0][0]);
    for j in 1..covariance.num_cols {
        covariance[j][j] = f64::sqrt(covariance[j][j]);
        for i in (0..j).rev() {
            covariance[j][i] /= covariance[j][j] * covariance[i][i];
        }
    }
}