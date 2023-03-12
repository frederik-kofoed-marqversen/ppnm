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