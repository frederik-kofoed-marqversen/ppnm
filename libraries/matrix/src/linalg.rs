use super::Matrix;

pub mod eig;
pub mod lstsq;
pub mod qr;
pub mod optimisation;

pub fn back_substitution(r: &Matrix<f64>, b: &mut Matrix<f64>) {
    for i in (0..b.num_rows).rev() {
        let mut sum = 0.0;
        for j in i+1..b.num_rows {
            sum += r[j][i] * b[0][j];
        }
        b[0][i] = (b[0][i] - sum) / r[i][i];
    }
}

/* pub fn gauss_elemination(a: &mut Matrix<f64>) {
    pseudo code at: https://en.wikipedia.org/wiki/Gaussian_elimination
} */

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::zip;
    
    #[test]
    fn test_back_substitution() {
        let a = Matrix::new(vec![vec![1.0, 0.0, 0.0], vec![-2.0, 1.0, 0.0], vec![1.0, 6.0, 1.0]]);
        let mut b = Matrix::new(vec![vec![4.0, -1.0, 2.0]]);
        back_substitution(&a, &mut b);
        
        assert!(
            zip(
                b.iter(), 
                [-24.0, -13.0, 2.0]
            ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-15))
        );
    }
}