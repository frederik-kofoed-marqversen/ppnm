use super::{Matrix, back_substitution};
use std::iter::zip;

pub fn decomp(mat: &mut Matrix<f64>, r: &mut Matrix<f64>) {
    let m = mat.num_cols;
    for i in 0..m {
        let ai = &mat[i];
        let norm = f64::sqrt(ai.iter().map(|x| x*x).sum());
        r[i][i] = norm;
        let qi: Vec<f64> = ai.iter().map(|x| x/norm).collect();
        mat[i].clone_from_slice(&qi[..]);
        for j in i+1..m {
            let aj = &mat[j];
            let inner_prod: f64 = zip(qi.iter(), aj.iter()).map(|x| x.0 * x.1).sum();
            r[j][i] = inner_prod;
            let new_col: Vec<f64> = zip(aj.iter(), qi.iter()).map(|(a, q)| a - inner_prod * q).collect();
            mat[j].clone_from_slice(&new_col[..]);
        }
    }
}

pub fn inverse(q: &Matrix<f64>, r: &Matrix<f64>) -> Matrix<f64> {
    let n = q.num_cols;
    let mut result = Matrix::zeros(n, n);
    for i in 0..n {
        let mut b = q.row(i);
        back_substitution(r, &mut b);
        result[i].clone_from_slice(&b[0]);
    }
    return result;
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::iter::zip;
    
    #[test]
    fn test_decomp() {
        let mut a = Matrix::new(vec![vec![2.0, 2.0, 1.0], vec![3.0, 4.0, 1.0]]);
        let mut r = Matrix::zeros(2, 2);
        decomp(&mut a, &mut r);
        
        assert!(
            zip(
                r.iter(), 
                [3.0, 0.0, 5.0, 1.0]
            ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-15))
        );
        assert!(
            zip(
                a.iter(), 
                [2.0/3.0, 2.0/3.0, 1.0/3.0, -1.0/3.0, 2.0/3.0, -2.0/3.0]
            ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-15))
        );
    }

    #[test]
    fn test_inverse() {
        let mut a = Matrix::new(vec![vec![1.0, 0.0, 3.0], vec![3.0, 2.0, 1.0], vec![1.0, 2.0, 0.0]]);
        let mut r = Matrix::zeros(3, 3);
        decomp(&mut a, &mut r);
        let a_inv = inverse(&a, &r);
        
        assert!(
            zip(
                a_inv.iter(), 
                [-0.2, 0.6, -0.6, 0.1, -0.3, 0.8, 0.4, -0.2, 0.2]
            ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-15))
        );
    }
}