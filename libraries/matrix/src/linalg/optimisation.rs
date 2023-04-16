use super::{Matrix, back_substitution, qr::decomp};

const STEP_LIM: f64 = 0.000000014901161193847656;  // f64::EPSILON.sqrt()

#[inline]
fn norm(v: &Vec<f64>) -> f64 {
    v.iter().fold(0.0, |sum, vi| sum+vi*vi).sqrt()
}

pub fn jacobian(f: &impl Fn(&Vec<f64>) -> Vec<f64>, x: &Vec<f64>) -> Matrix<f64> {
    let fx = f(x);
    let (n, m) = (fx.len(), x.len());
    let mut result = Matrix::zeros(n, m);

    let mut xdx = x.clone();
    for k in 0..m {
        let dx_k = STEP_LIM * f64::max(x[k].abs(), 1.0);
        xdx[k] += dx_k; // xdx = x + dx*e_k
        let fxdx = f(&xdx);
        for i in 0..n {
            result.set((fxdx[i] - fx[i])/dx_k, i, k);
        }
        xdx[k] = x[k];  // xdx = x
    }

    return result
}

pub fn newton(f: &impl Fn(&Vec<f64>) -> Vec<f64>, x: Vec<f64>, tol: Option<f64>) -> Vec<f64> {
    let tol: f64 = tol.unwrap_or(1e-3);
    let fx = f(&x);
    let (n, m) = (fx.len(), x.len());
    
    let mut x = Matrix::from_data(x, m, 1);
    loop {
        let fx = Matrix::from_data(f(x.data()), n, 1);
        let norm_fx = norm(fx.data());
        if norm_fx < tol {
            break
        }
        
        let mut jac = jacobian(f, x.data());
        let mut r = Matrix::zeros(m, m);
        decomp(&mut jac, &mut r);
        let mut dx = jac.transpose() * (-fx);
        back_substitution(&r, &mut dx);
        
        let mut lambda = 1.0;
        while norm(&f((&x+lambda*&dx).data())) > (1.0-lambda/2.0)*norm_fx && lambda > 1.0/32.0 {
            lambda /= 2.0;
        }

        x += &dx * lambda;
    }
    return x.data().clone()
}

#[cfg(test)]
mod tests{
    use std::iter::zip;
    use super::*;

    #[test]
    fn test_jacobian() {
        let f = |x: &Vec<f64>| -> Vec<f64> {vec![x[0]*x[0]*x[2], 5.0*x[0]]};
        let x = vec![3.0, 2.0, 1.0];
        assert!(zip(
            jacobian(&f, &x).iter(), 
            [6.0, 5.0, 0.0, 0.0, 9.0, 0.0]
        ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-6)));
    }

    #[test]
    fn test_newton() {
        let f = |x: &Vec<f64>| -> Vec<f64> {
            vec![
                -2.0*(1.0-x[0]) - 400.0*x[0]*(x[1]-x[0]*x[0]), 
                200.0*(x[1]-x[0]*x[0])
            ]
        };
        let x = vec![0.0, 0.0];
        
        assert!(zip(
            newton(&f, x, Some(1e-6)), 
            vec![1.0, 1.0]
        ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-6)));
    }
}