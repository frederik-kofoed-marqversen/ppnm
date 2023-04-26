use super::{Matrix, back_substitution, qr::decomp};
use std::iter::zip;

const STEP_LIM: f64 = 0.000000014901161193847656;  // f64::EPSILON.sqrt()

#[inline]
fn norm(v: &Vec<f64>) -> f64 {
    v.iter().fold(0.0, |sum, vi| sum+vi*vi).sqrt()
}
#[inline]
fn dot(u: &Vec<f64>, v: &Vec<f64>) -> f64 {
    zip(u, v).fold(0.0, |sum, (a, b)| sum + a*b)
}

pub fn jacobian(f: &impl Fn(&Vec<f64>) -> Vec<f64>, x: &Vec<f64>) -> Matrix<f64> {
    let fx = f(x);
    let (n, m) = (fx.len(), x.len());
    let mut jac = Matrix::zeros(n, m);

    let mut xdx = x.clone();
    for k in 0..m {
        let dx_k = STEP_LIM * f64::max(x[k].abs(), 1.0);
        xdx[k] += dx_k; // xdx = x + dx*e_k
        let fxdx = f(&xdx);
        for i in 0..n {
            jac.set((fxdx[i] - fx[i])/dx_k, i, k);
        }
        xdx[k] = x[k];  // xdx = x
    }

    return jac
}

pub fn gradient(f: &impl Fn(&Vec<f64>) -> f64, x: &Vec<f64>) -> Matrix<f64> {
    let ff = |x: &Vec<f64>| vec![f(x)];
    jacobian(&ff, x).transpose()
}

pub fn newton_root(f: &impl Fn(&Vec<f64>) -> Vec<f64>, x0: Vec<f64>, options: Option<(u32, f64)>) -> Result<Vec<f64>, Vec<f64>> {
    let (max_iter, acc) = options.unwrap_or((1000, 1e-3));
    let lambda_min = 2.0_f64.powi(-13);
    let fx = f(&x0);
    let (n, m) = (fx.len(), x0.len());
    
    let mut x = Matrix::from_data(x0, m, 1);
    for _ in 0..max_iter {
        let fx = Matrix::from_data(f(x.data()), n, 1);
        let norm_fx = norm(fx.data());
        if norm_fx < acc {
            return Ok(x[0].to_vec())
        }
        
        let mut jac = jacobian(f, x.data());
        let mut r = Matrix::zeros(m, m);
        decomp(&mut jac, &mut r);
        let mut dx = jac.transpose() * (-fx);
        back_substitution(&r, &mut dx);
        
        let mut lambda = 1.0;
        while norm(&f((&x+lambda*&dx).data())) > (1.0-lambda/2.0)*norm_fx && lambda > lambda_min {
            lambda /= 2.0;
        }

        x += &dx * lambda;
    }

    // Convergence not reached within max_iter iterations
    return Err(x.data().to_vec())
}

pub fn quasi_newton_min(f: &impl Fn(&Vec<f64>) -> f64, x0: Vec<f64>, options: Option<(u32, f64)>) -> Result<(Vec<f64>, f64, u32), (Vec<f64>, f64, u32)> {
    let (max_iter, acc) = options.unwrap_or((1000, 1e-3));
    let lambda_min = 2.0_f64.powi(-13);

    let dim = x0.len();
    let mut b = Matrix::idty(dim);  // inverse Hessian matrix
    let mut fx = f(&x0);
    let mut df = gradient(f, &x0);  // gradient
    let mut x = Matrix::from_data(x0, dim, 1);
    let mut iter: u32 = 0;

    while norm(df.data()) > acc {  // convergence
        iter += 1;
        if iter > max_iter {
            return Err((x[0].to_vec(), fx, max_iter))
        }
        
        // do linesearch along dx
        let dx = -&b * &df;  // Newton step
        let mut lambda = 1.0;
        while f((&x + lambda*&dx).data()) > fx && lambda > lambda_min {
            lambda /= 2.0;
        }

        let step = lambda * &dx;
        x += step.clone();
        fx = f(x.data());
        let old_df = df;
        df = gradient(f, x.data());
        
        // update b-matrix
        if lambda < lambda_min {  // reset b
            b = Matrix::idty(dim);
            continue
        } else {  // Symmetric Broyden’s update
            let y = &df - &old_df;
            let sy = dot(step.data(), y.data());
            if sy.abs() < 1e-6 {  // dont update
                continue
            }
            let u = &step - &b * &y;
            let gamma = dot(u.data(), y.data()) / (2.0 * sy);
            let a = (u - gamma*&step) / sy;
            let db = &a*step.transpose() + &step*a.transpose();
            b += db;
        }
    }

    return Ok((x[0].to_vec(), fx, iter))
}

pub fn dowhill_simplex(f: &impl Fn(&Vec<f64>) -> f64, mut points: Vec<Vec<f64>>, options: Option<(u32, f64)>) -> Result<(Vec<f64>, f64, u32), (Vec<f64>, f64, u32)> {
    let (max_iter, acc) = options.unwrap_or((1000, 1e-3));
    let (alpha, gamma, rho, sigma) = (1.0, 2.0, 0.5, 0.5);
    
    let dim = points[0].len();
    assert!(points.len() == dim + 1);

    let mut simplex: Vec<(f64, Matrix<f64>)> = Vec::with_capacity(dim + 1);
    for _ in 0..dim+1 {
        let x = Matrix::from_data(points.pop().unwrap(), dim, 1);
        let fx = f(x.data());
        simplex.push((fx, x));
    };

    let shrink = |simplex: &mut Vec<(f64, Matrix<f64>)>| {
        for i in 1..dim+1 {
            let x = &simplex[0].1 + sigma * (&simplex[i].1 - &simplex[0].1);
            let fx = f(x.data());
            simplex[i] = (fx, x);
        };
    };

    let convergence = |simplex: &Vec<(f64, Matrix<f64>)>| -> bool {
        let mut avg = 0.0;
        for i in 0..dim+1 {
            for j in i..dim+1 {
                avg += norm((&simplex[i].1 - &simplex[j].1).data());
            }
        }
        avg /= (dim + 1) as f64;
        return avg < acc
    };

    for it in 0..max_iter {
        simplex.sort_by(|a, b| (a.0).partial_cmp(&b.0).unwrap());
        if convergence(&simplex) {
            return Ok((simplex[0].1.data().to_vec(), simplex[0].0, it))
        };

        let centroid: Matrix<f64> = simplex[0..dim].iter().fold(Matrix::zeros(dim, 1), |sum, x| sum + &x.1) / dim as f64;
        let reflected = &centroid + alpha * (&centroid - &simplex[dim].1);
        let f_re = f(reflected.data());
        
        if simplex[0].0 <= f_re && f_re < simplex[dim-1].0 {
            simplex[dim] = (f_re, reflected);
            continue;
        };
        if f_re < simplex[0].0 {
            let expanded = &centroid + gamma * (&reflected - &centroid);
            let f_ex = f(expanded.data());
            if f_ex < f_re {
                simplex[dim] = (f_ex, expanded);
                continue;
            } else {
                simplex[dim] = (f_re, reflected);
                continue;
            }
        };
        if f_re < simplex[dim].0 {
            let contracted = &centroid + rho * (&reflected - &centroid);
            let f_co = f(contracted.data());
            if f_co < f_re {
                simplex[dim] = (f_co, contracted);
                continue
            } else {
                shrink(&mut simplex);
                continue;
            };
        } else {
            let contracted = &centroid + rho * (&simplex[dim].1 - &centroid);
            let f_co = f(contracted.data());
            if f_co < simplex[dim].0 {
                simplex[dim] = (f_co, contracted);
                continue
            } else {
                shrink(&mut simplex);
                continue;
            };
        };
    };
    return Err((simplex[0].1.data().to_vec(), simplex[0].0, max_iter))
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
    fn test_newton_root() {
        let f = |x: &Vec<f64>| -> Vec<f64> {
            vec![
                -2.0*(1.0-x[0]) - 400.0*x[0]*(x[1]-x[0]*x[0]), 
                200.0*(x[1]-x[0]*x[0])
            ]
        };
        let x = vec![0.0, 0.0];
        let result = newton_root(&f, x, Some((1000, 1e-6))).unwrap();
        
        assert!(zip(
            result, 
            vec![1.0, 1.0]
        ).fold(true, |acc, (item, test)| acc && ((item-test).abs() < 1e-6)));
    }

    #[test]
    #[should_panic]
    fn test_newton_root_error() {
        let f = |x: &Vec<f64>| -> Vec<f64> {vec![x[0]*x[0] + 1.0]};
        let x = vec![1.0];
        newton_root(&f, x, Some((1000, 1e-6))).unwrap();
    }

    #[test]
    fn test_min_1() {
        // Rosenbrock's valley function
        let f = |x: &Vec<f64>| -> f64 {(1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0]*x[0]).powi(2)};
        
        let x0 = vec![0.0, 2.0];
        let acc = 0.01;
        let max_iter = 10000;
        let (x, fx, iter) = quasi_newton_min(&f, x0, Some((max_iter, acc))).unwrap();
        
        let x_min = vec![1.0, 1.0];
        assert!(
            (x[0] - x_min[0]).abs() < acc &&
            (x[1] - x_min[1]).abs() < acc
        );
        assert!(
            (fx - 0.0).abs() < acc
        );
        assert!(
            iter < max_iter
        );
    }

    #[test]
    fn test_min_2() {
        // Himmelblau's function
        let f = |x: &Vec<f64>| -> f64 {(x[0]*x[0] + x[1] - 11.0).powi(2) + (x[0] + x[1]*x[1] - 7.0).powi(2)};
        
        let x0 = vec![0.0, 0.0];
        let acc = 0.01;
        let max_iter = 10000;
        let (x, fx, iter) = quasi_newton_min(&f, x0, Some((max_iter, acc))).unwrap();
        
        let x_min = vec![3.0, 2.0];
        assert!(
            (x[0] - x_min[0]).abs() < acc &&
            (x[1] - x_min[1]).abs() < acc
        );
        assert!(
            (fx - 0.0).abs() < acc
        );
        assert!(
            iter < max_iter
        );
    }

    #[test]
    fn test_downhill_simplex() {
        // Rosenbrock's valley function
        let f = |x: &Vec<f64>| -> f64 {(1.0 - x[0]).powi(2) + 100.0 * (x[1] - x[0]*x[0]).powi(2)};
        
        let x0 = vec![2.0, 2.0];
        let x1 = vec![2.0, 2.1];
        let x2 = vec![2.1, 2.1];

        let acc = 1e-6;
        let max_iter = 10000;
        let (x, fx, iter) = dowhill_simplex(&f, vec![x0, x1, x2], Some((max_iter, acc))).unwrap();
        
        dbg!(&x, &fx, &iter);

        let x_min = vec![1.0, 1.0];
        assert!(
            (x[0] - x_min[0]).abs() < acc &&
            (x[1] - x_min[1]).abs() < acc
        );
        assert!(
            (fx - 0.0).abs() < acc
        );
        assert!(
            iter < max_iter
        );
    }
}