extern crate matrix;
extern crate rand;
extern crate sfuns;

use matrix::Matrix;
use rand::Rng;
use sfuns::are_close;

use std::iter::zip;

fn times_j(a: &mut Matrix<f64>, p: usize, q: usize, theta: f64) {
    let (c, s) = (theta.cos(), theta.sin());
    for i in 0..a.num_rows {
        let (aip, aiq) = (a[p][i], a[q][i]);
        a[p][i] = c * aip - s * aiq;
        a[q][i] = s * aip + c * aiq;
    }
}

fn j_times(a: &mut Matrix<f64>, p: usize, q: usize, theta: f64) {
    let (c, s) = (theta.cos(), theta.sin());
    for j in 0..a.num_cols {
        let (apj, aqj) = (a[j][p], a[j][q]);
        a[p][j] = c * apj + s * aqj;
        a[q][j] = -s * apj + c * aqj;
    }
}

fn jacobi_cyclic(a: &mut Matrix<f64>, v: &mut Matrix<f64>) {
    let n = a.num_cols;
    let mut changed = false;
    while changed {
        changed = false;
        for p in 0..n-1 {
            for q in p+1..n {
                let (apq, app, aqq) = (a[q][p], a[p][p], a[q][q]);
                let theta = 0.5 * f64::atan2(2.0 * apq, aqq - app);
                let (c, s) = (theta.cos(), theta.sin());
                let new_app = c*c*app-2.0*s*c*apq+s*s*aqq;
                let new_aqq = s*s*app+2.0*s*c*apq+c*c*aqq;
                if new_app!=app || new_aqq!=aqq {
                    times_j(a, p, q, theta);
                    j_times(a, p, q, -theta);
                    times_j(v, p, q, theta);
                    changed = true;
                }
            }
        }
    }
}

fn random_matrix(rows: usize, cols: usize, rng: &Rng) -> Matrix<f64> {
    return Matrix::from_data((0..rows*cols).map(|_| rng.f64()).collect(), rows, cols)
}

fn mat_are_close(mat1: &Matrix<f64>, mat2: &Matrix<f64>) -> bool {
    for (a, b) in zip(mat1.iter(), mat2.iter()) {
        if !are_close(*a, *b) {return false};
    }
    return true;
}

fn main() {
    let rng = Rng::new(2347);

    let n = 15;
    let mut a = random_matrix(n, n, &rng);
    let mut v = Matrix::idty(n);
    let a_copy = a.clone();
    jacobi_cyclic(&mut a, &mut v);
    
    println!("Generated random {n} by {n} matrix A");
    println!("Perfomed cyclic Jacobi eigenvalue decomposition on A");
    println!("V^TV = I: {}", mat_are_close(&(&v.transpose() * &v), &Matrix::idty(n)));
    println!("VV^T = I: {}", mat_are_close(&(&v * &v.transpose()), &Matrix::idty(n)));
    println!("V D V^T = A: {}", mat_are_close(&(&(&v * &a) * &v.transpose()), &a_copy));
    println!("V^T A V = D: {}\n", mat_are_close(&(&(&v.transpose() * &a_copy) * &v), &a));

    let (r_max, dr) = (10.0, 0.3);
    let num_points = (r_max / dr) as usize - 1;
    let r: Vec<f64> = (0..num_points).map(|i| dr*(i+1) as f64).collect();
    let mut h = Matrix::zeros(num_points, num_points);
    for i in 0..num_points-1 {
        h[i][i] = -2.0;
        h[i+1][i] = 1.0;
        h[i][i+1] = 1.0;
    }
    h[num_points-1][num_points-1] = -2.0;
    let factor = -0.5/dr/dr;
    for elem in h.iter_mut() {
        *elem *= factor
    } 
    for i in 0..num_points {h[i][i] += -1.0/r[i];}
}