use super::Matrix;

fn are_close(a: f64, b: f64) -> bool {
    let acc = 1e-9;
    let eps = 1e-9;
    if f64::abs(a - b) < acc {return true}
    if f64::abs(a - b) < f64::max(a.abs(), b.abs()) * eps {return true}
    return false;
}

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
        a[j][p] = c * apj + s * aqj;
        a[j][q] = -s * apj + c * aqj;
    }
}

pub fn jacobi_cyclic(a: &mut Matrix<f64>, v: &mut Matrix<f64>) {
    let n = a.num_cols;
    let mut changed = true;
    while changed {
        changed = false;
        for p in 0..n-1 {
            for q in p+1..n {
                let (apq, app, aqq) = (a[q][p], a[p][p], a[q][q]);
                let theta = 0.5 * f64::atan2(2.0 * apq, aqq - app);
                let (c, s) = (theta.cos(), theta.sin());
                let new_app = c*c*app-2.0*s*c*apq+s*s*aqq;
                let new_aqq = s*s*app+2.0*s*c*apq+c*c*aqq;
                if !are_close(new_app, app) || !are_close(new_aqq, aqq) {
                    times_j(a, p, q, theta);
                    j_times(a, p, q, -theta);
                    times_j(v, p, q, theta);
                    changed = true;
                }
            }
        }
    }
}