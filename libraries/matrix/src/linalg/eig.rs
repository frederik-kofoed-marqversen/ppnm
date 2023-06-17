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

pub fn jacobi_cyclic(a: &mut Matrix<f64>) -> Matrix<f64> {
    // makes A diagonal -> D
    // returns transformation matrix such that A = V D V^T
    let n = a.num_cols;
    let mut v = Matrix::<f64>::idty(n);
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
                    times_j(&mut v, p, q, theta);
                    changed = true;
                }
            }
        }
    }
    return v
}

pub fn jacobi_cyclic_optimised(a: &mut Matrix<f64>) -> (Vec<f64>, Matrix<f64>) {
    // preserves lower triangle and diagonal of A 
    // At the end, upper triangle is zero
    let n = a.num_cols;
    let mut v = Matrix::<f64>::idty(n);
    let mut eigenvalues = Vec::with_capacity(n);
    for i in 0..n {eigenvalues.push(a[i][i])}
    
    let mut changed = true;
    while changed {
        changed = false;
        for p in 0..n-1 {
            for q in p+1..n {
                let (apq, app, aqq) = (a[q][p], eigenvalues[p], eigenvalues[q]);
                let theta = 0.5 * f64::atan2(2.0 * apq, aqq - app);
                let (c, s) = (theta.cos(), theta.sin());
                let new_app = c*c*app-2.0*s*c*apq+s*s*aqq;
                let new_aqq = s*s*app+2.0*s*c*apq+c*c*aqq;

                if !are_close(new_app, app) || !are_close(new_aqq, aqq) {
                    changed = true;  // do one more cycle after this one

                    for i in 0..p {  // update both collumn p and q up to i=p-1
                        let (aip, aiq) = (a[p][i], a[q][i]);
                        a[p][i] = c * aip - s * aiq;
                        a[q][i] = s * aip + c * aiq;
                    }
                    for j in q+1..n {  // update both row p and q from j=q+1
                        let (apj, aqj) = (a[j][p], a[j][q]);
                        a[j][p] = c * apj - s * aqj;
                        a[j][q] = s * apj + c * aqj;
                    }
                    for k in p+1..q {  // update the in-betweens
                        let (apk, akq) = (a[k][p], a[q][k]);
                        a[k][p] = c * apk - s * akq;
                        a[q][k] = c * akq + s * apk;
                    }
                    a[q][p] = 0.0;  // the zeroed out element
                    eigenvalues[p] = new_app;  // the diagonal values
                    eigenvalues[q] = new_aqq;
                    
                    times_j(&mut v, p, q, theta); // update v
                }
            }
        }
    }
    return (eigenvalues, v)
}

pub fn hessenberg(a: &mut Matrix<f64>) -> Matrix<f64> {
    // puts A to upper-Hessenberg form A <- H
    // returns transformation matrix such that A = V H V^T
    let n = a.num_cols;
    let mut v = Matrix::<f64>::idty(n);
    
    for p in 1..n-1 {
        for q in p+1..n {
            let theta = f64::atan2(-a.get(q, p-1), a.get(p, p-1));
            
            times_j(a, p, q, theta);
            j_times(a, p, q, -theta);
            times_j(&mut v, p, q, theta);
        }
    }

    return v
}


pub fn determinant_upper_hessenberg(h: &Matrix<f64>) -> f64 {
    let n = h.num_rows;
    assert!(h.num_cols == n, "Matrix is not square");
    let indices: Vec<usize> = (0..n).collect();
    return rec_det_hessenberg(h, &indices, &indices)
}

fn rec_det_hessenberg(h: &Matrix<f64>, rows: &[usize], cols: &[usize]) -> f64 {
    if rows.len() == 1 {
        return h.get(rows[0], cols[0])
    }
    let det1 = rec_det_hessenberg(h, &rows[1..], &cols[1..]);
    let det2 = rec_det_hessenberg(h, &[&rows[0..1], &rows[2..]].concat(), &cols[1..]);
    let det = h.get(rows[0], cols[0]) * det1 - h.get(rows[1], cols[0]) * det2;
    return det
}


#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_jacobi_evd() {
        let data: Vec<f64> = vec![4.0, 2.0, 2.0, 1.0];
        let mut mat = Matrix::from_data(data, 2, 2);
        let mat1 = mat.clone();
        let v = jacobi_cyclic(&mut mat);
        assert!(
            are_close(mat[0][0], 0.0)
            &&
            are_close(mat[1][1], 5.0)
        );
        assert_eq!(
            &v * mat * v.transpose(),
            mat1
        );
    }

    #[test]
    fn test_optimised_evd() {
        let data: Vec<f64> = vec![4.0, 2.0, 2.0, 1.0];
        let mut mat = Matrix::from_data(data, 2, 2);
        let (eigvals, v) = jacobi_cyclic_optimised(&mut mat);
        assert_eq!(
            eigvals,
            vec![0.0, 5.0]
        );
        assert_eq!(
            mat,
            Matrix::from_data(vec![4.0, 2.0, 0.0, 1.0], 2, 2)
        );
        assert_eq!(
            &v * v.transpose(),
            Matrix::idty(2)
        );
    }

    #[test]
    fn test_hessenberg() {
        let data: Vec<f64> = (1..=16).map(|n| n as f64).collect();
        let a = Matrix::from_data(data, 4, 4).transpose();
        let mut h = a.clone();
        let q = hessenberg(&mut h);
        
        assert!(
            h.get(2, 0).abs() + h.get(3, 0).abs() + h.get(3, 1).abs() < 1e-14
        );
        assert!(
            (&a - (&q * &h * q.transpose())).iter().fold(true, |acc, item| acc && item.abs() < 1e-14)
        );
        assert!(determinant_upper_hessenberg(&h).abs() < 1e-13);
    }

    #[test]
    fn test_hessenberg_determinant() {
        assert_eq!(determinant_upper_hessenberg(&Matrix::idty(5)), 1.0);
        
        let mat = Matrix::from_data(vec![1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0], 3, 3);
        assert_eq!(determinant_upper_hessenberg(&mat), 2.0);
    }
}