fn solve_tridiagonal_system(left_diag: &Vec<f64>, diag: &mut Vec<f64>, right_diag: &Vec<f64>, b: &mut Vec<f64>) {
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

fn binary_search_bin<T: PartialOrd>(sorted_array: &Vec<T>, item: T) -> usize {
    if item < sorted_array[0] || item > *sorted_array.last().unwrap() {panic!("item is out of bounds")}

    let (mut left, mut right) = (0, sorted_array.len() - 1);
    while right - left > 1 {
        let midpoint: usize = (left + right) / 2;
        if item > sorted_array[midpoint] {
            left = midpoint;
        } else {
            right = midpoint;
        }
    }
    return left;
}

#[derive(Debug)]
pub struct CubicSpline {
    x: Vec<f64>,
    y: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    d: Vec<f64>
}

impl CubicSpline {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        let (b, c, d) = Self::compute_coefficients(&x, &y);
        Self{x: x, y: y, b: b, c: c, d: d}
    }

    fn compute_coefficients(x: &Vec<f64>, y: &Vec<f64>) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let n = x.len();
        let mut b = Vec::<f64>::with_capacity(x.len());
        let mut d = Vec::<f64>::with_capacity(x.len());
        let mut q = Vec::<f64>::with_capacity(x.len());

        let mut dx = Vec::<f64>::with_capacity(x.len());
        let mut dydx = Vec::<f64>::with_capacity(x.len());
        dx.push(x[1] - x[0]);
        dydx.push((y[1] - y[0]) / dx[0]);
        
        d.push(2.0);
        q.push(1.0);
        b.push(3.0 * dydx[0]);
        for i in 1..n-1 {
            dx.push(x[i+1] - x[i]);
            dydx.push((y[i+1] - y[i]) / dx[i]);

            d.push(2.0 * (dx[i-1]/dx[i] + 1.0));
            q.push(dx[i-1]/dx[i]);
            b.push(3.0 * (dydx[i-1] + dydx[i] * dx[i-1]/dx[i]));
        }
        d.push(2.0);
        //q.push(0.0);
        b.push(3.0 * dydx[n-2]);

        let left_diag = vec![1.0; n]; // the left diagonal
        solve_tridiagonal_system(&left_diag, &mut d, &q, &mut b); // b is now the right spline coefficient

        d.pop();
        for i in 0..n-1 {
            // reuse q-vector as spline c coefficient
            q[i] = (-2.0*b[i] - b[i+1] + 3.0 * dydx[i]) / dx[i];
            // reuse d-vector as spline d coefficient
            d[i] = (b[i] + b[i+1] - 2.0 * dydx[i]) / (dx[i] * dx[i]);
        }

        return (b, q, d)
    }
    
    pub fn evaluate(&self, z: f64) -> f64 {
        let i = binary_search_bin(&self.x, z);
        let dx = z - self.x[i];
        return self.y[i] 
             + self.b[i] * dx 
             + self.c[i] * dx * dx
             + self.d[i] * dx * dx * dx
    }

    pub fn derivative(&self, z: f64) -> f64 {
        let i = binary_search_bin(&self.x, z);
        let dx = z - self.x[i];
        return self.b[i] 
             + self.c[i] * dx * 2.0
             + self.d[i] * dx * dx * 3.0
    }

    pub fn integral(&self, z: f64) -> f64 {
        let bin_index = binary_search_bin(&self.x, z);
        let mut result = 0.0;
        for (i, x) in self.x[1..=bin_index].iter().chain(std::iter::once(&z)).enumerate() {
            let dx = x - self.x[i];
            result += self.y[i] * dx 
                    + self.b[i] * dx * dx / 2.0 
                    + self.c[i] * dx * dx * dx / 3.0 
                    + self.d[i] * dx * dx * dx * dx / 4.0
        }
        return result;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_binary_search_1() {
        let array: Vec<f64> = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_bin(&array, 2.0),
            0
        );
    }

    #[test]
    fn test_binary_search_2() {
        let array: Vec<f64> = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_bin(&array, 2.1),
            0
        );
    }

    #[test]
    fn test_binary_search_3() {
        let array: Vec<f64> = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(
            binary_search_bin(&array, 7.0),
            4
        );
    }

    #[test]
    #[should_panic]
    fn test_binary_search_4() {
        let array: Vec<f64> = vec![2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
        binary_search_bin(&array, 10.0);
    }
}