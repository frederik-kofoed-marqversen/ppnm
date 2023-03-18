extern crate sfuns;
use sfuns::search::binary_search_leftmost;
use std::iter;

#[derive(Debug)]
struct QuadraticSpline {
    x: Vec<f64>,
    y: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
}

impl QuadraticSpline {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        let p = Self::compute_p(&x, &y);
        let c = Self::forward_recursion(0.0, &x, &p);
        let c = Self::backward_recursion(c.last().unwrap() / 2.0, &x, &p);
        let b = p.iter().enumerate().map(|(i, pi)| pi - c[i] * (x[i+1] - x[i])).collect();
        Self{x: x, y: y, b: b, c: c}
    }

    fn compute_p(x: &Vec<f64>, y: &Vec<f64>) -> Vec<f64> {
        let mut p = Vec::with_capacity(x.len()-1);
        for i in 0..x.len()-1 {
            p.push((y[i+1] - y[i]) / (x[i+1] - x[i]));
        }
        return p;
    }

    fn forward_recursion(c_0: f64, x: &Vec<f64>, p: &Vec<f64>) -> Vec<f64> {
        let mut c = Vec::with_capacity(x.len()-1);
        c.push(c_0);
        for i in 0..x.len()-2 {
            let dp = p[i+1] - p[i];
            let dx_i = x[i+1] - x[i];
            let dx_ip1 = x[i+2] - x[i+1];
            c.push((dp - c[i]*dx_i) / dx_ip1);
        }
        return c;
    }

    fn backward_recursion(c_nm1: f64, x: &Vec<f64>, p: &Vec<f64>) -> Vec<f64> {
        let n = x.len();
        let mut c = vec![0.0; n-1];
        c[n-2] = c_nm1;
        for i in (0..x.len()-2).rev() {
            let dp = p[i+1] - p[i];
            let dx_i = x[i+1] - x[i];
            let dx_ip1 = x[i+2] - x[i+1];
            c[i] = (dp - c[i+1]*dx_ip1) / dx_i;
        }
        return c;
    }
    
    pub fn evaluate(&self, z: f64) -> f64 {
        let i = binary_search_leftmost(&self.x, z) - 1;
        let dx = z - self.x[i];
        return self.y[i] 
             + self.b[i] * dx 
             + self.c[i] * dx * dx;
    }

    pub fn derivative(&self, z: f64) -> f64 {
        let i = binary_search_leftmost(&self.x, z) - 1;
        let dx = z - self.x[i];
        return self.b[i] 
             + self.c[i] * dx * 2.0;
    }

    pub fn integral(&self, z: f64) -> f64 {
        let bin_index = binary_search_leftmost(&self.x, z) - 1;
        let mut result = 0.0;
        for (i, x) in self.x[1..=bin_index].iter().chain(iter::once(&z)).enumerate() {
            let dx = x - self.x[i];
            result += self.y[i] * dx 
                    + self.b[i] * dx * dx / 2.0 
                    + self.c[i] * dx * dx * dx / 3.0 
        }
        return result;
    }
}

fn main() {
    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = vec![0.0; x.len()];
    let b = vec![0.0; x.len()-1];
    let c = vec![0.0; x.len()-1];
    let spline = QuadraticSpline::new(x, y);
    println!("Spline 1:\nCorrect b? {}\nCorrect c? {}\n", b==spline.b, c==spline.c);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = x.clone();
    let b = vec![1.0; x.len()-1];
    let c = vec![0.0; x.len()-1];
    let spline = QuadraticSpline::new(x, y);
    println!("Spline 2:\nCorrect b? {}\nCorrect c? {}\n", b==spline.b, c==spline.c);

    let x = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let y = x.iter().map(|x| x * x).collect();
    let b = vec![2.0, 4.0, 6.0, 8.0];
    let c = vec![1.0; x.len()-1];
    let spline = QuadraticSpline::new(x, y);
    println!("Spline 3:\nCorrect b? {}\nCorrect c? {}\n", b==spline.b, c==spline.c);
}