extern crate sfuns;
use sfuns::map_range;
use sfuns::search::binary_search_leftmost as search;

struct LinearSpline {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl LinearSpline {
    pub fn new(x: Vec<f64>, y: Vec<f64>) -> Self {
        return Self{ x: x, y: y};
    }
    
    pub fn evaluate(&self, z: f64) -> f64 {
        let i = search(&self.x, z) - 1;
        let dx = self.x[i+1] - self.x[i];
        let dy = self.y[i+1] - self.y[i];
        return self.y[i] + dy / dx * (z - self.x[i]);
    }

    pub fn integral(&self, z: f64) -> f64 {
        let mut result: f64 = 0.0;
        let i = search(&self.x, z) - 1;
        
        for j in 0..i {
            result += (self.y[j+1] + self.y[j]) * (self.x[j+1] - self.x[j]) / 2.0;
        }
        
        let dx = self.x[i+1] - self.x[i];
        let dy = self.y[i+1] - self.y[i];
        let y_z = self.y[i] + dy / dx * (z - self.x[i]);
        result += (y_z + self.y[i]) * (z - self.x[i]) / 2.0;
        
        return result;
    }

    pub fn derivative(&self, z: f64) -> f64 {
        let i = search(&self.x, z) - 1;
        let dx = self.x[i+1] - self.x[i];
        let dy = self.y[i+1] - self.y[i];
        return dy / dx;
    }
}

fn main() {
    let range = (0.0, 2.0 * std::f64::consts::PI);

    let mut xs = Vec::new();
    let mut ys = Vec::new();

    let num_samples = 10;
    for i in 0..=num_samples {
        let x = map_range(i as f64, &(0.0, num_samples as f64), &range);
        ys.push(f64::sin(x));
        xs.push(x);
    }

    let spline = LinearSpline::new(xs, ys);

    let num_points = 200;
    for i in 1..num_points {
        let x = map_range(i as f64, &(0.0, num_points as f64), &range);
        println!("{x} {} {} {}", spline.evaluate(x), spline.derivative(x), spline.integral(x));
    }
}