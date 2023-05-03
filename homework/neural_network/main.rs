extern crate matrix;
extern crate scientific;
extern crate sfuns;

use std::iter::zip;
use matrix::linalg::optimisation::dowhill_simplex;
use scientific::integration::integrate;
use std::cell::RefCell;
use sfuns::linspace;
use std::f64::consts::PI;

struct ArtificialNeuralNetwork {
    n: usize,
    activation_function: Box<dyn Fn(f64) -> f64>,
    derivative: Box<dyn Fn(f64) -> f64>,
    curvature: Box<dyn Fn(f64) -> f64>,
    antidirivative: Box<dyn Fn(f64) -> f64>,
    parameters: RefCell<Vec<f64>>,
}

#[allow(dead_code)]
impl ArtificialNeuralNetwork {
    fn new(hidden_neurons: usize, x_range: [f64; 2]) -> Self {
        let weights: Vec<f64> = (0..hidden_neurons).map(|n| (-1.0_f64).powi(n as i32)).collect();
        let mut shifts = linspace(x_range[0], x_range[1], hidden_neurons-1);
        shifts.push(x_range[1]);
        let scalings = vec![1.0; hidden_neurons];
        
        let mut parameters = Vec::new();
        parameters.extend(weights);
        parameters.extend(shifts);
        parameters.extend(scalings);

        Self {
            n: hidden_neurons,
            activation_function: Box::new(|x: f64| x * f64::exp(-x*x)),
            derivative: Box::new(|x: f64| (1.0 - 2.0*x*x) * f64::exp(-x*x)),
            curvature: Box::new(|x: f64| (2.0*x*x - 3.0) * 2.0 * x * f64::exp(-x*x)),
            antidirivative: Box::new(|x: f64| -0.5 * f64::exp(-x*x)),
            parameters: RefCell::new(parameters),
        }
    }

    fn response(&self, x: f64) -> f64 {
        let n = self.n;
        let parameters = self.parameters.borrow();
        (0..n).map(
            |i| (self.activation_function)((x-parameters[i+n])/parameters[i+2*n]) * parameters[i]
        ).sum()
    }

    fn derivative(&self, x: f64) -> f64 {
        let n = self.n;
        let parameters = self.parameters.borrow();
        (0..n).map(
            |i| (self.derivative)((x-parameters[i+n])/parameters[i+2*n]) * parameters[i] / parameters[i+2*n]
        ).sum()
    }

    fn curvature(&self, x: f64) -> f64 {
        let n = self.n;
        let parameters = self.parameters.borrow();
        (0..n).map(
            |i| (self.curvature)((x-parameters[i+n])/parameters[i+2*n]) * parameters[i] / parameters[i+2*n].powi(2)
        ).sum()
    }

    fn antidirivative(&self, x: f64) -> f64 {
        let n = self.n;
        let parameters = self.parameters.borrow();
        (0..n).map(
            |i| (self.antidirivative)((x-parameters[i+n])/parameters[i+2*n]) * parameters[i] * parameters[i+2*n]
        ).sum()
    }

    fn train(&mut self, cost_function: impl Fn(&Self) -> f64) -> f64 {
        let objective = |parameters: &Vec<f64>| -> f64 {
            *self.parameters.borrow_mut() = parameters.clone();
            cost_function(&self)
        };

        let mut simplex = vec![(*self.parameters.borrow()).clone()];
        for i in 0..simplex[0].len() {
            let mut point = simplex[0].clone();
            point[i] += point[i] / 2.0;
            simplex.push(point);
        }

        match dowhill_simplex(&objective, simplex, Some((10000, 1e-3))) {
            Err((_, result, _)) => {
                eprintln!("Optimisation did not converge => may result in inadequate network response!");
                return result
            }
            Ok((_, result, _)) => return result,
        };
    }
}

fn main() {
    // training data
    let f = |x: f64| f64::cos(5.0*x - 1.0) * f64::exp(-x*x);
    let range = [-1.0, 1.0];
    let num_points = 10;
    
    let xs: Vec<f64> = linspace(range[0], range[1], num_points);
    let ys: Vec<f64> = xs.iter().map(|x| f(*x)).collect();

    println!("# Part A - training data");
    for (x, y) in zip(&xs, &ys) {
        println!("{x} {y}");
    }
    println!("\n");

    // train network
    let mut network = ArtificialNeuralNetwork::new(4, range);
    let cost_function = |ann: &ArtificialNeuralNetwork| -> f64 {
        zip(&xs, &ys).map(|(x, y)| (ann.response(*x) - y).powi(2)).sum() 
    };
    network.train(cost_function);

    // print network response
    println!("# Part A");
    for x in linspace(range[0], range[1], 200) {
        println!("{x:.4} {:.6}", network.response(x));
    }
    println!("\n");




    // PART C
    let range = [-2.0*PI, 2.0*PI];
    let x0: f64 = 0.0;
    let y0: Vec<f64> = vec![1.0, 0.0];
    
    let differential_equation = |_r: f64, y: Vec<f64>| -> f64 {
        let (f, _, ddf) = (y[0], y[1], y[2]);
        return ddf + f
    };

    let mut network = ArtificialNeuralNetwork::new(6, range.clone());
    let cost_function = |ann: &ArtificialNeuralNetwork| -> f64 {
        let (alpha, beta) = (1.0, 1.0);
        
        let int = integrate(
            &|x| {
                let y = vec![ann.response(x), ann.derivative(x), ann.curvature(x)];
                return differential_equation(x, y).powi(2)
            },
            range[0], 
            range[1],
            Some((1e-4, 0.0))
        ).0;
        let p1 = alpha * (ann.response(x0) - y0[0]).powi(2);
        let p2 = beta * (ann.derivative(x0) - y0[1]).powi(2);

        return int + p1 + p2;
    };
    network.train(cost_function);

    // print network response
    println!("# Part C");
    for x in linspace(range[0], range[1], 200) {
        println!("{x:.4} {:.6}", network.response(x));
    }
    println!("\n");

}