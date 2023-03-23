extern crate sfuns;
extern crate scientific;
use scientific::spline::CubicSpline;
use sfuns::map_range;

fn main() {
    let gamma = 1.0;
    let x0 = 2.0;
    let range = (x0 - 10.0, x0 + 10.0);

    let mut xs = Vec::new();
    let mut ys = Vec::new();

    let lorentz = |x: f64| -> f64 {(std::f64::consts::PI * gamma * (1.0 + ((x - x0) / gamma).powf(2.0))).powf(-1.0)};

    let num_samples = 20;
    for i in 0..num_samples {
        let x = map_range(i as f64, &(0.0, (num_samples-1) as f64), &range);
        let y = lorentz(x.clone());

        println!("{x} {y}");
        xs.push(x);
        ys.push(y);
    }
    println!("\n");

    let spline = CubicSpline::new(xs, ys);

    let num_points = 200;
    for i in 0..num_points {
        let x = map_range(i as f64, &(0.0, num_points as f64), &range);
        println!("{x} {} {} {} {}", lorentz(x.clone()), spline.evaluate(x), spline.derivative(x), spline.integral(x));
    }
}