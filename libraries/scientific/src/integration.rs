use super::rand::Rng;

pub fn integrate(f: &impl Fn(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>) -> (f64, f64, u32) {
    match (a.is_infinite(), b.is_infinite()) {
        (true, true) => {
            if a == b {panic!("Bad integral limits")};
            let ff = |t: f64| (f)(t/(1.0-t*t)) * (1.0+t*t)/(1.0-t*t).powi(2);
            return recursive_adaptive_integrator(&ff, -1.0, 1.0, precision);
        },
        (false, true) => {
            if b < 0.0 {
                let (int, error, evals) = integrate(f, b, a, precision); 
                return (-int, error, evals);
            };
            let ff = |t: f64| (f)(a + (1.0-t)/t) / (t*t);
            return recursive_adaptive_integrator(&ff, 0.0, 1.0, precision);
        },
        (true, false) => {
            if a > 0.0 {
                let (int, error, evals) = integrate(f, b, a, precision); 
                return (-int, error, evals);
            };
            let ff = |t: f64| (f)(b - (1.0-t)/t) / (t*t);
            return recursive_adaptive_integrator(&ff, 0.0, 1.0, precision);
        },
        (false, false) => {
            return recursive_adaptive_integrator(f, a, b, precision);
        },
    };
}

pub fn recursive_adaptive_integrator(f: &impl Fn(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>) -> (f64, f64, u32) {
    let (abs, rel) = precision.unwrap_or((0.01, 0.01));
    let (f2, f3) = ((f)(a+2.0/6.0*(b-a)), (f)(a+4.0/6.0*(b-a)));
    let (int, variance, evals) = integrate_recursive(&f, a, b, abs, rel, f2, f3);
    return (int, variance.sqrt(), evals+2)
}

fn integrate_recursive(f: &impl Fn(f64) -> f64, a: f64, b: f64, abs:f64, rel:f64, f2: f64, f3: f64) -> (f64, f64, u32) {
    let (f1, f4) = ((f)(a+1.0/6.0*(b-a)), (f)(a+5.0/6.0*(b-a)));
    let int = (2.0*f1+f2+f3+2.0*f4) / 6.0 * (b-a); // higher order rule: trapezeum
    let int2 = (    f1+f2+f3+    f4) / 4.0 * (b-a); // lower order rule: rectangular
    let error = (int - int2).abs();
    if error <= f64::max(abs, rel*int.abs()) {
        return (int, error.powi(2), 2);
    } else {
        let (int1, variance1, evals1) = integrate_recursive(f, a, (a+b)/2.0, abs/f64::sqrt(2.0), rel, f1, f2);
        let (int2, variance2, evals2) = integrate_recursive(f, (a+b)/2.0, b, abs/f64::sqrt(2.0), rel, f3, f4);
        return (int1+int2, variance1+variance2, evals1+evals2+2);
    }
}

pub fn clenshaw_curtis(f: &impl Fn(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>) -> (f64, f64, u32) {
    let ff = |theta: f64| {
        (f)((a+b)/2.0 + (b-a)/2.0*f64::cos(theta)) * f64::sin(theta) * (b-a)/2.0
    };
    return recursive_adaptive_integrator(&ff, 0.0, std::f64::consts::PI, precision);
}

pub fn plain_monte_carlo(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, num_points: u64, rng: &mut Rng) -> (f64, f64) {
    let dim = a.len();
    let volume = std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));
    let (mut sum, mut sum2) = (0.0, 0.0);
    let mut x = vec![0.0; dim];
    for _ in 0..num_points {
        for k in 0..dim {x[k] = a[k] + rng.f64()*(b[k] - a[k])}
        let fx = f(&x);
        sum += fx;
        sum2 += fx * fx;
    }
    let mean = sum / num_points as f64;
    let sigma = (sum2/num_points as f64 - mean*mean).sqrt();
    
    return (mean * volume, sigma * volume / f64::sqrt(num_points as f64))
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    
    #[test]
    fn test_integrate_1() {
        let f = |x: f64| 4.0 * f64::sqrt(1.0 - x*x);
        let i = integrate(&f, 0.0, 1.0, Some((0.01, 0.0))).0;
        assert!((i - PI).abs() < 0.01);
    }

    #[test]
    fn test_integrate_2() {
        let f = |x: f64| 4.0 * f64::sqrt(1.0 - x*x);
        let i = integrate(&f, 1.0, 0.0, Some((0.01, 0.0))).0;
        assert!((i - (-PI)).abs() < 0.01);
    }

    #[test]
    fn test_integrate_3() {
        let f = |x: f64| f64::ln(x)/f64::sqrt(x);
        let i = integrate(&f, 0.0, 1.0, None).0;
        assert!((i - (-4.0)).abs() < 0.01);
    }

    #[test]
    fn test_clenshaw_curtis_1() {
        let f = |x: f64| 1.0/f64::sqrt(x);
        let i = clenshaw_curtis(&f, 0.0, 1.0, None).0;
        assert!((i - 2.0).abs() < 0.01);
    }

    #[test]
    fn test_clenshaw_curtis_2() {
        let f = |x: f64| f64::ln(x)/f64::sqrt(x);
        let i = clenshaw_curtis(&f, 0.0, 1.0, Some((0.01, 0.0))).0;
        assert!((i - (-4.0)).abs() < 0.01);
    }

    #[test]
    fn test_improper_1() {
        let f = |x: f64| x.powi(-2);
        let i = integrate(&f, 1.0, f64::INFINITY, None).0;
        assert!((i - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_improper_2() {
        let f = |x: f64| f64::exp(-x*x);
        let i = integrate(&f, f64::NEG_INFINITY, f64::INFINITY, None).0;
        assert!((i - PI.sqrt()).abs() < 0.01);
    }

    #[test]
    fn test_improper_3() {
        // 4! = gamma(5 - 1)
        let f = |x: f64| x.powi(5-1) * f64::exp(-x);
        let i = integrate(&f, 0.0, f64::INFINITY, None).0;
        assert!((i - 24.0).abs() < 0.01);
    }

    #[test]
    fn test_plain_monte_carlo_1() {
        let mut rng = Rng::new(1234);
        let f = |x: &Vec<f64>| {
            if x.iter().fold(0.0, |sum, val| sum + val*val).sqrt() <= 1.0 {
                return 1.0
            } else {
                return 0.0
            }
        };
        let (int, sigma) = plain_monte_carlo(&f, vec![0.0, 0.0], vec![2.0, 2.0], 1e6 as u64, &mut rng);
        assert!(
            (int - PI / 4.0).abs() < sigma
        );
    }

    #[test]
    fn test_plain_monte_carlo_2() {
        let mut rng = Rng::new(1234);
        let f = |x: &Vec<f64>| {
            if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
                return x[0] + x[1]
            } else {
                return 0.0
            }
        };
        let (int, sigma) = plain_monte_carlo(&f, vec![0.0, 0.0], vec![2.0, 2.0], 1e6 as u64, &mut rng);
        assert!(
            (int - 13.0 / 20.0).abs() < sigma
        );
    }
}