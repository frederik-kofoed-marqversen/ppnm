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
}