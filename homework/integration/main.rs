extern crate sfuns;

use sfuns::map_range;
use std::f64::consts::PI;
use std::io::Write;

pub fn integrate(f: &mut impl FnMut(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>) -> f64 {
    let (abs, rel) = precision.unwrap_or((0.01, 0.01));
    let (f2, f3) = ((f)(a+2.0/6.0*(b-a)), (f)(a+4.0/6.0*(b-a)));
    return integrate_recursive(f, a, b, abs, rel, f2, f3);
}

fn integrate_recursive(f: &mut impl FnMut(f64) -> f64, a: f64, b: f64, abs:f64, rel:f64, f2: f64, f3: f64) -> f64 {
    let (f1, f4) = ((f)(a+1.0/6.0*(b-a)), (f)(a+5.0/6.0*(b-a)));
    let int = (2.0*f1+f2+f3+2.0*f4) / 6.0 * (b-a); // higher order rule: trapezeum
    let int2 = (    f1+f2+f3+    f4) / 4.0 * (b-a); // lower order rule: rectangular
    let error = (int - int2).abs();
    if error <= f64::max(abs, rel*int.abs()) {
        return int;
    } else {
        return integrate_recursive(f, a, (a+b)/2.0, abs/f64::sqrt(2.0), rel, f1, f2)
             + integrate_recursive(f, (a+b)/2.0, b, abs/f64::sqrt(2.0), rel, f3, f4);
    }
}

pub fn clenshaw_curtis(f: &mut impl FnMut(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>) -> f64 {
    let mut ff = |theta: f64| {
        (f)((a+b)/2.0 + (b-a)/2.0*f64::cos(theta)) * f64::sin(theta) * (b-a)/2.0
    };
    return integrate(&mut ff, 0.0, std::f64::consts::PI, precision);
}

pub fn integrate_improper(f: &mut impl FnMut(f64) -> f64, a: f64, b: f64, precision: Option<(f64, f64)>) -> f64 {
    match (a.is_infinite(), b.is_infinite()) {
        (true, true) => {
            let mut ff = |t: f64| (f)(t/(1.0-t*t)) * (1.0+t*t)/(1.0-t*t).powi(2);
            return integrate(&mut ff, -1.0, 1.0, precision);
        },
        (false, true) => {
            let mut ff = |t: f64| (f)(a + (1.0-t)/t) / (t*t);
            return integrate(&mut ff, 0.0, 1.0, precision);
        },
        (true, false) => {
            let mut ff = |t: f64| (f)(b - (1.0-t)/t) / (t*t);
            return integrate(&mut ff, 0.0, 1.0, precision);
        },
        (false, false) => return integrate(f, a, b, precision),
    };
}

fn erf(z: f64) -> f64 {
    integrate(&mut |x: f64| f64::exp(-x*x) * 2.0 / f64::sqrt(PI), 0.0, z, Some((0.0001, 0.0)))
}

fn main() -> std::io::Result<()> {
    // write erf.data
    let mut file = std::fs::File::create("erf.data")?;
    let num_points = 100;
    for i in 0..100 {
        let z = map_range(i as f64, &(0.0, num_points as f64), &(-3.0, 3.0));
        file.write(format!("{z}\t{}\n", erf(z)).as_bytes())?;
    }
    
    // count integrand evaluations
    let mut file = std::fs::File::create("Out.txt")?;

    // f(x)=1/√(x)
    file.write("Integrand evaluations for f(x)=1/√(x)\n".as_bytes())?;
    let mut count: usize = 0;
    let mut f = |x| {
        count += 1;
        return 1.0/f64::sqrt(x);
    };
    integrate(&mut f, 0.0, 1.0, None);
    let string = format!("\tintegrate(f, 0, 1): {}\n", count);
    file.write(string.as_bytes())?;
    
    let mut count: usize = 0;
    let mut f = |x| {
        count += 1;
        return 1.0/f64::sqrt(x);
    };
    clenshaw_curtis(&mut f, 0.0, 1.0, None);
    let string = format!("\tclenshaw_curtis(f, 0, 1): {}\n", count);
    file.write(string.as_bytes())?;

    // f(x)=ln(x)/√(x)
    file.write("\nIntegrand evaluations for f(x)=ln(x)/√(x)\n".as_bytes())?;
    let mut count: usize = 0;
    let mut f = |x| {
        count += 1;
        return f64::ln(x)/f64::sqrt(x);
    };
    integrate(&mut f, 0.0, 1.0, None);
    let string = format!("\tintegrate(f, 0, 1): {}\n", count);
    file.write(string.as_bytes())?;
    
    let mut count: usize = 0;
    let mut f = |x| {
        count += 1;
        return f64::ln(x)/f64::sqrt(x);
    };
    clenshaw_curtis(&mut f, 0.0, 1.0, None);
    let string = format!("\tclenshaw_curtis(f, 0, 1): {}\n", count);
    file.write(string.as_bytes())?;



    return Ok(());
}

#[cfg(test)]
mod tests {
    #![allow(non_snake_case)]
    use super::*;
    
    #[test]
    fn integral_A1() {
        let mut f = |x: f64| f64::sqrt(x);
        let i = integrate(&mut f, 0.0, 1.0, None);
        assert!((i - 2.0/3.0).abs() < 0.01);
    }

    #[test]
    fn integral_A2() {
        let mut f = |x: f64| 1.0/f64::sqrt(x);
        let i = integrate(&mut f, 0.0, 1.0, None);
        assert!((i - 2.0).abs() < 0.01);
    }
    
    #[test]
    fn integral_A3() {
        let mut f = |x: f64| 4.0 * f64::sqrt(1.0 - x*x);
        let i = integrate(&mut f, 0.0, 1.0, Some((0.01, 0.0)));
        assert!((i - std::f64::consts::PI).abs() < 0.01);
    }

    #[test]
    fn integral_A4() {
        let mut f = |x: f64| f64::ln(x)/f64::sqrt(x);
        let i = integrate(&mut f, 0.0, 1.0, None);
        assert!((i - (-4.0)).abs() < 0.01);
    }

    #[test]
    fn integral_B1() {
        let mut f = |x: f64| 1.0/f64::sqrt(x);
        let i = clenshaw_curtis(&mut f, 0.0, 1.0, None);
        assert!((i - 2.0).abs() < 0.01);
    }

    #[test]
    fn integral_B2() {
        let mut f = |x: f64| f64::ln(x)/f64::sqrt(x);
        let i = clenshaw_curtis(&mut f, 0.0, 1.0, Some((0.01, 0.0)));
        assert!((i - (-4.0)).abs() < 0.01);
    }

    #[test]
    fn improper_polynomial_integral() {
        // x**(-2)
        let mut f = |x: f64| x.powi(-2);
        let i = integrate_improper(&mut f, 1.0, f64::INFINITY, None);
        assert!((i - 1.0).abs() < 0.01);
    }

    #[test]
    fn gaussian_integral() {
        // e**(-x**2)
        let mut f = |x: f64| f64::exp(-x*x);
        let i = integrate_improper(&mut f, f64::NEG_INFINITY, f64::INFINITY, None);
        assert!((i - PI.sqrt()).abs() < 0.01);
    }

    #[test]
    fn gamma_function() {
        // 4! = gamma(5 - 1)
        let mut f = |x: f64| x.powi(5-1) * f64::exp(-x);
        let i = integrate_improper(&mut f, 0.0, f64::INFINITY, None);
        assert!((i - 24.0).abs() < 0.01);
    }
}