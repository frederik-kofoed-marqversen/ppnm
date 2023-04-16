use std::f64::consts::PI;
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
        return (int, error, 2);
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

pub fn monte_carlo_sampling(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, samples: u32, rng: &mut Rng) -> (f64, f64) {
    assert!(samples > 0);
    
    let dim = a.len();
    let volume = std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));
    let (mut sum, mut sum2) = (0.0, 0.0);
    let mut x = vec![0.0; dim];
    for _ in 0..samples {
        for i in 0..dim {x[i] = a[i] + rng.f64()*(b[i] - a[i])}
        let fx = f(&x);
        sum += fx;
        sum2 += fx * fx;
    }
    let mean = sum / samples as f64;
    let variance = sum2/samples as f64 - mean*mean;
    
    return (mean * volume, volume * f64::sqrt(variance / samples as f64))
}

fn corput_number(n: u32, base: u8) -> f64 {
    let mut n = n as f64;
    let base = base as f64;
    
    let mut q = 0.0;
    let mut bk = 1.0/base as f64;
    while n > 0.0 {
        q += (n % base) * bk;
        n /= base;
        bk /= base;
    }
    return q
}

fn halton_number(n: u32, dimension: usize) -> Vec<f64> {
    // prime numbers makes sure that all pairs are coprime
    const BASE: [u8; 19] = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67];
    assert!(dimension < BASE.len());
    return BASE[..dimension].iter().map(|b| corput_number(n, *b)).collect()
}

#[inline]
fn lattice_number(n:u32, alphas: &Vec<f64>) -> Vec<f64> {
    return alphas.iter().map(|a| (n as f64 * a).fract()).collect()
}

#[inline]
fn lattice_alphas(dimension: usize) -> Vec<f64> {
    (0..dimension).map(|i| (PI + i as f64).sqrt().fract()).collect()
}

pub fn low_descrepancy_sampling(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, samples: u32) -> (f64, f64) {
    assert!(samples > 0);
    
    let dimension = a.len();
    let volume = std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));
    let (mut sum1, mut sum2) = (0.0, 0.0);
    
    let alphas = lattice_alphas(dimension);
    let mut x = vec![0.0; dimension];
    for n in 0..samples/2 {
        for (i, random) in halton_number(n, dimension).iter().enumerate() {x[i] = a[i] + random*(b[i] - a[i])}
        let fx = f(&x);
        sum1 += fx;

        for (i, random) in lattice_number(n, &alphas).iter().enumerate() {x[i] = a[i] + random*(b[i] - a[i])}
        let fx = f(&x);
        sum2 += fx;
    }
    let int1 = volume * sum1 / (samples/2) as f64;
    let int2 = volume * sum2 / (samples/2) as f64;
    
    return (int1, (int1 - int2).abs())
}

pub fn recursive_stratified_sampling(f: &impl Fn(&Vec<f64>) -> f64, a: Vec<f64>, b: Vec<f64>, samples: u32, rng: &mut Rng) -> (f64, f64) {
    assert!(samples > 0);

    let dimension = a.len();
    let min_samples = (dimension * 32) as u32;
    let estimation_samples = samples / 5;

    if samples <= min_samples + estimation_samples {
        // cannot subdivide without subsiding min_sample
        return monte_carlo_sampling(f, a, b, samples, rng)
    }
    
    let (mut mean_left, mut mean_right) = (vec![0.0; dimension], vec![0.0; dimension]);
    let (mut variance_left, mut variance_right) = (vec![0.0; dimension], vec![0.0; dimension]);
    let (mut samples_left, mut samples_right) = (vec![0; dimension], vec![0; dimension]);
    
    // sample and store results
    let mut x = vec![0.0; dimension];
    for _ in 0..estimation_samples {
        for i in 0..dimension {x[i] = a[i] + rng.f64()*(b[i] - a[i])}
        let fx = f(&x);

        for i in 0..dimension {
            if x[i] > (a[i] + b[i]) / 2.0 {
                mean_right[i] += fx;
                variance_right[i] += fx*fx;
                samples_right[i] += 1;
            } else {
                mean_left[i] += fx;
                variance_left[i] += fx*fx;
                samples_left[i] += 1;
            }
        }
    }

    // make sums into actual means and variances
    for i in 0..dimension {
        if samples_left[i] == 0 || samples_right[i] == 0 {panic!("{} {} {} {}", i, samples_left[i], samples_right[i], estimation_samples)};
        mean_left[i] /= samples_left[i] as f64; // Test if samples_left/right is 0!!!
        mean_right[i] /= samples_right[i] as f64;
        
        variance_left[i] = variance_left[i] / samples_left[i] as f64 - mean_left[i]*mean_left[i];
        variance_right[i] = variance_right[i] / samples_right[i] as f64 - mean_right[i]*mean_right[i];
    }
    
    // find dimension which minimises total variance
    let mut index = 0;
    let mut min_variance = f64::INFINITY;
    for i in 0..dimension {
        let variance = variance_left[i]/(4*samples_left[i]) as f64 + variance_right[i]/(4*samples_right[i]) as f64;
        if variance < min_variance {
            min_variance = variance;
            index = i;
        }
    }

    // compute new left and right volume boundaries and sample sizes
    let mut a2 = a.clone();
    let mut b2 = b.clone();
    a2[index] = (a[index] + b[index]) / 2.0;
    b2[index] = (a[index] + b[index]) / 2.0;
    
    let remaining_samples = samples - estimation_samples;
    let mut sub_samples_left = (variance_left[index].sqrt() / (variance_left[index].sqrt() + variance_right[index].sqrt()) * remaining_samples as f64) as u32;
    let mut sub_samples_right = remaining_samples - sub_samples_left;

    // calculate integrals on each subdivision
    let volume = 0.5 * std::iter::zip(&a, &b).fold(1.0, |prod, (ai, bi)| prod * (bi - ai));

    let (int_left, error_left) = if sub_samples_left <= samples_left[index] {
        sub_samples_left = samples_left[index];
        (mean_left[index] * volume, volume * f64::sqrt(variance_left[index] / samples_left[index] as f64))
    } else {
        dbg!(sub_samples_left);
        recursive_stratified_sampling(f, a, b2, sub_samples_left, rng)
    };
    let (int_right, error_right) = if sub_samples_right <= samples_right[index] {
        sub_samples_right = samples_right[index];
        (mean_right[index] * volume, volume * f64::sqrt(variance_right[index] / samples_right[index] as f64))
    } else {
        dbg!(sub_samples_right);
        recursive_stratified_sampling(f, a2, b, sub_samples_right, rng)
    };

    // calculate total mean and variance
    let int = int_left + int_right;
    let error = f64::sqrt(error_left*error_left/(4*sub_samples_left) as f64 + error_right*error_right/(4*sub_samples_right) as f64);

    return (int, error)
}


#[cfg(test)]
mod tests {
    use super::*;
    
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
    fn test_monte_carlo_methods_1() {
        let samples = 1000 as u32;
        let expected = PI / 4.0;
        let mut rng = Rng::new(1234);
        let f = |x: &Vec<f64>| {
            if x.iter().fold(0.0, |sum, val| sum + val*val).sqrt() <= 1.0 {
                return 1.0
            } else {
                return 0.0
            }
        };

        let (int, err) = monte_carlo_sampling(&f, vec![0.0, 0.0], vec![1.0, 1.0], samples, &mut rng);
        assert!((int - expected).abs() < err);
        
        let (int, err) = low_descrepancy_sampling(&f, vec![0.0, 0.0], vec![1.0, 1.0], samples);
        assert!((int - expected).abs() < err);

        let (int, err) = recursive_stratified_sampling(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, &mut rng);
        dbg!(&expected, &int, &err);
        assert!((int - expected).abs() < err);
    }

    /* #[test]
    fn test_monte_carlo_methods_2() {
        let samples = 1e4 as u32;
        let expected = 13.0 / 20.0;
        let mut rng = Rng::new(1234);
        let f = |x: &Vec<f64>| {
            if x[0] >= 0.0 && x[1] <= 1.0 && x[1] >= x[0]*x[0] {
                return x[0] + x[1]
            } else {
                return 0.0
            }
        };
        
        let (int, err) = monte_carlo_sampling(&f, vec![0.0, 0.0], vec![1.0, 1.0], samples, &mut rng);
        assert!((int - expected).abs() < err);
        
        let (int, err) = low_descrepancy_sampling(&f, vec![0.0, 0.0], vec![1.0, 1.0], samples);
        assert!((int - expected).abs() < err);

        let (int, err) = recursive_stratified_sampling(&f, vec![0.0, 0.0], vec![2.0, 2.0], samples, &mut rng);
        dbg!(&expected, &int, &err);
        assert!((int - expected).abs() < err);
    } */
}