use std::f64::consts::PI;

pub fn factorial(num: u64) -> u64 {(1..=num).product()}

pub fn gamma(x: f64) -> f64{
    // single precision error function (Abramowitz and Stegun, from Wikipedia)
    if x < 0.0 {return PI / f64::sin(PI * x) / gamma(1.0 - x)};
    if x < 9.0 {return gamma(x + 1.0) / x};
    let lngamma = x * f64::ln(x + 1.0 / (12.0 * x - 1.0 / x / 10.0)) - x + f64::ln(2.0 * PI / x) / 2.0;
    return f64::exp(lngamma);
}

pub fn lngamma(x: f64) -> f64 {
    // single precision gamma function (formula fom wikipedia)
    if x < 0.0 {panic!("Cannnot compute lngamma for negative numbers");}
    if x < 9.0 {return lngamma(x + 1.0) - f64::ln(x);}
    let lngamma = x * f64::ln(x + 1.0 / (12.0 * x - 1.0 / x / 10.0)) - x + f64::ln(2.0 * PI / x) / 2.0;
    return lngamma
}

pub fn erf(x: f64) -> f64{
    // single precision error function (Abramowitz and Stegun, from Wikipedia)
    if x < 0.0 {return -erf(-x)};
    let a = [0.254829592, -0.284496736, 1.421413741, -1.453152027, 1.061405429];
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let sum = t * (a[0] + t * (a[1] + t * (a[2] + t * (a[3] + t * a[4]))));/* the right thing */
    return 1.0 - sum * f64::exp(- x * x);
}

pub fn are_close(a: f64, b: f64) -> bool {
    let acc = 1e-9;
    let eps = 1e-9;
    if f64::abs(a - b) < acc {return true}
    if f64::abs(a - b) < f64::max(a.abs(), b.abs()) * eps {return true}
    return false;
}