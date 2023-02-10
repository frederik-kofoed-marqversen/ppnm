use std::f64::consts::PI;

pub fn lngamma(x: f64) -> f64 {
    // single precision gamma function (formula fom wikipedia)
    if x < 0.0 {
        return f64::ln(PI / f64::sin(PI * x)) - lngamma(1.0 - x);  // Euler's reflection formula
    }
    if x < 9.0 {
        return lngamma(x + 1.0) - f64::ln(x);  // Recurrence relation
    }
    
    let lngamma = x * f64::ln(x + 1.0 / (12.0 * x - 1.0 / x / 10.0)) - x + f64::ln(2.0 * PI / x) / 2.0;
    return lngamma
}