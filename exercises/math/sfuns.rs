use std::f64::consts::PI;

pub fn gamma(x: f64) -> f64 {
    // single precision gamma function (formula fom wikipedia)
    if x < 0.0 {
        return PI / f64::sin(PI * x) / gamma(1.0 - x);  // Euler's reflection formula
    }
    if x < 9.0 {
        return gamma(x + 1.0) / x;  // Recurrence relation
    }
    
    let lngamma = x * f64::ln(x + 1.0 / (12.0 * x - 1.0 / x / 10.0)) - x + f64::ln(2.0 * PI / x) / 2.0;
    return f64::exp(lngamma)
}