fn main() {
    let sqrt2 = (2 as f32).sqrt();
    println!("sqrt(2) = {sqrt2}");
    println!("sqrt2 * sqrt2 = {} (should be equal to 2)", sqrt2 * sqrt2);

    println!("gamma(1) = {}", f64::exp(sfuns::lngamma(1.0)));
    println!("gamma(2) = {}", f64::exp(sfuns::lngamma(2.0)));
    println!("gamma(3) = {}", f64::exp(sfuns::lngamma(3.0)));
    println!("gamma(10) = {}", f64::exp(sfuns::lngamma(10.0)));

    println!("ln(gamma(10000)) = {}", sfuns::lngamma(10000.0));
} 
