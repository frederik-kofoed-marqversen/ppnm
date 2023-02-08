fn main() {
    let sqrt2 = (2 as f32).sqrt();
    println!("sqrt(2) = {sqrt2}");
    println!("sqrt2 * sqrt2 = {} (should be equal to 2)", sqrt2 * sqrt2);

    println!("gamma(1) = {}", sfuns::gamma(1.0));
    println!("gamma(2) = {}", sfuns::gamma(2.0));
    println!("gamma(3) = {}", sfuns::gamma(3.0));
    println!("gamma(10) = {}", sfuns::gamma(10.0));
} 
