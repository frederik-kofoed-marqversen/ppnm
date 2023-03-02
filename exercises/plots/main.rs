extern crate sfuns;
use std::io::Write;

fn main() -> std::io::Result<()> {
    let mut file = std::fs::File::create("target/factorials.data")?;
    let mut factorial = 1;
    file.write_all(format!("0  {factorial}\n").as_bytes())?;
    for n in 1..5 {
        factorial *= n;
        file.write_all(format!("{n}  {factorial}\n").as_bytes())?;
    }

    let mut file = std::fs::File::create("target/gamma.data")?;
    let num_point = 500.0;
    let range = (-5.0, 5.0);
    for n in 0..=num_point as i32 {
        let x = f64::from(n) / num_point * (range.1 - range.0) + range.0;
        let gamma = sfuns::gamma(x);
        file.write_all(format!("{x}  {gamma}\n").as_bytes())?;
    }

    let mut file = std::fs::File::create("target/lngamma.data")?;
    let num_point = 500.0;
    let range = (0.0, 20.0);
    for n in 0..num_point as i32 {
        let x = f64::from(n) / num_point * (range.1 - range.0) + range.0;
        let lngamma = sfuns::lngamma(x);
        file.write_all(format!("{x}  {lngamma}\n").as_bytes())?;
    }

    Ok(())
}