extern crate matrix;

use matrix::Matrix;
use matrix::linalg::lstsq;
use std::iter::zip;
use std::io::{Write, BufRead};

fn parse_line(line: &str) -> std::io::Result<Vec<f64>> {
    let mut items = line.split([':', ',']).filter(|x| !x.is_empty()).map(|x| x.trim());
    items.next(); // first item is the desciptor
    return Ok(items.map(
        |x| x.parse().expect(&format!("Could not interpret '{x}' as an f64"))
    ).collect())
}

fn main() -> std::io::Result<()> {
    let in_file = std::fs::File::open("data.txt")?;
    let mut lines = std::io::BufReader::new(in_file).lines();

    let ts = Matrix::new(vec![parse_line(&lines.next().unwrap()?)?]);
    let mut ln_ys = Matrix::new(vec![parse_line(&lines.next().unwrap()?)?]);
    let mut ln_dys = Matrix::new(vec![parse_line(&lines.next().unwrap()?)?]);

    let mut out_file = std::fs::File::create("out.txt")?;
    for ((t, y), dy) in zip(zip(ts.iter(), ln_ys.iter()), ln_dys.iter()) {
        out_file.write(format!("{t}\t{y}\t{dy}\n").as_bytes())?;
    }
    out_file.write(b"\n\n")?;

    for (y, dy) in zip(ln_ys.iter_mut(), ln_dys.iter_mut()) { // take ln of ys and dys
        *dy = *dy / *y;
        *y = f64::ln(*y);
    }

    let ln_fns: Vec<&dyn Fn(f64) -> f64> = vec![&|_| 1.0, &|x| -x];

    let (c, mut sigma) = lstsq::fit(&ln_fns, &ts, &ln_ys, &ln_dys);
    lstsq::correlations(&mut sigma);
    
    let half_life_true = 3.6319;
    let half_life = f64::ln(2.0) / c[0][1];
    let half_life_max = f64::ln(2.0) / (c[0][1] - sigma[1][1]);
    let half_life_min = f64::ln(2.0) / (c[0][1] + sigma[1][1]);
    println!("\nHalf life (lambda/day): {half_life} in [{half_life_min}, {half_life_max}]");
    println!(
        "The correct val is {half_life_true} which is {} within the bounds.\n", 
        if half_life_min < half_life_true && half_life_true < half_life_max {"INDEED"} else {"UNFORTUNATELY NOT"}
    );
    
    let fit = |t| f64::exp(c[0][0] - c[0][1] * t);
    let fit_max = |t| f64::exp(c[0][0] + sigma[0][0] - (c[0][1] - sigma[1][1]) * t);
    let fit_min = |t| f64::exp(c[0][0] - sigma[0][0] - (c[0][1] + sigma[1][1]) * t);

    let num_points = 100.0;
    for i in 0..=num_points as u32 {
        let t = i as f64 / num_points * 16.0;
        let y = fit(t);
        out_file.write(format!("{t}\t{y}\n").as_bytes())?;
    }
    out_file.write(b"\n\n")?;
    
    for i in 0..=num_points as u32 {
        let t = i as f64 / num_points * 16.0;
        let y = fit_max(t);
        out_file.write(format!("{t}\t{y}\n").as_bytes())?;
    }
    out_file.write(b"\n")?;
    for i in 0..=num_points as u32 {
        let t = i as f64 / num_points * 16.0;
        let y = fit_min(t);
        out_file.write(format!("{t}\t{y}\n").as_bytes())?;
    }
    out_file.write(b"\n\n")?;

    return Ok(())
}