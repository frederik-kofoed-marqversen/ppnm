extern crate vec3d;
use vec3d::Vec3d;
fn main(){
    let u = Vec3d::new(1.0, 1.0, 1.0);
    let v = 2.0 * &u;
    println!("Length of vector = {}", u.norm());
    println!("vector u = {u}");
    println!("vector v = {v}");
    println!("v - u = {}", &v - &u);
    println!("v * u = {}", v.dot(&u));
    println!("v x u = {}", v.cross(&u));

    let w = Vec3d::new(1.0, -1.0, 1.0);
    println!("vector w = {w}");
    println!("v x w = {}", v.cross(&w));
    println!("(v x w) * v = {}", v.cross(&w).dot(&v));

    let x = Vec3d::new(0.1, 0.1, 0.1);
    println!("vector x = {}", x);
    let mut y = Vec3d::new(0.0, 0.0, 0.0);
    for _i in 0..8 {y = &y + &x};
    println!("x + x + ... + x = {}", y);

    println!("8 * x = {}", &x * 8.0);
    println!("They are {}close!", if Vec3d::are_close(&(&x * 8.0), &y) {""} else {"not "});

    println!("7 * x = {}", &x * 7.0);
    println!("They are {}close!", if Vec3d::are_close(&(&x * 7.0), &y) {""} else {"not "});
}