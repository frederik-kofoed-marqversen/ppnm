fn main(){
/*
    // T ASK 1
    println!("Task 1");
    let mut i : u32 = 1;
    while i + 1 > i {i += 1;}
    println!("My max un-signed int = {i}");
    println!("std::u32::MAX = {}", u32::MAX);
    
    let mut i : i32 = 1;
    while i + 1 > i {i += 1;}
    println!("My max signed int = {i}");
    println!("std::i32::MAX = {}", i32::MAX);

    let mut i : i32 = 1;
    while i - 1 < i {i -= 1;}
    println!("My min int = {i}");
    println!("std::i32::MIN = {}", i32::MIN);
*/
    // TASK 2
    println!("\n Task 2");
    let mut x : f32 = 1.0;
    while 1.0 + x != 1.0 {x /= 2.0;};
    x *= 2.0;
    println!("My f32 epsilon = {}", x);
    println!("std::f32::EPSILON = {}", f32::EPSILON);

    let mut y : f64 = 1.0;
    while 1.0 + y != 1.0 {y /= 2.0;};
    y *= 2.0;
    println!("My f64 epsilon = {}", y);
    println!("std::f64::EPSILON = {}", f64::EPSILON);
 
    // TASK 3
    println!("\n Task 3");
    let n = 1e6 as u32;
    let epsilon = f64::EPSILON;
    let tiny = epsilon / 2.0;
    let (mut sum_a, mut sum_b) = (0.0, 0.0);

    sum_a += 1.0;
    for _i in 0..n {
        sum_a += tiny;
        sum_b += tiny;
    }
    sum_b += 1.0;
    println!("sum_a-1 = {} should be {}", sum_a-1.0, (n as f64)*tiny);
    println!("sum_b-1 = {} should be {}", sum_b-1.0, (n as f64)*tiny);

    // TASK 4
    println!("\n Task 4");
    let d1 = 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1 + 0.1;
    let d2 = 8.0 * 0.1;
    println!("d1={d1}");
    println!("d2={d2}");
    println!("d1==d2 ? => {}", d1 == d2);

    fn are_close(a: f64, b: f64) -> bool {
        let acc = 1e-9;
        let eps = 1e-9;
        if f64::abs(a - b) < acc {return true}
        if f64::abs(a - b) < (a.abs() + b.abs()) * eps {return true}
        return false;
    }
    println!("d1~=d2 ? => {}", are_close(d1, d2));
}