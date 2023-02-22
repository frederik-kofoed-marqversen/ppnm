use std::thread;

#[derive(Clone,Copy,Debug)]
struct Data { 
    a : u32,
    b : u32,
    sum : f64,
}

impl Data {
    fn new() -> Self {Self {a: 0, b:0, sum:0.0}}
}

fn harmonic(obj: &mut Data) {
	obj.sum = 0.0;
    for i in obj.a..obj.b {
        obj.sum += 1.0/i as f64;
    }
}

fn main() {
    let mut num_terms : u32 = 1e8 as u32;
    let mut num_threads : u32 = 1;

    let args: Vec<String> = std::env::args().collect();
    for arg in &args[1..] {  // ignore first argument since this will always be ./main.bin
        let words: Vec<&str> = arg.split(":").collect();
        match words[0] {
            "-terms" => num_terms = words[1].parse::<u32>().unwrap(),
            "-threads" => num_threads = words[1].parse::<u32>().unwrap(),
            _ => println!("Command '{}' not found.", words[1]),
        }
    }

    println!("num_terms={num_terms} num_threads={num_threads}");

    let mut intervals: Vec<Data> = Vec::with_capacity(num_threads as usize);
    for i in 0..num_threads {
        let mut obj = Data::new();
        obj.a = 1 + (i * num_terms) / num_threads;
        obj.b = 1 + ((i + 1) * num_terms) / num_threads;
        intervals.push(obj);
    }

    let mut children = vec![];
    for mut interval in intervals.iter().cloned() {
        children.push(thread::spawn(move || {
            harmonic(&mut interval);
            interval
        }));
    }

    for (i, child) in children.into_iter().enumerate() {
        intervals[i] = child.join().unwrap();
    }
    
    let mut sum: f64 = 0.0;
    for interval in intervals {
        // println!("{:?}", interval);
        sum += interval.sum;
    }

    println!("Total sum = {sum}");

}