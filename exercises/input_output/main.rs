use std::io::{BufRead, Write, Read};

fn main() -> std::io::Result<()> {
    let mut input_file = String::new();
    let mut output_file = String::new();

    let args: Vec<String> = std::env::args().collect();

    for arg in &args[1..] {  // ignore first argument since this will always be ./main.bin
        let words: Vec<&str> = arg.split(":").collect();
        let command = words[0];
        let argument = words[1].to_string();
        match command {
            "-input" | "-i" => input_file = argument,
            "-output" | "-o" => output_file = argument,
            "-numbers" => println!("Parsing '-numbers' argument \n{}", parse_string(&argument, vec![','])),
            _ => println!("Command '{command}' not found."),
        }
    }

    println!("Parsing input from stdin \nPress Enter to Exit");
    let mut lines = std::io::stdin().lock().lines();
    while let Some(line) = lines.next() {
        let input = line.unwrap();
        if input.len() == 0 {break};
        println!("{}", parse_string(&input, vec![' ', '\t', '\n']));
        println!("Press Enter to Exit");
    }

    if input_file.len() > 0 {
        let mut file = std::fs::File::open(&input_file)?;
        let mut contents = String::new();
        file.read_to_string(&mut contents)?;
        let output = parse_string(&contents, vec![' ', '\t', '\n']);
        
        if output_file.len() > 0 {
            println!("Parsing file '{input_file}' and saving to '{output_file}'");
            let mut file = std::fs::File::create(output_file)?;
            file.write_all(output.as_bytes())?;
        } else {
            println!("Parsing file '{input_file}'");
            println!("{}", output);
        }
    }

    return Ok(())
}

fn parse_string(string_to_parse: &str, split_delimiters: Vec<char>) -> String {
    let mut result = String::new();
    for entry in string_to_parse.split(&split_delimiters[..]).filter(|&x| !x.is_empty()) {
        let number = entry.parse::<f64>().expect(&format!("Could not interpret '{entry}' as an f64"));
        result += &format!("{number} {} {}\n", f64::sin(number), f64::cos(number));
    }
    return result;
}