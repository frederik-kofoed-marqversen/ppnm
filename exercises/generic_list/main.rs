use std::ptr::NonNull;
use std::alloc::{self, Layout};

#[derive(Debug)]
struct MyVec<T> {
    ptr: NonNull<T>,
    cap: usize,
    len: usize,
}

impl<T> MyVec<T> {
    pub fn new() -> Self {
        assert!(std::mem::size_of::<T>() != 0, "do not support zero-sized types"); // this is apparently important
        MyVec {
            ptr: NonNull::dangling(),
            len: 0,
            cap: 0,
        }
    }

    fn grow(&mut self) {
        let (new_cap, new_layout) = if self.cap == 0 {
            (1, Layout::array::<T>(1).unwrap())
        } else {
            let new_cap = 2 * self.cap;
            let new_layout = Layout::array::<T>(new_cap).unwrap();
            (new_cap, new_layout)
        };

        let new_ptr = if self.cap == 0 {
            unsafe { alloc::alloc(new_layout) }
        } else {
            let old_layout = Layout::array::<T>(self.cap).unwrap();
            let old_ptr = self.ptr.as_ptr() as *mut u8;
            unsafe { alloc::realloc(old_ptr, old_layout, new_layout.size()) }
        };

        // If allocation fails we must abort using handle_alloc_error()
        // new_ptr will in that case be Null-type which is caught by NonNull
        self.ptr = match NonNull::new(new_ptr as *mut T) {
            Some(p) => p,
            None => alloc::handle_alloc_error(new_layout),
        };
        self.cap = new_cap;
    }

    pub fn push(&mut self, elem: T) {
        if self.len == self.cap { self.grow(); }
        unsafe {
            std::ptr::write(self.ptr.as_ptr().add(self.len), elem);
        }
        self.len += 1;
    }

    pub fn pop(&mut self) -> Option<T> {
        if self.len == 0 {
            None
        } else {
            self.len -= 1;
            unsafe {
                Some(std::ptr::read(self.ptr.as_ptr().add(self.len)))
            }
        }
    }

    pub fn read(&self, idx: usize) -> T {
        assert!(idx < self.len); // make sure idx is good
        unsafe {
            std::ptr::read(self.ptr.as_ptr().add(idx))
        }
    }
}

// In proper Rust Drop is implemented such that when Vec comes out of scope its memory is freed
impl<T> Drop for MyVec<T> {
    fn drop(&mut self) {
        if self.cap != 0 {
            while let Some(_) = self.pop() { }
            let layout = Layout::array::<T>(self.cap).unwrap();
            unsafe {
                alloc::dealloc(self.ptr.as_ptr() as *mut u8, layout);
            }
        }
    }
}

use std::io::Read;
fn main() -> std::io::Result<()> {
    let mut input_file = String::new();

    let args: Vec<String> = std::env::args().collect();
    for arg in &args[1..] {  // ignore first argument since this will always be ./main.bin
        let words: Vec<&str> = arg.split(":").collect();
        match words[0] {
            "-input" | "-i" => input_file = words[1].to_string(),
            _ => panic!("Command '{}' not found.", words[0]),
        }
    }

    assert!(input_file.len() > 0);

    let mut file = std::fs::File::open(&input_file)?;
    let mut contents = String::new();
    file.read_to_string(&mut contents)?;
    
    let mut numbers = MyVec::new();
    for entry in contents.split([' ', '\t', '\n']).filter(|&x| !x.is_empty()) {
        numbers.push(entry.parse::<f64>().unwrap());
    }

    for i in 0..numbers.len {
        println!("number: {:+1.4e}", numbers.read(i));
    }

    return Ok(())
}