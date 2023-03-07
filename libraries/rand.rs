// Hashing functions taken from the fastrand crate

use std::cell::Cell;

#[derive(Debug)]
pub struct Rng(Cell<u64>);

impl Rng {
    pub fn new(seed: u64) -> Self {
        let rng = Rng(Cell::new(0));
        rng.set_state(seed);
        rng
    }

    pub fn set_state(&self, seed: u64) {
        self.0.set(seed);
    }

    fn step(&self) -> u64 {
        let s = self.0.get().wrapping_add(0xA0761D6478BD642F);
        self.0.set(s);
        s
    }

    pub fn u64(&self) -> u64 {
        let s = self.step();
        let t = u128::from(s) * u128::from(s ^ 0xE7037ED1A0B428DB);
        (t as u64) ^ (t >> 64) as u64
    }

    pub fn f64(&self) -> f64 {
        let b = 64;
        let f = std::f64::MANTISSA_DIGITS - 1;
        f64::from_bits((1 << (b - 2)) - (1 << f) + (self.u64() >> (b - f))) - 1.0
    }
}