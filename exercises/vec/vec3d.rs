use std::{ops, fmt};

fn are_close(a: f64, b: f64) -> bool {
    let acc = 1e-9;
    let eps = 1e-9;
    if f64::abs(a - b) < acc {return true}
    if f64::abs(a - b) < (a.abs() + b.abs()) * eps {return true}
    return false;
}

pub struct Vec3d {
    x : f64,
    y : f64,
    z : f64,
}

impl Vec3d {
    pub fn new(x0: f64, y0: f64, z0: f64) -> Self {
        Self {x: x0, y:y0, z:z0}
    }
    
    pub fn norm(&self) -> f64 {
        (self.x.powf(2.0) + self.y.powf(2.0) + self.z.powf(2.0)).sqrt()
    }

    pub fn dot(&self, rhs: &Vec3d) -> f64 {
        self.x * rhs.x + self.y * rhs.y + self.z * rhs.z
    }

    pub fn cross(&self, rhs: &Vec3d) -> Vec3d{
        Self {
            x: self.y * rhs.z - self.z * rhs.y,
            y: self.z * rhs.x - self.x * rhs.z,
            z: self.x * rhs.y - self.y * rhs.x,
        }
    }

    pub fn are_close(u: &Vec3d, v: &Vec3d) -> bool {
        if !are_close(u.x, v.x) {return false}
        if !are_close(u.y, v.y) {return false}
        if !are_close(u.z, v.z) {return false}
        return true;
    }
}

impl ops::Add<&Vec3d> for &Vec3d {
    type Output = Vec3d;
    fn add(self, rhs: &Vec3d) -> Vec3d {
        Vec3d {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

impl ops::Neg for &Vec3d {
    type Output = Vec3d;
    fn neg(self) -> Vec3d {
        Vec3d {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl ops::Sub<&Vec3d> for &Vec3d {
    type Output = Vec3d;
    fn sub(self, rhs: &Vec3d) -> Vec3d {
        self + &(-rhs)
    }
}

impl ops::Mul<f64> for &Vec3d {
    type Output = Vec3d;
    fn mul(self, c: f64) -> Vec3d {
        Vec3d {
            x: self.x * c,
            y: self.y * c,
            z: self.z * c,
        }
    }
}

impl ops::Mul<&Vec3d> for f64 {
    type Output = Vec3d;
    fn mul(self, u: &Vec3d) -> Vec3d {
        u * self
    }
}

impl fmt::Display for Vec3d {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[{} {} {}]", self.x, self.y, self.z)
    }
}