use std::iter::zip;

#[derive(Debug)]
pub struct RungeKuttaStepper {
    a: Vec<Vec<f64>>,
    b: Vec<f64>,
    b_star: Vec<f64>,
    c: Vec<f64>,
    s: usize,
}

fn in_place_add(a: &mut Vec<f64>, b: Vec<f64>) {
    zip(a.iter_mut(), b).for_each(|(ai, bi)| *ai += bi);
}

impl RungeKuttaStepper {
    pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>, b_star: Vec<f64>, c: Vec<f64>) -> Self {
        // Method consistency
        assert!(c[0] == 0.0);
        assert!(b.iter().sum::<f64>() - 1.0 < f64::EPSILON);
        assert!(b_star.iter().sum::<f64>() - 1.0 < f64::EPSILON);
        
        // Size compatibility
        let s = b.len();
        assert!(b_star.len() == s);
        assert!(c.len() == s);
        assert!(a.len() == s);
        for (i, row) in a.iter().enumerate() {assert!(row.len() >= i)}
        
        Self {a: a, b: b, b_star: b_star, c: c, s: s}
    }

    fn compute_ks(&self, f: fn(f64, Vec<f64>) -> Vec<f64>, x: f64, y: &Vec<f64>, h: f64) -> Vec<Vec<f64>> {
        let mut ks: Vec<Vec<f64>> = Vec::with_capacity(self.s);
        for n in 0..self.s {
            let xn = x + self.c[n] * h;
            let mut yn = y.clone();
            for i in 0..n {
                let term: Vec<f64> = ks[i].iter().map(|val| self.a[n][i] * h * val).collect();
                in_place_add(&mut yn, term)
            }
            let kn = (f)(xn, yn);
            ks.push(kn);
        }

        return ks;
    }
    
    pub fn step(&self, f: fn(f64, Vec<f64>) -> Vec<f64>, x: f64, mut y: Vec<f64>, h: f64) -> (Vec<f64>, Vec<f64>) {
        let mut error = vec![0.0; y.len()];
        let ks = self.compute_ks(f, x, &y, h);
        for (i, ki) in ks.iter().enumerate() {
            let term = ki.iter().map(|val| self.b[i] * h * val).collect();
            in_place_add(&mut y, term);
            
            let term = ki.iter().map(|val| (self.b[i] - self.b_star[i]) * h * val).collect();
            in_place_add(&mut error, term);
        }
        return (y, error);
    }
}

impl RungeKuttaStepper {
    pub fn rk12() -> Self {
        let a = vec![
            vec![],
            vec![1.0]
        ];
        let b = vec![0.5, 0.5];
        let b_star = vec![1.0, 0.0];
        let c = vec![0.0, 1.0];
        return Self::new(a, b, b_star, c);
    }

    pub fn rk23() -> Self {
        let a = vec![
            vec![],
            vec![1.0/2.0],
            vec![0.0,       3.0/4.0],
            vec![2.0/9.0,   1.0/3.0,    4.0/9.0]
        ];
        let b = vec![2.0/9.0, 1.0/3.0, 4.0/9.0, 0.0];
        let b_star = vec![7.0/24.0, 1.0/4.0, 1.0/3.0, 1.0/8.0];
        let c = vec![0.0, 1.0/2.0, 3.0/4.0, 1.0];
        return Self::new(a, b, b_star, c);
    }

    pub fn rk34() -> Self {
        let a = vec![
            vec![], 
            vec![1.0/5.0],
            vec![3.0/40.0,      9.0/40.0],
            vec![3.0/10.0,      -9.0/10.0,  6.0/5.0],
            vec![-11.0/54.0,    5.0/2.0,    -70.0/27.0,   -35.0/27.0],
            vec![1631.0/55296.0,175.0/512.0,575.0/13824.0,44275.0/110592.0, 253.0/4096.0],
        ];
        let b = vec![
            37.0/378.0, 0.0, 250.0/621.0, 125.0/594.0, 0.0, 512.0/1771.0
        ];
        let b_star = vec![
            2825.0/27648.0, 0.0, 18575.0/48384.0, 13525.0/55296.0, 277.0/14336.0, 1.0/4.0
        ];
        let c = vec![
            0.0, 1.0/5.0, 3.0/10.0, 3.0/5.0, 1.0, 7.0/8.0
        ];
        return Self::new(a, b, b_star, c);
    }
    
    pub fn rk45() -> Self {
        let a = vec![
            vec![], 
            vec![1.0/4.0],
            vec![3.0/32.0,      9.0/32.0],
            vec![1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0],
            vec![439.0/216.0,   -8.0,           3680.0/513.0,   -845.0/4104.0],
            vec![-8.0/27.0,     2.0,            -3544.0/2565.0, 1859.0/4104.0,   -11.0/40.0],
        ];
        let b = vec![
            16.0/135.0,     0.0,            6656.0/12825.0, 28561.0/56430.0, -9.0/50.0, 2.0/55.0
        ];
        let b_star = vec![
            25.0/216.0,     0.0,            1408.0/2565.0,  2197.0/4104.0,  -1.0/5.0,   0.0
        ];
        let c = vec![
            0.0,            1.0/4.0,        3.0/8.0,        12.0/13.0,      1.0,        1.0/2.0
        ];
        return Self::new(a, b, b_star, c);
    }
}

#[derive(Debug)]
pub struct AdaptiveStepSizeDriver {
    stepper: RungeKuttaStepper,
    abs: f64, 
    rel: f64,
}

impl AdaptiveStepSizeDriver {
    pub fn new(stepper: RungeKuttaStepper) -> Self {
        Self{stepper: stepper, abs: 0.01, rel: 0.01}
    }

    pub fn run(&self, f: fn(f64, Vec<f64>) -> Vec<f64>, a: f64, ya: Vec<f64>, b: f64, mut h: f64) -> (Vec<f64>, Vec<Vec<f64>>) {
        assert!(a < b);

        let mut x = a;
        let mut y = ya;
        let mut x_list = vec![x];
        let mut y_list = vec![y.clone()];
        loop {
            if x >= b {return (x_list, y_list)}  // job done
            if x+h > b {h = b-x;}  // last step should end at b
            let (yh, mut err) = self.stepper.step(f, x, y.clone(), h);
            let tol: Vec<f64> = yh.iter().map(|ti| f64::max(self.abs, ti.abs() * self.rel) * h/(b-a)).collect();
            err.iter_mut().for_each(|ei| *ei=ei.abs());
            let mut ok = true;
            for (ti, ei) in zip(&tol, &err) {
                if !(ei<ti) {ok = false; break;}
            }
            if ok {
                x += h;
                y = yh;
                x_list.push(x);
                y_list.push(y.clone());
            }
            let factor = zip(&tol, &err).map(|(ti, ei)| ti/ei).fold(f64::INFINITY, |a, b| a.min(b));
            h *= f64::min(factor.powf(0.25) * 0.95, 2.0);  // adapt stepsize
        }
    }
}