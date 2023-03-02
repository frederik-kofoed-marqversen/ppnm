use std::ops;

#[derive(Debug, Clone)]
pub struct Matrix {
    pub data: Vec<Vec<f64>>,
    pub num_cols: usize,
    pub num_rows: usize,
}

impl Matrix {
    pub fn new(data: Vec<Vec<f64>>) -> Self {
        Self {num_cols: data.len(), num_rows: data[0].len(), data: data,}
    }

    pub fn zeros(num_rows: usize, num_cols: usize) -> Self {
        Self::new(vec![vec![0.0; num_rows]; num_cols])
    }

    pub fn ones(num_rows: usize, num_cols: usize) -> Self {
        Self::new(vec![vec![1.0; num_rows]; num_cols])
    }

    pub fn idty(dim: usize) -> Self {
        let mut result = Self::zeros(dim, dim);
        for i in 0..dim {result.data[i][i] = 1.0};
        return result;
    }

    pub fn zeros_like(matrix: &Self) -> Self {
        Self::zeros(matrix.num_rows, matrix.num_cols)
    }

    pub fn norm_sq(&self) -> f64{
        self.data.iter().map(|col| col.iter().map(|x| x*x).sum::<f64>()).sum()
    }

    pub fn mul(a: &Self, b: &Self) -> Self {
        let mut c = Self::zeros(a.num_rows, b.num_cols);
        for (k, ak) in a.data.iter().enumerate() {
            for (j, bj) in b.data.iter().enumerate() {
                let bjk = bj[k];
                let cj = &mut c.data[j];
                for (i, aki) in ak.iter().enumerate() {
                    cj[i] += aki * bjk;
                }
            }
        }
        return c
    }

    pub fn transpose(&self) -> Self {
        let mut result: Self = Self::zeros(self.num_cols, self.num_rows);
        for j in 0..self.num_rows {
            for i in 0..self.num_cols {
                result.data[j][i] = self.data[i][j];
            }
        }
        return result;
    }
}

impl ops::Mul<&Matrix> for &Matrix {
    type Output = Matrix;
    fn mul(self, rhs: &Matrix) -> Matrix {Matrix::mul(self, rhs)}
}