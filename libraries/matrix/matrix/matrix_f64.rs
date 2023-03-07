use super::Matrix;

impl Matrix<f64> {
    pub fn zeros(num_rows: usize, num_cols: usize) -> Self {
        let data: Vec<f64> = vec![0.0; num_rows * num_cols];
        Self::new_from_data(data, num_rows, num_cols)
    }

    pub fn ones(num_rows: usize, num_cols: usize) -> Self {
        let data: Vec<f64> = vec![1.0; num_rows * num_cols];
        Self::new_from_data(data, num_rows, num_cols)
    }

    pub fn idty(dim: usize) -> Self {
        let mut result = Self::zeros(dim, dim);
        for i in 0..dim {result.data[(dim + 1) * i] = 1.0};
        return result;
    }

    pub fn zeros_like(matrix: &Self) -> Self {
        Self::zeros(matrix.num_rows, matrix.num_cols)
    }

    pub fn norm(&self) -> f64 {
        f64::sqrt(self.iter().map(|x| x.0*x.0).sum())
    }
}