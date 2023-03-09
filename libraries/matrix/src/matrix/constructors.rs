use super::Matrix;

impl<T: From<u8> + Copy> Matrix<T> {
    pub fn zeros(num_rows: usize, num_cols: usize) -> Self {
        let data: Vec<T> = vec![T::from(0); num_rows * num_cols];
        Self::from_data(data, num_rows, num_cols)
    }

    pub fn ones(num_rows: usize, num_cols: usize) -> Self {
        let data: Vec<T> = vec![T::from(1); num_rows * num_cols];
        Self::from_data(data, num_rows, num_cols)
    }

    pub fn idty(dim: usize) -> Self {
        let mut result = Self::zeros(dim, dim);
            for i in 0..dim {result.data[(dim + 1) * i] = T::from(1)};
            return result;
    }

    pub fn zeros_like(matrix: &Self) -> Self {
        Self::zeros(matrix.num_rows, matrix.num_cols)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_idty_1() {
        let mat1: Matrix<f64> = Matrix::idty(3);
        let data: Vec<f64> = vec![1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let mat2 = Matrix::from_data(data, 3, 3);
        assert_eq!(
            mat1,
            mat2
        );
    }

    #[test]
    fn test_idty_2() {
        let mat1: Matrix<u8> = Matrix::idty(3);
        let data: Vec<u8> = vec![1, 0, 0, 0, 1, 0, 0, 0, 1];
        let mat2 = Matrix::from_data(data, 3, 3);
        assert_eq!(
            mat1,
            mat2
        );
    }

    #[test]
    fn test_idty_3() {
        let mat1: Matrix<usize> = Matrix::idty(3);
        let data: Vec<usize> = vec![1, 0, 0, 0, 1, 0, 0, 0, 1];
        let mat2 = Matrix::from_data(data, 3, 3);
        assert_eq!(
            mat1,
            mat2
        );
    }

    #[test]
    fn test_ones_1() {
        let mat1: Matrix<f32> = Matrix::ones(2, 3);
        let data: Vec<f32> = vec![1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
        let mat2 = Matrix::from_data(data, 2, 3);
        assert_eq!(
            mat1,
            mat2
        );
    }
}