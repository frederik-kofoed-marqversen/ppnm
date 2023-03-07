use super::Matrix;
use std::ops::Mul;

macro_rules! matrix_multiply_scalar {
    ($MatrixType:ty, $T:tt) => {
        impl<$T: Mul<Output = $T> + Copy> Mul<$T> for $MatrixType {
            type Output = Matrix<$T>;
            fn mul(self, scalar: $T) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for element in self.data.iter() {
                    data.push(scalar * *element);
                }
                return Matrix::new_from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
  }
  matrix_multiply_scalar!(Matrix<T>, T);
  matrix_multiply_scalar!(&Matrix<T>, T);

  /* macro_rules! scalar_multiply_matrix {
    ($LHS:ty, $T:tt) => {
        impl<$T: Mul<Output = $T> + Copy> Mul<$LHS> for $T {
            type Output = Matrix<$T>;
            fn mul(self, matrix: $T) -> Matrix<$T> {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for element in matrix.data.iter() {
                    data.push(self * *element);
                }
                return Matrix::new_from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
  }
  scalar_multiply_matrix!(Matrix<T>, T);
  scalar_multiply_matrix!(&Matrix<T>, T); */