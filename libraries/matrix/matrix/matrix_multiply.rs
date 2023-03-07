use super::Matrix;
use std::ops::{Add, Mul};
use std::iter::Sum;

// MULTIPLY TWO MATRICES
macro_rules! matrix_multiply {
    ($LHS:ty, $RHS:ty, $T:tt) => {
        impl<$T: Add<Output = $T> + Mul<Output = $T> + Copy + Sum> Mul<$RHS> for $LHS {
            type Output = Matrix<$T>;
            fn mul(self, other: $RHS) -> Self::Output {
                if self.num_cols != other.num_rows {
                    panic!("Uncompatible dimensions!");
                }
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for j in 0..other.num_cols {
                    let other_column_j = &other.data[j * other.num_rows .. (j + 1) * other.num_rows];
                    for i in 0..self.num_rows {
                        data.push(
                            other_column_j.iter().enumerate().map(|(k, element)| *element * self.get(i, k)).sum()
                        );
                    }
                }
                return Matrix::new_from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
}
matrix_multiply!(&Matrix<T>, &Matrix<T>, T);
matrix_multiply!(Matrix<T>, Matrix<T>, T);
matrix_multiply!(&Matrix<T>, Matrix<T>, T);
matrix_multiply!(Matrix<T>, &Matrix<T>, T);