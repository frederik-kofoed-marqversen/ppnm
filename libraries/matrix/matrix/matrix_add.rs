use super::Matrix;
use std::ops::{Add, Sub, Neg};
use std::iter::zip;

// ADD TWO MATRICES
macro_rules! matrix_add {
    ($LHS:ty, $RHS:ty, $T:tt ) => {
        impl<$T: Add<Output = $T> + Copy> Add<$RHS> for $LHS {
            type Output = Matrix<$T>;
            fn add(self, other: $RHS) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for (a, b) in zip(self.iter(), other.iter()) {
                    data.push(*a.0 + *b.0);
                }
                return Matrix::new_from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
}
matrix_add!(&Matrix<T>, &Matrix<T>, T);
matrix_add!(Matrix<T>, Matrix<T>, T);
matrix_add!(&Matrix<T>, Matrix<T>, T);
matrix_add!(Matrix<T>, &Matrix<T>, T);


// SUBTRACT TWO MATRICES
macro_rules! matrix_subtract {
    ($LHS:ty, $RHS:ty, $T:tt ) => {
        impl<$T: Sub<Output = $T> + Copy> Sub<$RHS> for $LHS {
            type Output = Matrix<$T>;
            fn sub(self, other: $RHS) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for (a, b) in zip(self.iter(), other.iter()) {
                    data.push(*a.0 - *b.0);
                }
                return Matrix::new_from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
}
matrix_subtract!(&Matrix<T>, &Matrix<T>, T);
matrix_subtract!(Matrix<T>, Matrix<T>, T);
matrix_subtract!(&Matrix<T>, Matrix<T>, T);
matrix_subtract!(Matrix<T>, &Matrix<T>, T);

// NEGATE A MATRIX
macro_rules! matrix_negate {
    ($RHS:ty, $T:tt ) => {
        impl<$T: Neg<Output = $T> + Copy> Neg for $RHS {
            type Output = Matrix<$T>;
            fn neg(self) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for (scalar, _, _) in self.iter() {
                    data.push(-*scalar);
                }
                return Matrix::new_from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
}
matrix_negate!(&Matrix<T>, T);
matrix_negate!(Matrix<T>, T);