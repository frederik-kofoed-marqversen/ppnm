use super::Matrix;
use std::ops::{Add, Sub, Neg, AddAssign};
use std::iter::zip;

// ADD TWO MATRICES
macro_rules! matrix_add {
    ($LHS:ty, $RHS:ty, $T:tt) => {
        impl<$T: Add<Output = $T> + Copy> Add<$RHS> for $LHS {
            type Output = Matrix<$T>;
            fn add(self, other: $RHS) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for (a, b) in zip(self.iter(), other.iter()) {
                    data.push(*a + *b);
                }
                return Matrix::from_data(data, self.num_rows, self.num_cols);
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
    ($LHS:ty, $RHS:ty, $T:tt) => {
        impl<$T: Sub<Output = $T> + Copy> Sub<$RHS> for $LHS {
            type Output = Matrix<$T>;
            fn sub(self, other: $RHS) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for (a, b) in zip(self.iter(), other.iter()) {
                    data.push(*a - *b);
                }
                return Matrix::from_data(data, self.num_rows, self.num_cols);
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
    ($RHS:ty, $T:tt) => {
        impl<$T: Neg<Output = $T> + Copy> Neg for $RHS {
            type Output = Matrix<$T>;
            fn neg(self) -> Self::Output {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for scalar in self.iter() {
                    data.push(-*scalar);
                }
                return Matrix::from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
}
matrix_negate!(&Matrix<T>, T);
matrix_negate!(Matrix<T>, T);

impl<T: Copy + AddAssign<T>> AddAssign<Matrix<T>> for Matrix<T> {
    fn add_assign(&mut self, other: Matrix<T>) {
        for (item, scalar) in zip(self.iter_mut(), other.iter()) {
            *item += *scalar;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mat_neg() {
        let mat = Matrix::from_data(vec![1, 2, 3, 4], 2, 2);
        assert_eq!(
            -mat,
            Matrix::from_data(vec![-1, -2, -3, -4], 2, 2)
        );
    }

    #[test]
    fn test_mat_sub() {
        let mat1 = Matrix::from_data((1..=6).collect(), 3, 2);
        let mat2 = Matrix::from_data(vec![1; 6], 3, 2);
        assert_eq!(
            mat1 - mat2,
            Matrix::from_data((0..6).collect(), 3, 2)
        );
    }

    #[test]
    fn test_mat_add() {
        let mat1 = Matrix::from_data((0..6).collect(), 3, 2);
        let mat2 = Matrix::from_data(vec![1; 6], 3, 2);
        assert_eq!(
            mat1 + mat2,
            Matrix::from_data((1..=6).collect(), 3, 2)
        );
    }
}