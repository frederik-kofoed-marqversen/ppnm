use super::Matrix;
use std::ops::{Add, Mul};
use std::iter::Sum;

impl<T: Copy + Sum + Mul<Output = T> + Add<Output = T>> Matrix<T> {
    pub fn mat_mul(lhs: &Self, rhs: &Self) -> Self {
        if lhs.num_cols != rhs.num_rows {
            panic!("Non-compatible dimensions!");
        }
        let mut data: Vec<T> = Vec::with_capacity(lhs.num_rows * rhs.num_cols);
        for j in 0..rhs.num_cols {
            let rhs_column_j = &rhs.data[j * rhs.num_rows .. (j + 1) * rhs.num_rows];
            for i in 0..lhs.num_rows {
                data.push(
                    rhs_column_j.iter().enumerate().map(|(k, element)| *element * lhs.get(i, k)).sum()
                );
            }
        }
        return Matrix::from_data(data, lhs.num_rows, rhs.num_cols);
    }

    pub fn pow(matrix: &Self, power: u32) -> Self {
        if power == 1 {return matrix.clone()}
        matrix * Self::pow(matrix, power - 1)
    }
}

macro_rules! matrix_multiply {
    ($LHS:ty, $RHS:ty, $T:tt) => {
        impl<$T: Add<Output = $T> + Mul<Output = $T> + Copy + Sum> Mul<$RHS> for $LHS {
            type Output = Matrix<$T>;
            fn mul(self, other: $RHS) -> Matrix<$T> {
                Matrix::<$T>::mat_mul(&self, &other)
            }
        }
    };
}
matrix_multiply!(&Matrix<T>, &Matrix<T>, T);
matrix_multiply!(Matrix<T>, Matrix<T>, T);
matrix_multiply!(&Matrix<T>, Matrix<T>, T);
matrix_multiply!(Matrix<T>, &Matrix<T>, T);



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mul_1() {
        let mat = Matrix::new(vec![vec![1, 3], vec![2, 4]]);
        assert_eq!(
            &mat * &mat,
            Matrix::from_data(vec![7, 15, 10, 22], 2, 2)
        );
    }

    #[test]
    fn test_mul_2() {
        let mat1 = Matrix::new(vec![vec![1, 1], vec![1, -1]]);
        let mat2 = Matrix::new(vec![vec![1, 0], vec![0, -1]]);
        assert_eq!(
            (&mat1 * &mat2) * &mat1,
            Matrix::from_data(vec![0, 2, 2, 0], 2, 2)
        );
    }

    #[test]
    fn test_mul_3() {
        let mat1 = Matrix::new(vec![vec![1, 2], vec![1, 2], vec![1, 2]]);
        let mat2 = Matrix::from_data((0..6).collect(), 3, 2);
        assert_eq!(
            mat1 * mat2,
            Matrix::from_data(vec![3, 6, 12, 24], 2, 2)
        );
    }

    #[test]
    #[should_panic(expected = "Non-compatible dimensions!")]
    fn test_mul_4() {
        let mat1 = Matrix::new(vec![vec![1, 2], vec![1, 2], vec![1, 2]]);
        let mat2 = Matrix::from_data((0..6).collect(), 2, 3);
        let _ = mat2 * mat1;
    }

    #[test]
    fn test_pow_1() {
        let mat = Matrix::new(vec![vec![1, 3], vec![2, 4]]);
        assert_eq!(
            Matrix::pow(&mat, 3),
            (&mat * &mat) * &mat
        );
    }
}