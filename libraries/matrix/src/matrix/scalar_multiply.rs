use super::Matrix;
use std::ops::{Mul, MulAssign};

macro_rules! matrix_multiply_scalar_right {
    ($MatrixType:ty, $T:tt) => {
        impl<$T: Mul<Output = $T> + Copy> Mul<$T> for $MatrixType {
            type Output = Matrix<$T>;
            fn mul(self, scalar: $T) -> Matrix<$T> {
                let mut data: Vec<$T> = Vec::with_capacity(self.num_rows * self.num_cols);
                for element in self.iter() {
                    data.push(scalar * *element);
                }
                return Matrix::from_data(data, self.num_rows, self.num_cols);
            }
        }
    };
}
matrix_multiply_scalar_right!(Matrix<T>, T);
matrix_multiply_scalar_right!(&Matrix<T>, T);

macro_rules! matrix_multiply_scalar_left {
    ($MatrixType:ty, $T:ty) => {
        impl Mul<$MatrixType> for $T {
            type Output = Matrix<$T>;
            fn mul(self, matrix: $MatrixType) -> Matrix<$T> {
                matrix * self
            }
        }
    };
}
macro_rules! matrix_multiply_scalar_left_impl {
    ($($T: ty),*) => {$(
        matrix_multiply_scalar_left!(Matrix<$T>, $T);
        matrix_multiply_scalar_left!(&Matrix<$T>, $T);
    )*}
}
matrix_multiply_scalar_left_impl!(u8, u16, u32, u64, usize, i8, i16, i32, i64, isize, f32, f64);

impl<T: Copy + MulAssign<T>> MulAssign<T> for Matrix<T> {
    fn mul_assign(&mut self, scalar: T) {
        for item in self.iter_mut() {
            *item *= scalar;
        }
    }
}