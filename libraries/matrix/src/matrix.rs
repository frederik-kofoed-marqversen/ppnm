mod constructors;
mod matrix_add;
mod matrix_multiply;
mod scalar_multiply;

#[derive(Debug, Clone, PartialEq)]
pub struct Matrix<T> {
    data: Vec<T>,
    pub num_cols: usize,
    pub num_rows: usize,
}

impl<T: Copy> Matrix<T> {
    pub fn new(mut mat: Vec<Vec<T>>) -> Self {
        let num_cols = mat.len();
        let num_rows = mat[0].len();
        let mut data = Vec::with_capacity(num_cols * num_rows);
        for col in mat.iter_mut() {
            while col.len() > 0 {
                data.push(col.remove(0));
            }
        }
        Self {
            data: data,
            num_cols: num_cols,
            num_rows: num_rows,
        }
    }

    pub fn from_data(data: Vec<T>, num_rows: usize, num_cols: usize) -> Self {
        if data.len() != num_cols * num_rows {panic!("Non-compatible dimensions!");}
        Self {
            data: data, 
            num_cols: num_cols, 
            num_rows: num_rows,
        }
    }

    pub fn get(&self, row_index: usize, col_index: usize) -> T {
        self.data[self.num_rows * col_index + row_index]
    }

    pub fn iter(&self) -> impl Iterator<Item=&T>{
        self.data.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut T> {
        self.data.iter_mut()
    }

    pub fn transpose(&self) -> Self {
        let mut data = Vec::with_capacity(self.num_cols * self.num_rows);
        for i in 0..self.num_rows {
            for j in 0..self.num_cols {
                data.push(self.get(i, j));
            }
        }
        return Self::from_data(data, self.num_cols, self.num_rows);
    }
}

use std::ops::{Index, IndexMut};
impl<T> Index<usize> for Matrix<T> {
    type Output = [T];
    
    fn index(&self, col_index: usize) -> &[T] {
        &self.data[col_index * self.num_rows..(col_index + 1) * self.num_rows]
    }
}

impl<T> IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, col_index: usize) -> &mut [T] {
        &mut self.data[col_index * self.num_rows..(col_index + 1) * self.num_rows]
    }
}

/* impl<T: Copy + std::fmt::Display> std::fmt::Display for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[\n")?;
        for i in 0..self.num_rows {
            write!(f, " [")?;
            for j in 0..self.num_cols {
                write!(f, "{},", self.get(i, j))?;
            }
            write!(f, "]\n")?;
        }
        write!(f, "]")
    }
} */

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index() {
        let mat = Matrix::from_data(vec![1, 2, 3, 4, 5, 6], 2, 3);
        assert_eq!(
            mat[1],
            [3, 4]
        );
    }

    /* #[test]
    fn test_index() {
        let mat = Matrix::from_data(vec![1, 2, 3, 4, 5, 6], 2, 3);
        assert_eq!(
            mat[1..4],
            [2, 3, 4]
        );
    }

    #[test]
    fn test_index_mut() {
        let mut mat = Matrix::from_data(vec![1, 2, 3, 4, 5, 6], 2, 3);
        mat[1..4] = [0; 3];
        assert_eq!(
            mat[..],
            [1, 0, 0, 0, 5, 6]
        );
    } */
}