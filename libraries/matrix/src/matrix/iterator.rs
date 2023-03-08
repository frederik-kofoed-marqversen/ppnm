/* use super::Matrix;

impl<T> Matrix<T> {
    pub fn iter<'a>(&'a self) -> MatrixIterator<'a, T> {
        MatrixIterator{0, 0, 0, &self.data, self.num_rows, self.num_cols}
    }
    
    pub fn iter_mut<'a>(&'a mut self) -> MatrixIteratorMut<'a, T> {
        MatrixIteratorMut{0, 0, &self.data, self.num_rows, self.num_cols}
    }
} */

pub struct MatrixIterator<'a, T> {
    col_index: usize,
    row_index: usize,
    index: usize,
    data: &'a [T],
    num_cols: usize,
    num_rows: usize,
}

impl<'a, T> Iterator for MatrixIterator<'a, T> {
    type Item = (&'a T, usize, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.row_index == self.num_rows {
            self.col_index += 1;
            self.row_index = 0;
        }
        if self.col_index == self.num_cols {
            return None
        } else {
            let result = Some((&self.data[self.index], self.row_index, self.col_index));
            self.index += 1;
            self.row_index += 1;
            return result
        }
    }
}

pub struct MatrixIteratorMut<'a, T> {
    col_index: usize,
    row_index: usize,
    data: &'a [T],
    num_cols: usize,
    num_rows: usize,
}

impl<'a, T: Copy> std::iter::Iterator for MatrixIteratorMut<'a, T> {
    type Item = (&'a mut T, usize, usize);
    fn next(&mut self) -> Option<Self::Item> {
        if self.col_index == self.num_cols {
            self.row_index += 1;
            self.col_index = 0;
        }
        return if self.row_index == self.num_rows {
            None
        } else {
            let col_index = self.col_index;
            self.col_index += 1;
            let data = std::mem::replace(&mut self.data, &mut []);
            if let Some((v, rest)) = data.split_first_mut() {
            self.data = rest;
            Some((v, self.row_index, col_index))
            } else {
            None
            }
        };
    }
}