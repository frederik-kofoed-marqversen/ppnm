pub struct MatrixIterator<'a, T> {
    col_index: usize,
    row_index: usize,
    index: usize,
    data: &'a [T],
    num_cols: usize,
    num_rows: usize,
}

impl<'a, T> MatrixIterator<'a, T> {
    pub fn new(data: &'a [T], num_rows: usize, num_cols: usize) -> MatrixIterator<'a, T> {
        MatrixIterator {
            col_index: 0, 
            row_index: 0, 
            index: 0, 
            data: data, 
            num_rows: num_rows, 
            num_cols: num_cols,
        }
    }
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