mod matrix_f64;
mod iterator;
mod matrix_add;
mod matrix_multiply;
mod scalar_multiply;

#[derive(Debug, Clone)]
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

    fn new_from_data(data: Vec<T>, num_rows: usize, num_cols: usize) -> Self {
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

    pub fn iter<'a>(&'a self) -> iterator::MatrixIterator<'a, T> {
        iterator::MatrixIterator::new(&self.data, self.num_rows, self.num_cols)
    }

    pub fn transpose(&self) -> Self {
        let mut data = Vec::with_capacity(self.num_cols * self.num_rows);
        for i in 0..self.num_rows {
            for j in 0..self.num_cols {
                data.push(self.get(i, j));
            }
        }
        return Self::new_from_data(data, self.num_cols, self.num_rows);
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