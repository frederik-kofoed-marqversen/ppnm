mod matrix;
use matrix::Matrix;

fn main() {
    let mat1 = Matrix::new(vec![vec![0, 2], vec![1, 3]]);
    let mat2 = Matrix::new(vec![vec![0, 2], vec![1, 3]]);
    dbg!(mat1 * mat2);
    let mat = Matrix::ones(3, 3);
    dbg!(&mat * 2.0);
    //dbg!(2.0 * &mat);
}