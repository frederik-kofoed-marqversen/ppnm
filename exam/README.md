I have been appointed problem 21: "Hessenberg factorization of a real square matrix using Jacobi transformations".

In Derivations.pdf I derive the formula for the angle used for each of the Jacobian rotation matrices and also show how the algorithm correctly and in a finite number of steps computes the Hessenberg decomposition of a general matrix. Also I derive a recursive formula for the determinant of a Hessenberg matrix and compare the scaling with the  and with QR factorisation.

My implementation of the algorithm is in the file located at "../libraries/matrix/src/linalg/eig.rs". In the current folder I test the method on a randomly generated 5x5 matrix. The results of the test are put into "Test.txt".

Finally I run the method on a number of different sized matrices and time the execution. The timing is performed in parallel by sending the jobs into the background. For comparison I also time my QR factorisation algorithm on the same set of matrices (RNG seeds are the same). The results are plotted in "Plot.svg".

My self-evaluation of this exam is 10 out of 10.