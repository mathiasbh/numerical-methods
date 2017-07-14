Name: Mathias Bojsen-Hansen
Date: 28-06-2016 - 30-06-2016
--

Available projects: number 9, or 11.
Choosing project 9: Inverse iteration algorithm for eigenvalues (and eigenvectors)
--

EXERSICE A: Implement inverse iteration and test if it works for some matrix

EXERCISE B: Modify the previous code to Rayleigh quotient iteration with varying eigenvalue.

EXERCISE C: Test and compare effectiveness (steps and timing) of A and B, i.e
	    without and with changing the eigenvalue shift as the program runs.
	    See how the amount of steps increases for how bad the initial
	    eigenvalue guess is for the two methods.
--

The idea behind inverse iteration is to find an eigenvector and -value.
The algoritm will find the nearest associated eigenvector for a user specified
eigenvalue (approximately known) by solving the linear expression (A-s*I)b_(i+1) = b_(i)
where s, is the specified eigenvalue (constant in exercise A, variable in exercise B),
A is a matrix, and b is a vector.

Use the program by giving the function a matrix to solve, an approximate guess on an eigenvalue
and an initial guess on eigenvector. The program returns the eigenvalue while the eigenvector is
writes over the vector guess (passed through parameter).

Test using matrix A = (first matrix)
    1.5000    1.0000   -0.5000
    2.5000    0.7500   -1.2500
    1.5000    0.7500   -0.2500

with known eigenvalues and eigenvectors:
eigval = 0.5,  vec = ( 1  0  2) 
eigval = 2  ,  vec = ( 1  1  1)
eigval =-0.5,  vec = (-1  2  0)


Test using matrix B = (second matrix)
    2   3   1
    0   4   0
    2   0   1

with known eigenvalues and eigenvectors:
eigval = 4,  vec = ( 9  4  6)
eigval = 3,  vec = ( 1  0  1)
eigval = 0,  vec = (-1  0  2)

(OBS: The eigenvectors are not normalized to unity in result, but I will make sure they are correct.)

Note to OutputC_steps.pdf figure:
     The figure represents the number of steps to reach a certain tolerance for the two methods as a
     function of how far the guess is from the true eigenvalue.

     We see inverse iteration quickly jumps to the specified maximum of steps after about 200 shifts away from
     the true eigenvalue. By implementing rayleigh iteration the number of steps only increases from 33 to 50
     by going from true eigenvalue to 1000 shifts away. This method is much more efficient.

The two files print_functions.o and qrgivens.o are from folder E03_linearEquation
