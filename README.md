# MPI_2D_Laplace_Solver
Two dimensional laplace equation solver using MPI
# Details
Uses an n x n matrix (n must be even), filled with randomized doubles from 0-49. 
Only uses 4 processes.
Uses a five-point stencil, assigning average of four neighbors (up,down,left,right) to each element. If element lies on the border of the matrix and is missing neighbors, the missing neighbors are considered as 25. Calculates average difference betweeen iterations, running until difference is less than 0.1.


