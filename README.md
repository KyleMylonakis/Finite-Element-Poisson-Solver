# Finite-Element-Poisson-Solver
Using Finite Element Methods to solve a two dimensional Poisson problem.

The writeup is the included PDF file. 

The extra files are the python code used to run the simulation. The main piece of code to run, which includes the algorithm creating the assembly matrix, is titled "compute mass and stiffness matrix.py".

The other python files create the triangulation data structure and define the generic basis elements; these files are simply imported in the usual way in python.

SciPy and NumPy are required.

Moreover the first for loop in the main file indexes over n, which is the inverse of the lattice constant. In stage n
in the for loop, it generates an approximate solution to the variational problem for a uniform triangulation of the square with 2n^2 triangles.

To see the exponent of the rate of error decrease as the mesh size is included have the for loop iterate over more than one element, and uncomment the code at the end of the main file and comment out the following code

ax = Axes3D(plt.gcf())
ax.scatter(lagrange_points[:,0],lagrange_points[:,1],full_solution,zdir = 'z')
plt.show()

located at lines 230-232 instead.
