The writeup is the included PDF file. 

The extra files are the python code used
 to run the simulation. The main piece of code with the assembly algorith is 
titled compute mass and stiffness matrix.py. 


The other python files create the triangulation data structure and define the 
generic basis elements; these files are simply imported in the usual way in
 python



If you want to run the file you will need to uncomment some print and plot
 statements in the main file. 
Moreover the first for loop in the main file 
indexes over the inverse of the lattice constant, in the sense that at stage n
in the for loop, 
it generates an approximate solution to the variational 
problem for a uniform triangulation of the square with 2n^2 triangles.

