import numpy as np
from math import sqrt

lagrange_test_points = [ [0,0], [0.5,0] ,[1,0], [0.5,0.5],[0,1],[0,0.5] ]
coefficients = []

for i in range(0,6):
    
    b = [] # Initialize the b vector
    for j in range (0,6):
        if j == i:
            b.append(1)
        else:
            b.append(0)
            
    b = np.array(b)
    coefficient_matrix = []
    for j in range(0,6):
        row = []
        row.append(1)
        row.append(lagrange_test_points[j][0])
        row.append(lagrange_test_points[j][1])
        row.append(lagrange_test_points[j][0]*lagrange_test_points[j][1])
        row.append(lagrange_test_points[j][0]**2)
        row.append(lagrange_test_points[j][1]**2)
        coefficient_matrix.append(row)
    coefficient_matrix =np.array(coefficient_matrix)
    #print coefficient_matrix
    
    coefficients.append(np.linalg.solve(coefficient_matrix,b))

coefficients = np.array(coefficients)
# print coefficients

# Define the shape functions on arbitrary elements
def shape_function(x,y,n):
    return coefficients[n][0] + coefficients[n][1]*x + coefficients[n][2]*y + coefficients[n][3]*x*y + coefficients[n][4]*x*x + coefficients[n][5]*y*y
    
def grad_shape_function(x,y,n):
    gradient = np.array([ coefficients[n][1] + coefficients[n][3]*y + 2*x*coefficients[n][4], coefficients[n][2] + coefficients[n][3]*x + 2*y*coefficients[n][5] ])
    return gradient
      