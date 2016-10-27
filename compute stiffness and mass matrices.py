import numpy as np
import triangulate_a_square as tri
import compute_generic_basis_quad_2d_lagrange as gen
from scipy.integrate import dblquad
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import log
import matplotlib.patches as mpatches
fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

elltwo_error_box = []
energy_error_box = []
H1_error_box = []
sup_error_box = []
node_spacing_box = []


def forcing_function(x,y):
    return 200*(x**2 + y**2 -x - y)
def analytic_soln(x,y):
    return 100*(x**2 - x)*(y**2-y)




for n in range(10,11):
    #n=2 # Set as an integer
    #print tri.InverseMeshConstant(n)[0]
    
    lagrange_points = []        
    
    # Initialize the mass and stiffness matrix    
    mass_matrix = []
    stiffness_matrix = []
    
    for i in range(0,4*n**2 + 4*n +1):
        row = []
        for j in range(0,4*n**2 + 4*n +1):
            row.append(0)
        mass_matrix.append(row)
        stiffness_matrix.append(row)
            
            
                    
    for triangle in tri.InverseMeshConstant(n)[0]:
            
        # Compute the linear mapping part of the affine mapping taking our generic triangle to the informed triangle
        jacobian = np.array( [ [ triangle[1][0] - triangle[0][0], triangle[2][0] - triangle[0][0] ] , [ triangle[1][1] - triangle[0][1], triangle[2][1] - triangle[0][1] ] ] )
        #print jacobian
        
        # Compute where the lagrange test points are sent by the affine mapping
        
        for test_points in gen.lagrange_test_points:
            lagrange_points.append( np.transpose(np.add(np.dot(jacobian, np.transpose(test_points)) , np.transpose([triangle[0][0], triangle[0][1]] ) ) ))
    
    lagrange_points = np.array(lagrange_points)
    actual_lagrange_points = map(tuple,lagrange_points)
    lagrange_points = np.array(list(set(actual_lagrange_points)))
    #print lagrange_points
    
    
    # This tells us if the k-th lagrange point is on the boundary            
    
    lagrange_boundary = []
    for i in lagrange_points:
        if i[0] == 1 or i[1] == 1 or i[0] == 0 or i[1] == 0:
            lagrange_boundary.append(1)
        else:
            lagrange_boundary.append(0)

    # Now we begin the assembly procedure
    
    for triangle in tri.InverseMeshConstant(n)[0]:
        
        # Compute the linear mapping part of the affine mapping taking our generic triangle to the informed triangle
        jacobian = np.array( [ [ triangle[1][0] - triangle[0][0], triangle[2][0] - triangle[0][0] ] , [ triangle[1][1] - triangle[0][1], triangle[2][1] - triangle[0][1] ] ] )
        #print jacobian
        
        # Compute where the lagrange test points are sent by the affine mapping
        
        transformed_test_points = []
        for test_points in gen.lagrange_test_points:
            transformed_test_points.append( np.transpose(np.add(np.dot(jacobian, np.transpose(test_points)) , np.transpose([triangle[0][0], triangle[0][1]] ) ) ))
        
        # determine determinant of the Jacobian of the affine mapping
        scalar_jacobian = abs(np.linalg.det(jacobian))
        
        # Invert the jacobian of the affine mapping
        inverse_jacobian = np.linalg.inv(jacobian)
        
        
        # Compute the Mass submatrix for the triangle
        mass_submatrix = []
        for i in range(0,6):
            row = []
            for j in range(0,6):
                row.append( dblquad(lambda x,y: gen.shape_function(x,y,i) * gen.shape_function(x,y,j) * scalar_jacobian, 0 , 1, lambda x:0, lambda x: 1-x )[0] )
            mass_submatrix.append(row)
            
        mass_submatrix = np.array(mass_submatrix)
        #print mass_submatrix
        for i in range(0,6):
            for j in range(0,6):
                #lagrange_points = list(set(map(tuple,lagrange_points)))
                lagrange_points = lagrange_points.tolist()
                temp_i = np.transpose(np.add(np.dot(jacobian, np.transpose(gen.lagrange_test_points[i])) , np.transpose([triangle[0][0], triangle[0][1]]) ) ).tolist()
                temp_j = np.transpose(np.add(np.dot(jacobian, np.transpose(gen.lagrange_test_points[j])) , np.transpose([triangle[0][0], triangle[0][1]]) ) ).tolist()
                #print lagrange_points.index(temp_j)
                if lagrange_boundary[lagrange_points.index(temp_i)] == 0 and lagrange_boundary[lagrange_points.index(temp_j)] == 0:
                    mass_matrix[lagrange_points.index(temp_i)][lagrange_points.index(temp_j)] = mass_matrix[lagrange_points.index(temp_i)][lagrange_points.index(temp_j)] + mass_submatrix[i][j]
                lagrange_points = np.array(lagrange_points)
                
        mass_matrix = np.array(mass_matrix)
        #print mass_matrix
                
        # Compute the Stiffness submatrix for the triangle
        stiffness_submatrix = []
        for i in range(0,6):
            row = []
            for j in range(0,6):
                row.append( dblquad(lambda x,y:np.dot( np.dot(gen.grad_shape_function(x,y,i),inverse_jacobian ), np.dot(gen.grad_shape_function(x,y,j),inverse_jacobian ))* scalar_jacobian, 0 , 1, lambda x:0, lambda x: 1-x )[0] )
            stiffness_submatrix.append(row)
            
        stiffness_submatrix = np.array(stiffness_submatrix)
        #print stiffness_submatrix
        for i in range(0,6):
            for j in range(0,6):
                #lagrange_points = list(set(map(tuple,lagrange_points)))
                lagrange_points = lagrange_points.tolist()
                temp_i = np.transpose(np.add(np.dot(jacobian, np.transpose(gen.lagrange_test_points[i])) , np.transpose([triangle[0][0], triangle[0][1]]) ) ).tolist()
                temp_j = np.transpose(np.add(np.dot(jacobian, np.transpose(gen.lagrange_test_points[j])) , np.transpose([triangle[0][0], triangle[0][1]]) ) ).tolist()
                #print lagrange_points.index(temp_j)
                if lagrange_boundary[lagrange_points.index(temp_i)] == 0 and lagrange_boundary[lagrange_points.index(temp_j)] == 0:
                    stiffness_matrix[lagrange_points.index(temp_i)][lagrange_points.index(temp_j)] = stiffness_matrix[lagrange_points.index(temp_i)][lagrange_points.index(temp_j)] + stiffness_submatrix[i][j]
                lagrange_points = np.array(lagrange_points)
                
        stiffness_matrix = np.array(stiffness_matrix)
        #print stiffness_matrix
    
    # We need to remove the rows and columns of the stiffness matrix which are all zero (these correspond to boundary nodes)
    
    blank_row = []
    for i in range(0,4*n**2 + 4*n +1):
        blank_row.append(0)
    
    temp_stiffness_matrix = []
    actual_stiffness_matrix = []
    temp_mass_matrix = []
    actual_mass_matrix = []
    
    for i in range(0,4*n**2 + 4*n + 1):
        if  not (np.array_equal(stiffness_matrix[i] , blank_row)):
            temp_stiffness_matrix.append(stiffness_matrix[i])
        if not (np.array_equal(mass_matrix[i], blank_row)):
            temp_mass_matrix.append(mass_matrix[i])
    
    temp_stiffness_matrix = np.transpose(np.array(temp_stiffness_matrix))
    temp_mass_matrix = np.transpose(np.array(temp_mass_matrix))
    
    blank_row = []
    for j in range(0,len(temp_stiffness_matrix[0])):
        blank_row.append(0)
    
    for j in range(0,4*n**2 + 4*n +1):
        #print temp_stiffness_matrix[:,j]
        if not (np.array_equal(temp_stiffness_matrix[j], blank_row)):
            actual_stiffness_matrix.append(temp_stiffness_matrix[j])
        if not (np.array_equal(temp_mass_matrix[j], blank_row)):
            actual_mass_matrix.append(temp_mass_matrix[j])
    
    actual_mass_matrix = np.array(actual_mass_matrix)
    actual_stiffness_matrix = np.array(actual_stiffness_matrix)
    
            
    
    # Now we solve the linear algebra system Ku = Mf
    
    # Initialize the forcing term f
    forcing = []
    nonzero_lagrange_points = []
    for i in range(0,len(lagrange_points)):
        if lagrange_boundary[i] == 0:
            forcing.append(-forcing_function(lagrange_points[i][0],lagrange_points[i][1]))
            nonzero_lagrange_points.append(lagrange_points[i])
    lagrange_points = np.array(lagrange_points)
    nonzero_lagrange_points = np.array(nonzero_lagrange_points)
    
    #print actual_stiffness_matrix
    
    # Solve for u
    solution = np.linalg.solve(actual_stiffness_matrix,np.dot(actual_mass_matrix,forcing))
    #print solution
    
    full_solution = []
    very_temp = 0
    for i in range(0,len(lagrange_points)):
        if lagrange_boundary[i] == 1:
            full_solution.append(0)
        else:
            full_solution.append(solution[very_temp])
            very_temp = very_temp +1
    full_solution = np.array(full_solution)
    del very_temp
    
    #print solution
    #print full_solution
    
    # I need to find the L2 Error:
    elltwo_error = 0
    
    for i in range(0,len(lagrange_points)):
        elltwo_error += ((analytic_soln(lagrange_points[i][0],lagrange_points[i][1]) - full_solution[i])**2/float(2*n**2))
    elltwo_error = elltwo_error**0.5
    elltwo_error_box.append(elltwo_error)
    print "The L2 Error is ", elltwo_error
        
    
    # Sup Error
    sup_bowl = []
    sup_error = 0
    for i in range(0,len(nonzero_lagrange_points)):
        sup_bowl.append(abs(analytic_soln(nonzero_lagrange_points[i][0],nonzero_lagrange_points[i][1]) - solution[i]))
    sup_error = max(sup_bowl)
    sup_error_box.append(sup_error)
    print "The sup Error is ", sup_error
    
    node_spacing_box.append(1/float(n))

ax = Axes3D(plt.gcf())
ax.scatter(lagrange_points[:,0],lagrange_points[:,1],full_solution,zdir = 'z')
plt.show()

log_node_spacing_box = []
for i in node_spacing_box:
    log_node_spacing_box.append(-log(i))
    
log_elltwo_error_box = []
for i in elltwo_error_box:
    log_elltwo_error_box.append(log(i))
    
log_sup_error_box = []
for i in sup_error_box:
    log_sup_error_box.append(log(i))

log_energy_error_box = []
for i in energy_error_box:
    log_energy_error_box.append(log(i))
    

log_H1_error_box = []
for i in H1_error_box:
    log_H1_error_box.append(log(i))

sup_rate = np.polyfit(np.array(log_node_spacing_box), np.array(log_sup_error_box), 1)[0]  
elltwo_rate = np.polyfit(np.array(log_node_spacing_box), np.array(log_elltwo_error_box), 1)[0]

    
#print sup_rate
#print elltwo_rate

#plt.plot(log_node_spacing_box,log_sup_error_box)    
#plt.plot(log_node_spacing_box,log_elltwo_error_box)

#blue_patch = mpatches.Patch(color='blue', label='Sup Error')
#green_patch = mpatches.Patch(color='green', label = 'L2 Error')
#plt.legend(handles=[blue_patch, green_patch])

#plt.show()

#plt.xlabel('log(1/n)')
#plt.ylabel('log(error)')
#plt.title('log-log Error Plot')


#plt.show()