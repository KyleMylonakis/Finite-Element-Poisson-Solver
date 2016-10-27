import numpy as np

def InverseMeshConstant(n):
    
    # Pick any integer greater than one. The total number of triangles will be 2n^2
    h = 1/float(n) # Lattice constant (well up to scaling by root 2)

    verticies = []
    triangulation = []
    boundary = []
    
    for i in range(0,n+1):
        for j in range(0,n+1):
            verticies.append([h*j,h*i])
    verticies  = np.array(verticies)   
    #print verticies
    
    # Create the triangulation
    for i in range(0,n):
        for j in range(0,n):
            triangulation.append([verticies[(n+1)*i+j],verticies[(n+1)*i+j+1],verticies[(n+1)*i+j+n+1]])
            triangulation.append([verticies[(n+1)*i+j+1],verticies[(n+1)*i+j+n+2],verticies[(n+1)*i+j+n+1]])
    
    # This tells us if the k-th vertex is on the boundary            
   # for i in verticies:
   #     if i[0] == 1 or i[1] == 1 or i[0] == 0 or i[1] == 0:
   #         boundary.append(1)
   #     else:
   #         boundary.append(0)
    
    #boundary = np.array(boundary)                              
    triangulation = np.array(triangulation)
    #print len(triangulation)
   
    return [triangulation, verticies, boundary] 
    
    
    #print triangulation
