from scipy.cluster.vq import kmeans
import numpy as np


#it's two "triangles" connected by a single edge, clustered probably because of 
# a similar connection between groups, but it produces two different outputs
matrix = [[0, 1, 1, 0, 0, 0],
          [1, 0, 1, 0, 0, 0],
          [1, 1, 0, 1, 0, 0],
          [0, 0, 1, 0, 1, 1],
          [0, 0, 0, 1, 0, 1],
          [0, 0, 0, 1, 1, 0]]
 
for i in range(10):
    print kmeans(np.array(matrix), 2)
