# TEST THE INTERSECTION OF A LINE WITH A CIRCLE PYTHON
# Import file with function
from intersect_line_circle import intersect_line_circle
import numpy as np
# Define number of lines/LORs to generate
n = 10
# Define radius of circle
r = 20
# Generate random coefficients for the line equation
A = np.random.uniform(-1, 1, n)
B = np.random.uniform(-1, 1, n)
C = np.random.uniform(-1, 1, n)

# Call function
for i in range(n):
    pts = intersect_line_circle(A[i], B[i], C[i], r, show=True,verbose=True)
    print(pts)