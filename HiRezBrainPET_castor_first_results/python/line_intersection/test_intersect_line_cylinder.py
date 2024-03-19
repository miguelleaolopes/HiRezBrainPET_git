
# TEST THE INTERSECTION OF A LINE WITH A CYLINDER PYTHON
# Import file with function
from intersect_line_cylinder import intersect_line_cylinder
from intersect_line_circle import intersect_line_circle
import numpy as np
import h5py
import tqdm

# Define number of lines/LORs to generate
n = 10
# Define radius of cylinder
r = 155
# Define length of cylinder
l = 1000

mat_file = '../Files/220409181025_4997kLORs.mat'
file = h5py.File(mat_file, 'r')
data = []
for key in file.keys():
    data.append(file[key])
LORs = data[2][()]

print('Shape of LORs:', LORs.shape)

# Define points for the line
p = LORs[0,:,:].T
# Define direction for the line
d = LORs[1,:,:].T

# Get A, B, C from line equation in 2D
A = d[:,1]
B = -d[:,0]
C = -A*p[:,0] - B*p[:,1]

# Call function
show_plot = False
verbose_comment = False
count = 0
for i in tqdm.tqdm(range(int(A.shape[0]/100))):
# for i in [0,1]:
    pts = intersect_line_cylinder(p[i], d[i], r, l, show=show_plot,verbose=verbose_comment)
    pts2d = intersect_line_circle(A[i], B[i], C[i], r, show=show_plot,verbose=verbose_comment)
    if pts is None:
        # print('No intersection for line', i)
        count += 1

# pts = intersect_line_cylinder([0,0,0], [1,1,1], np.sqrt(2), 4, show=show_plot,verbose=verbose_comment)
# pts = intersect_line_cylinder(p[1], d[1], 155, 600, show=show_plot,verbose=verbose_comment)

# print('Point of intersection 3D:', pts)
# print('Point of intersection 2D:', pts2d)

print('Number of lines with no intersection:', count)