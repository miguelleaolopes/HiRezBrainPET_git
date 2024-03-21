# Load generic libraries
import numpy as np
from tqdm import tqdm
# Load specific libraries
import h5py
import datetime
from fractions import Fraction
from rich.console import Console
# Load functions of intersection
from line_intersection.intersect_line_cylinder import intersect_line_cylinder

### Parameters that can be changed
# Path to input file - mat_file
# Path to output files - file_diff_msec, file_diff_msec_h5, file_datetime
# Fraction of events to choose - nvalue
# Option to save files - save
###

##############################################
# Load data from .mat file and save data. Deprecated if LoadSave_TimeLORs.py is used
##############################################

mat_file = '../Files/220409181025_4997kLORs.mat'
file = h5py.File(mat_file, 'r')
Console(style='bold green').print('File {} loaded successfully'.format(mat_file))

# Print the keys of the file
Console(style='purple').print('List of keys in the file:', list(file.keys()))

# Load the data from all keys
data = [file[key] for key in file.keys()]
# Set the data from the keys
EvtTime, LORs, RunName, opt = data[1][()], data[2][()], data[3][()], data[4]

##############################################
# Set times of events
##############################################

# Check if conversion is needed
file_diff_msec = '../Files/time_diff_msec.npy'
file_diff_msec_h5 = '../Files/time_diff_msec.h5'
file_datetime = '../Files/time_datetime.npy'
try:
    time_diff_msec = np.load(file_diff_msec, allow_pickle=True)
    time_datetime = np.load(file_datetime, allow_pickle=True)
    Console(style='purple').print('File {} exists'.format(file_diff_msec))
    Console(style='purple').print('File {} exists'.format(file_datetime))
except FileNotFoundError:
    Console(style='red').print('File {} does not exist'.format(file_diff_msec))

    Console(style='purple').print('Converting time...')
    # Remove the extra dimension
    time_full = np.squeeze(EvtTime)

    # Convert the event time to a datetime object
    time_datetime = [datetime.datetime.fromordinal(int(t)) + datetime.timedelta(days=t % 1) - datetime.timedelta(days=366) for t in tqdm(time_full, desc='Converting to datetime')]
    # Time diference between the first and each event
    time_diff = [t - time_datetime[0] for t in time_datetime]

    # Difference in seconds
    time_diff_sec = [t.total_seconds() for t in time_diff]

    # Difference in milliseconds
    time_diff_msec = [t.total_seconds() * 1000 for t in time_diff]

    # Save the time difference in a files
    np.save(file_diff_msec, time_diff_msec)
    np.save(file_datetime, time_datetime)
    # Save the time difference in a hdf5 file
    hf = h5py.File(file_diff_msec_h5, 'w')
    hf.create_dataset('time_diff_msec', data=time_diff_msec)
    hf.close()

## Print the time of the first event
Console(style='bold').print('Time of the first event:', time_datetime[0])
Console(style='bold').print('Time of the last event: ', time_datetime[-1])
Console(style='bold').print('Time elapsed:           ', time_datetime[-1] - time_datetime[0])

##############################################
# Set LORs
##############################################

# Shape of LORs
Console().print('Shape of LORs [(Pm,Dn); (x,y,z); nEvents]:', LORs.shape)

nEvents = LORs.shape[2]

# Since there is a lot of Events where are gone chose only a fraction 
nvalue = 4
fract = Fraction(1, nvalue)
# First fraction of events
Console().print('Choosing first {} events...'.format(str(fract)))

# Loop over the fractions
for i in range(nvalue):
    Console().print('Fraction {}:'.format(i+1), int(nEvents*fract*i), 'to', int(nEvents*fract*(i+1)))
    # list of indexs for each fraction
    indexs = np.arange(int(nEvents*fract*i), int(nEvents*fract*(i+1)))

    print('Number of events chosen:', len(indexs))

    # Getting points along the LORs of each event but only for the chosen events, since there is a lot of them
    LORs_small = LORs[:,:,indexs]
    time_diff_msec_small = np.array(time_diff_msec)[indexs]
    Console().print('Shape of LORs_small:', LORs_small.shape)
    nEvents_small = LORs_small.shape[2]

    # Get the points of intersection of the LORs with the cylinder of radius r and length l
    Console().print('Getting the points of intersection of the LORs with the cylinder of radius r and length l...')
    Console().print('Number of points before intersection:', nEvents_small)
    intersects_all = [intersect_line_cylinder(LORs_small[0,:,i], LORs_small[1,:,i], radius=155, height=300, show=False, verbose=False) for i in tqdm(range(nEvents_small), desc='Intersecting LORs with cylinder')]

    # Find valid intersections 
    valid_indices = [i for i, x in enumerate(intersects_all) if x is not None and len(x) == 2]
    pts_intersect_small = np.array([intersects_all[i] for i in valid_indices])
    times_intersect_small = np.array([time_diff_msec_small[i] for i in valid_indices])
    Console().print('Number of points after intersection (removing no intersections):', pts_intersect_small.shape)


    ##############################################
    # Save the points of intersection
    ##############################################

    save = True
    if save:    
        # Save the points of intersection along with time, x, y, and z coordinates to a text file
        data_to_save = np.column_stack((times_intersect_small, pts_intersect_small[:, 0, 0], pts_intersect_small[:, 0, 1], pts_intersect_small[:, 0, 2], pts_intersect_small[:, 1, 0], pts_intersect_small[:, 1, 1], pts_intersect_small[:, 1, 2]))
        # Change name depending on the fraction of number of points
        np.savetxt('../Files/intersect_{}th_points_with_time.txt'.format(i+1), data_to_save, header='time (msec), x1, y1, z1, x2, y2, z2', delimiter=',')
        Console(style='bold green').print('File saved: intersect_{}th_points_with_time.txt'.format(i+1))


        # # Create an HDF5 file and store the data
        # with h5py.File('../Files/intersect_{}th_points_with_time.format(i+1).h5', 'w') as hf:
        #     # Store times_intersect_small as a dataset
        #     hf.create_dataset('times', data=times_intersect_small)

        #     # Store pts_intersect_small as a dataset
        #     # hf.create_dataset('points', data=pts_intersect_small)
        #     hf.create_dataset('points_xm', data=pts_intersect_small[:, 0, 0])
        #     hf.create_dataset('points_ym', data=pts_intersect_small[:, 0, 1])
        #     hf.create_dataset('points_zm', data=pts_intersect_small[:, 0, 2])
        #     hf.create_dataset('points_xd', data=pts_intersect_small[:, 1, 0])
        #     hf.create_dataset('points_yd', data=pts_intersect_small[:, 1, 1])
        #     hf.create_dataset('points_zd', data=pts_intersect_small[:, 1, 2])
        # Console(style='bold green').print('File saved: intersect_{}th_points_with_time.format(i+1).h5')
