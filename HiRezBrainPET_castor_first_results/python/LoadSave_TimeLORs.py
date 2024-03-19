# Load generic libraries
import numpy as np
from tqdm import tqdm
# Load specific libraries
import h5py
import datetime
from fractions import Fraction

##############################################
# Load data from .mat file and save data
##############################################

mat_file = '../Files/220409181025_4997kLORs.mat'
file = h5py.File(mat_file, 'r')

# Save the data from all keys
data = []
for key in file.keys():
    print(key)
    print(file[key])
    data.append(file[key])

# Save the data from the keys
EvtTime = data[1][()]
LORs = data[2][()]
RunName =data[3][()]
opt = data[4]
##############################################
# Set times of events
##############################################
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

###################################################
# Save the time difference in a npy file and hdf5
###################################################
np.save('../Files/time_diff_msec.npy', time_diff_msec)
np.save('../Files/time_diff_sec.npy', time_diff_sec)
np.save('../Files/time_datetime.npy', time_datetime)

hf = h5py.File('../Files/time_diff_msec.h5', 'w')
hf.create_dataset('time_diff_msec', data=time_diff_msec)
hf.close()
hf = h5py.File('../Files/time_diff_sec.h5', 'w')
hf.create_dataset('time_diff_sec', data=time_diff_sec)
hf.close()

###################################################
# Save the LORs in a npy file and hdf5
###################################################
np.save('../Files/LORs.npy', LORs)
hf = h5py.File('../Files/LORs.h5', 'w')
hf.create_dataset('LORs', data=LORs)
hf.close()