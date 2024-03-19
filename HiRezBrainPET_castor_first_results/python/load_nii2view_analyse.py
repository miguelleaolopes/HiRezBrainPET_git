import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import pickle
import glob

# Path to the folder containing the nii files
nii_file_listmode = "./benchmark_pet_list-mode_challenger_it2.nii"
nii_file_sensitivity = "./benchmark_pet_list-mode_challenger_sensitivity.nii"

# Dictionary to store the data for each file
data_dict = {}

# Control to show the plots
show_plot_slices = False
show_plot_3D = False

# Loop over the files
for nii_file in [nii_file_listmode, nii_file_sensitivity]:
    print(f"Reading {nii_file}")

    # Load nii file
    nii_image = nib.load(nii_file)
    print(f"Header: {nii_image.header}")

    # Get the image data
    img_data = nii_image.get_fdata()
    print(f"Image data shape: {img_data.shape}")

    ####################################################
    if (show_plot_slices == True):
        # Create a plot
        fig, (ax1,ax2) = plt.subplots(1,2, figsize=(5,6), gridspec_kw={'width_ratios': [20,1]})
        plt.subplots_adjust(bottom=0.15, right=0.8, top=0.9, wspace=0.15) # Adjust the margin to make room for the slider

        # Plot the middle slice
        current_slice = img_data.shape[2]//2
        im = ax1.imshow(img_data[:,:,current_slice], cmap="gray")
        ax1.set_title(f"Slice {current_slice+1}/{img_data.shape[2]}")
        plt.colorbar(im,cax=ax2)

        # Create a slider to scroll through the images
        slice_slider_ax = plt.axes([0.1, 0.05, 0.8, 0.03]) # Define the position and size of the slider [left, bottom, width, height]
        slice_slider = Slider(slice_slider_ax, "Slice", 0, img_data.shape[2]-1, valinit=current_slice, valstep=1)

        # Function to update the image when the slider is moved
        def update(val):
            current_slice = int(slice_slider.val)
            im.set_data(img_data[:,:,current_slice])
            ax1.set_title(f"Slice {current_slice+1}/{img_data.shape[2]}")
            fig.canvas.draw_idle()

        # Connect the slider to the update function
        slice_slider.on_changed(update)

        # Show the plot
        plt.show()

    ####################################################
    if (show_plot_3D == True):
        # Create a 3D plot
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Determine the non-zero range for colormap normalization
        nonzero_data = img_data[np.nonzero(img_data)]
        vmax = np.max(nonzero_data)
        vmin = np.percentile(nonzero_data, 10)

        # Set the aspect ratio to 'equal' to avoid stretching
        ax.set_box_aspect([img_data.shape[1], img_data.shape[0], img_data.shape[2]])

        # Plot each slice individually
        for z in range(img_data.shape[2]):
            masked_img_data = np.where(img_data < vmin, np.nan, img_data)
            slice_data = masked_img_data[:, :, z]
            x, y = np.meshgrid(np.arange(slice_data.shape[1]), np.arange(slice_data.shape[0]))
            colors = plt.cm.gray((slice_data - vmin) / (vmax - vmin))  # Normalize grayscale values
            ax.plot_surface(x, y, np.full_like(slice_data, z), facecolors=colors, edgecolor='none')

        # Add a colorbar
        norm = plt.Normalize(vmin=vmin, vmax=vmax)
        sm = plt.cm.ScalarMappable(cmap='gray', norm=norm)
        cbar = fig.colorbar(sm, ax=ax, shrink=0.7)

        # Set the axis labels
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Slice')
        ax.set_title("3D image")

        # Show the plot
        plt.show()
    
    ####################################################
    # Store the data in the dictionary
    data_dict[nii_file] = {"nii_image": nii_image}

# Save the data dictionary to a file
with open("data_dict.pkl", "wb") as f:
    pickle.dump(data_dict, f)

####################################################
# This nii image has the information about the duration and start time of the acquisition !!! only if nFrames > 1
# => in the extensions of the header (nii_image.header.extensions[0].get_content()[(0x____, 0x____)].value)
# Print the extensions  to identify the number of the wanted extensions:
# print(nii_image.header.extensions)
# See also https://nipy.org/nibabel/nibabel_images.html#extensions