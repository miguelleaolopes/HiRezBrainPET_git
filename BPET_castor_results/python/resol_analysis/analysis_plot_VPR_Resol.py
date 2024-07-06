import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pandas as pd
from itertools import cycle
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# # Load the data
# data = pd.read_csv('example/output_1test_vs01.csv')

# # Extracting rod sizes
# rod_sizes = [1.1, 1.0, 0.95, 0.9, 0.85, 0.8]

# # Setting up two subplots: one for VPR and one for Resolution
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))

# # Looping through each rod size to plot VPR and Resolution on separate subplots
# for rod in rod_sizes:
#     vpr_column = f"VPR {rod}"
#     resol_column = f"Resol {rod}"

#     # Plotting VPR on the first subplot
#     ax1.plot(data['Iteration'], data[vpr_column], marker='o', label=f"VPR {rod}")

#     # Plotting Resolution on the second subplot
#     ax2.plot(data['Iteration'], data[resol_column], marker='o', label=f"Resol {rod}")

# # Setting titles, labels, and layout for VPR plot
# ax1.set_title('VPR Values by Rod Size Across Iterations')
# ax1.set_xlabel('Iteration')
# ax1.set_ylabel('VPR')
# ax1.grid(True)
# ax1.legend(loc='upper right')

# # Setting titles, labels, and layout for Resolution plot
# ax2.set_title('Resolution Values by Rod Size Across Iterations')
# ax2.set_xlabel('Iteration')
# ax2.set_ylabel('Resolution')
# ax2.grid(True)
# ax2.legend(loc='upper right')

# # Adjust layout for clarity
# plt.tight_layout()

# # Show the plots
# plt.show()
def plot_vpr_resol_it_subsets(files, labels, subsets, rod_sizes, figsize=(14, 8), fontsize=12, legend_loc = ["lower left", "upper left"], save_bool=False, save_path='output'):

    plt.figure(figsize=figsize)
    markers = cycle(['o', 'x', '^', 's', 'p', '*', 'D', 'h'])
    colors = cycle(['black', 'red', 'green', 'blue', 'orange', 'magenta'][:len(rod_sizes)])

    ax_vpr = plt.subplot(2, 1, 1)
    ax_resol = plt.subplot(2, 1, 2)

    legend_handles_colors = [] 
    legend_handles_files = []  

    for file_path, subsets, label, marker in zip(files, subsets, labels, markers):
        data = pd.read_csv(file_path)
        data['Iteration'] = data['Iteration'].astype(int)  
        data.sort_values('Iteration', inplace=True)  

        for rod, color in zip(rod_sizes, colors):
            vpr_column = f"VPR {rod}"
            resol_column = f"Resol {rod}"
            
            # Cycle through iterations / (50/subsets)
            # print(subsets)
            # print(len(data['Iteration']))
            # print(step)
            step = int(50 / subsets)
            x = data['Iteration'][step-1::step]/len(data['Iteration'])*10
            ax_vpr.plot(x, data[vpr_column][step-1::step], marker=marker, color=color, label=f"{label} Resol {rod}")
            ax_resol.plot(x, data[resol_column][step-1::step], marker=marker, color=color, label=f"{label} Resol {rod}")
            
        legend_handles_files.append(Line2D([0], [0], marker=marker, color='black', markerfacecolor='black', markersize=7, label=label))

    # Create color legend handles only once, outside the main loop
    for rod, color in zip(rod_sizes, colors):
        legend_handles_colors.append(mpatches.Patch(color=color, label=f"Rod size {rod}"))

    # Set titles, labels, and legends for VPR subplot
    ax_vpr.set_title('VPR by Rod Size Across Iterations/Subsets', fontsize=fontsize)
    ax_vpr.set_xlabel('Iteration/Subsets', fontsize=fontsize)
    ax_vpr.set_ylabel('VPR', fontsize=fontsize)
    ax_vpr.tick_params(axis='x', labelsize=fontsize-2)
    ax_vpr.tick_params(axis='y', labelsize=fontsize-2)
    ax_vpr.grid(True)

    # Set titles, labels, and legends for Resolution subplot
    ax_resol.set_title('Resolvability by Rod Size Across Iterations/Subsets', fontsize=fontsize)
    ax_resol.set_xlabel('Iteration/Subsets', fontsize=fontsize)
    ax_resol.set_ylabel('Resolution', fontsize=fontsize)
    ax_resol.tick_params(axis='x', labelsize=fontsize-2)
    ax_resol.tick_params(axis='y', labelsize=fontsize-2)
    ax_resol.grid(True)

    # Add legend for rod sizes and markers
    rod_legend = ax_vpr.legend(handles=legend_handles_colors, loc=legend_loc[0])
    marker_legend = ax_resol.legend(handles=legend_handles_files, labels=labels, loc=legend_loc[1])
    # Add both legends to the plot
    ax_vpr.add_artist(rod_legend)
    ax_resol.add_artist(marker_legend)

    # Adjust layout for clarity
    plt.tight_layout()
    if save_bool: plt.savefig(save_path)
    plt.show()

# ------------------------------------------------------------
# IDEAS:
# position of the legend as option adding the default as it is now
# ------------------------------------------------------------

# Function to plot multiple CSV files on the same plot for VPR and Resol
def plot_vpr_resol(files, labels, rod_sizes, figsize=(14, 8), fontsize=12, legend_loc = ["lower left", "upper left"], save_bool=False, save_path='output'):

    plot = plt.figure(figsize=figsize)
    # plt.style.use('seaborn-v0_8-talk')

    # Defining a list of markers to cycle through
    markers = cycle(['o', 'x', 'D', 's', 'p', '*', '^', 'h'])

    # Defining a list of colors to cycle through
    colors = cycle(['black', 'red', 'green', 'blue', 'orange', 'magenta'][:len(rod_sizes)])
    # Creating subplots for VPR and Resolution
    ax_vpr = plt.subplot(2, 1, 1)
    ax_resol = plt.subplot(2, 1, 2)

    legend_handles_colors = []  # Store legend handles for colors
    legend_handles_files = []  # Store legend handles for files

    for file_path, label, marker in zip(files, labels, markers):
        # Load data
        data = pd.read_csv(file_path)
        data['Iteration'] = data['Iteration'].astype(int)  # Convert Iteration column to int
        data.sort_values('Iteration', inplace=True)  # Sort data by Iteration

        # Plotting VPR for each rod size
        for rod, color in zip(rod_sizes, colors):
            vpr_column = f"VPR {rod}"
            resol_column = f"Resol {rod}"
            
            ax_vpr.plot(data['Iteration'], data[vpr_column], marker=marker, markersize=6, color=color, label=f"{label} VPR {rod}")
            
            # Plot Resolution on the second subplot
            ax_resol.plot(data['Iteration'], data[resol_column], marker=marker, markersize=6, color=color, label=f"{label} Resol {rod}")
            
        legend_handles_files.append(Line2D([0], [0], marker=marker, color='black', markerfacecolor='black', markersize=7, label=label))

    # Create color legend handles only once, outside the main loop
    for rod, color in zip(rod_sizes, colors):
        legend_handles_colors.append(mpatches.Patch(color=color, label=f"Rod size {rod}"))

    # Set titles, labels, and legends for VPR subplot
    ax_vpr.set_title('VPR by Rod Size Across Iterations', fontsize=fontsize)
    ax_vpr.set_xlabel('Iteration', fontsize=fontsize)
    ax_vpr.set_ylabel('VPR', fontsize=fontsize)
    ax_vpr.tick_params(axis='x', labelsize=fontsize-2)
    ax_vpr.tick_params(axis='y', labelsize=fontsize-2)
    ax_vpr.grid(True)
    # plt.xticks(fontsize=fontsize)
    # plt.yticks(fontsize=fontsize)
    # Set titles, labels, and legends for Resolution subplot
    ax_resol.set_title('Resolvability by Rod Size Across Iterations', fontsize=fontsize)
    ax_resol.set_xlabel('Iteration', fontsize=fontsize)
    ax_resol.set_ylabel('Resolvability', fontsize=fontsize)
    ax_resol.tick_params(axis='x', labelsize=fontsize-2)
    ax_resol.tick_params(axis='y', labelsize=fontsize-2)
    ax_resol.grid(True)

    # Add legend for rod sizes and markers
    rod_legend = ax_vpr.legend(handles=legend_handles_colors, loc=legend_loc[0], fontsize=fontsize)
    marker_legend = ax_resol.legend(handles=legend_handles_files, labels=labels, loc=legend_loc[1], fontsize=fontsize)
    # Add both legends to the plot
    ax_vpr.add_artist(rod_legend)
    ax_resol.add_artist(marker_legend)

    # Adjust layout for clarity
    plt.tight_layout()
    if save_bool: plt.savefig(save_path)
    plt.show()

    return plot

# List of file paths
# file_paths = ['example/output_4test_rs.csv', 'example/output_5test_rs.csv']
# labels = ['Test 4', 'Test 5']  # Labels corresponding to each file
# rod_sizes = [1.1, 1.0, 0.95, 0.9, 0.85, 0.8]  # Rod sizes
# rod_sizes_1 = [1.1, 1.0, 0.95]  # Less rod sizes
# rod_sizes_2 = [0.9, 0.85, 0.8]  # Less rod sizes
# rod_sizes_3 = [0.8]  # Less rod sizes

# # Call the function with the list of file paths, labels, and rod sizes
# # plot_vpr_resol(file_paths, labels, rod_sizes_2)


# path_1 = '../../Results/rsmall_z2s/vs015_MLEM_projs/'
# path_2 = '../../Results/rsmall_z2s/vs015_opts_dri/D95/'
# path_3 = '../../Results/rsmall_z2s/vs015_opts_dri/'
# # path_4 = '../HiRezBrainPET_git/BPET_castor_results/Results/rsmall_z2s/vs015_MLEM_dri_g12_3s/'

# folder_1 = ['jos','dri','inc','sid']
# folder_2 = ['b01', 'b03', 'b05', 'b08', 'b1', 'b10', 'b100']
# folder_3 = ['MLEM', 'OSL_F/b1', 'OSL_P/b1', 'D95/b1']
# name_1 = 'vs015_MLEM_projs_parabola_z2030'
# name_2 = 'vs015_opts_dri_D95_betas_parabola_z2030'
# name_3 = 'vs015_opts_dri_parabola_z2030'

# output_paths = [[], [], []]
# paths = [path_1, path_2, path_3]
# folders = [folder_1, folder_2, folder_3]
# names = [name_1, name_2, name_3]

# for i, (path, folder, name) in enumerate(zip(paths, folders, names)):
#     for f in folder:
#         output_paths[i].append(path + f + '/output_' + name + '.csv')

# plot_vpr_resol(output_paths[1], folder_1, rod_sizes)
# plot_vpr_resol(output_paths[1], folder_1, rod_sizes_1)
# plot_vpr_resol(output_paths[1], folder_1, rod_sizes_2)
# plot_vpr_resol(output_paths[1], folder_1, rod_sizes_3)
# plot_vpr_resol(output_paths[2], folder_3, rod_sizes)
# plot_vpr_resol(output_paths[1], folder_2, rod_sizes_3)



path_4 = '../HiRezBrainPET_git/BPET_castor_results/Results/rsmall_z2s/vs015_MLEM_jos_ittest/'
folder_4 = ['1test', '2test', '3test', '4test', '5test', '6test', '7test', '8test']
labels_4 = ['10:50', '20:25', '25:20', '50:10', '100:5', '250:2', '500:1', '2:50,40:10']
name_4 = 'vs015_MLEM_jos_ittest'
subsets = [50,25,20,10,5,2,1,10]

# folder_4 = ['2test', '4test']
# labels_4 = ['20:25', '50:10']
# subsets = [25,10]

# output_paths_its = [path_4 + f + '/output_' + name_4 + '.csv' for f in folder_4]
# plot_vpr_resol_it_subsets(output_paths_its, labels_4, subsets, rod_sizes_3)
  